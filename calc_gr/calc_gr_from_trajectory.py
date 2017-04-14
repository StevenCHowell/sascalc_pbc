import logging
import os
import sys
import time
import numpy as np
from scipy.spatial.distance import pdist
import sasmol.sasmol as sasmol

sys.path.append('./')
try:
    import gr as fortran_gr
except ImportError:
    pass

FORMAT = "%(asctime)-15s: %(message)s"
logging.basicConfig(format=FORMAT, level=logging.DEBUG)


def slow_update_gr(coor, box_length, gr, dr):

    n_atoms = len(coor)

    for i in xrange(n_atoms - 1):
        for j in xrange(i + 1, n_atoms):

            dx = coor[i, 0] - coor[j, 0]
            dy = coor[i, 1] - coor[j, 1]
            dz = coor[i, 2] - coor[j, 2]

            dx -= box_length * ((dx / box_length).round())
            dy -= box_length * ((dy / box_length).round())
            dz -= box_length * ((dz / box_length).round())

            r = np.sqrt((dx * dx) + (dy * dy) + (dz * dz))

            if (r < box_length / 2.0):
                ig = int(r / dr)  # round down
                gr[ig] += 2


def calc_gr(coor, box_length, gr, dr):

    dx = pdist(coor[:, :1])
    dy = pdist(coor[:, 1:2])
    dz = pdist(coor[:, 2:])

    dx -= box_length * ((dx / box_length).round())
    dy -= box_length * ((dy / box_length).round())
    dz -= box_length * ((dz / box_length).round())

    r = np.sqrt((dx * dx) + (dy * dy) + (dz * dz))

    r_i = (r[r < box_length / 2.0] / dr).astype(int)  # round down
    r_i_unique = np.unique(r_i, return_counts=True)
    gr[r_i_unique[0]] += 2 * r_i_unique[1]


def main(pdb_fname, run_log_fname, stride=1, sigma=1, dcd_fname=None,
         n_skip=0, n_bins=333, gr_fname=''):

    assert op.exists(pdb_fname), 'No such pdb file: {}'.format(pdb_fname)
    box_mol = sasmol.SasMol(0)
    box_mol.read_pdb(pdb_fname)

    if dcd_fname:
        assert op.exists(dcd_fname), 'No such dcd file: {}'.format(dcd_fname)
        dcd_file = box_mol.open_dcd_read(dcd_fname)
        n_frames = dcd_file[2]
    else:
        n_frames = 1

    run_log = np.loadtxt(run_log_fname)

    if len(run_log) != n_frames:
        if len(run_log) == n_frames + 1:
            run_log = run_log[:-1]
            logging.warning('dcd file  had one more frame than the log file, '
                         'discarding last line\ndcd_fname:\t{}\n'
                         'run_log_fname:\t{}'.format(dcd_fname, run_log_fname))
        else:
            logging.error('mismatch between dcd and log file \n'
                          'dcd_fname:\t{}\nrun_log_fname:\t{}'.format(
                              dcd_fname, run_log_fname))

    n_atoms = box_mol.natoms()
    n_gr = n_frames - n_skip  # number of frames, or g(r) curves, to averages

    box_length = run_log[:, 1] * sigma
    print('box_length: (min, max) = ({}, {})'.format(box_length.min(),
                                                     box_length.max()))

    gr_all = np.zeros((n_gr, n_bins))  # one g(r) for ecah dcd frame
    gr = np.zeros((n_bins, 2))

    # using the same r_grid for each frame
    dr = box_length[n_skip:].max() / (2.0 * n_bins)  # delg in F&S
    print(dr)
    bin_index = np.arange(n_bins) + 1
    rho = (n_atoms / (box_length ** 3.0)).reshape(-1, 1)  # frame specific density

    r = (bin_index - 0.5) * dr
    bin_volume = 4.0 / 3.0 * np.pi * ((r + dr / 2) ** 3 - (r - dr / 2) ** 3)
    n_ideal = bin_volume * rho  # expected n for ideal gas

    if n_skip:
        for i in xrange(n_skip):
            # read and throw away these coordinates
            box_mol.read_dcd_step(dcd_file, i)

    tic = time.time()
    for i in xrange(n_skip, n_frames):
        sys.stdout.flush()

        try:
            box_mol.read_dcd_step(dcd_file, i)  #, no_print=True)
        except NameError:
            print('calculating g(r) for {}'.format(pdb_fname))

        coor = box_mol.coor()[0] * sigma

        # slow_update_gr(coor, box_length[i], gr_all[i-n_skip], dr)
        # calc_gr(coor, box_length[i], gr_all[i-n_skip], dr)
        fortran_gr.calc_gr(coor, box_length[i], gr_all[i-n_skip], dr)

        gr_all[i-n_skip] /= n_ideal[i]  # normalize expected n for ideal gas

    toc = time.time() - tic
    print('\nrun time: {} seconds'.format(toc))
    box_mol.close_dcd_read(dcd_file[0])
    gr[:, 0] = r
    gr[:, 1] = np.mean(gr_all, axis=0) / n_atoms  # normalize by frames and atoms

    if not gr_fname:
        gr_fname = 'gr_{}.dat'.format(pdb_fname[:-4])

    gr = gr[gr[:, 0] < box_length.min()/2]
    np.savetxt(gr_fname, gr, fmt='%.14f')

    return gr


def plot_gr(gr, stride=1, show=False, plot_fname=''):
    import matplotlib.pyplot as plt

    gr_dat = np.loadtxt('argon_85K_gr.dat')

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.set_ylabel('g(r)')
    ax1.set_xlabel('r')
    scale_factor = 1.0
    ax1.plot(gr[:, 0], scale_factor * gr[:, 1], color='red', lw=2)
    ax1.plot(gr_dat[:, 0], gr_dat[:, 1], color='blue', lw=2)

    if not plot_fname:
        plot_fname = 'steve_gr.png'

    plt.savefig(plot_fname)

    if show:
        plt.show()
    else:
        plt.close('all')


if __name__ == '__main__':
    import os.path as op
    sigma = 3.405

    test = False

    if test:
        run_path = './run2_output'
        pdb_fname = 'run2.pdb'
        dcd_fname = 'run2_last100.dcd'
        xst_fname = 'box_length_last100.txt'

        pdb_fname = op.join(run_path, pdb_fname)
        dcd_fname = op.join(run_path, dcd_fname)
        xst_fname = op.join(run_path, xst_fname)

        n_skip = 0

        gr_fname = 'gr_fortran.dat'; plot_fname = 'gr_fortran.png'
        # gr_fname = 'gr_python.dat'; plot_fname = 'gr_python.png'

    else:
        run_path = '../../simulations/lj_sphere_monomer/runs/p_0p14/run2_output'
        pdb_fname = 'run2.pdb'
        dcd_fname = 'run2.dcd'
        xst_fname = 'box_length.txt'

        pdb_fname = op.join(run_path, pdb_fname)
        dcd_fname = op.join(run_path, dcd_fname)
        xst_fname = op.join(run_path, xst_fname)

        n_skip = 1000
        # n_skip = 25000
        # gr_fname = 'gr_1000_25001.dat'; plot_fname = 'gr_1000_25001.png'
        gr_fname = 'gr_last_frame.dat'; plot_fname = 'gr_last_frame.png'

    gr = main(pdb_fname, xst_fname, sigma=sigma, dcd_fname=dcd_fname,
              n_skip=n_skip, gr_fname=gr_fname)

    plot_gr(gr, plot_fname=plot_fname)
