'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import sys,os,math,locale,string
import numpy,scipy
import random
from sassie.sasmol import sasmol as sasmol
import sassie.interface.input_filter as input_filter

sys.path.append('./')
import mmc_sq_filter	### OPEN -- this is imported but not called yet
import nenergy

#
#	MMC_SQ
#
#	1/26/2010	--	initial coding			:	jc
#	6/10/2010	--	adapted to npt			:	jc
#	10/24/2010	--	at-npt				:	jc
#	07/10/2012	--	re-factoring began		:	jc
#
#	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	MMC_SQ is the module that performs NPT molecular Monte Carlo simulations 
	of ensembles of bio-macromolecules.  There are options to allow the 
	molecules to have internal degrees of freedom (flexibility between
	segments) and for the simulations to adapt the pressure to match
	experimental scattering data.

	The program reads in a single starting structure and creates a box
	with a user defined number of copies of the initial molecule at a given
	starting density on an fcc lattice.  Options to restart from pre-existing
	configuations are provided.

	The user will also supply a temperature, pressure, and a set of 
	Monte Carlo parameters (number of steps, translational move size,
	box-volume move size.  The user will specify the number of coordinates
	to save.

	For internal moves the user will define flexible residues between
	rigid domains and provide moveset parameters (delta theta, internal 
	alignment basis).  Note that the initial coordinates are linearized
	to determine the surface residues that will be considered in the 
	calculation of non-bonding forces.
	
	INPUT/OUTPUT: (see get_parameters() method below for full details)	

	In general:

	INPUT: 
	
		PDB file with initial single molecular structure
		MC parameters (defined above)
		IDP parameters (defined above)
		Scattering data

	OUTPUT:
	
		PDB of final coordinates of all molecules
		DCD binary file with ensemble coordinates
		Trajectory data (synthetic I(Q) -vs- experimental I(Q)

'''

def stop_here():

	print '>>>> STOPPING HERE <<<<'
	sys.exit()

def get_parameters():

	'''
		This method will be used in the GUI and left here for stand-alone access

	'''

	svariables = {}


	svariables['runname'] = ('lrun_0','string') 				# prefix to name output files
	
	svariables['initial_pdb'] = ('argon.pdb','string')			# PDB structure file for initial single molecular structure
	#svariables['initial_pdb'] = ('lysozyme_md.pdb','string')			# PDB structure file for initial single molecular structure

	svariables['scattering_data_file'] = ('scat.dat','string') 		# experimental scattering data file

	svariables['start_type']= ('1','int')					# 1--> make fcc coords, 2--> restart from molecular coords
	#svariables['start_type']= ('2','int')					# 1--> make fcc coords, 2--> restart from molecular coords
	
	svariables['number_of_molecules'] = ('2048','int')			# number of replicated molecules in system

	svariables['number_of_steps'] = ('1000','int')			# number of MC steps for each cycle
	svariables['number_of_cycles'] = ('0','int')				# number of independent cycles to converge to scattering data (if == 0; don't compare)

	svariables['dcd_save_frequency'] = ('1000','int')			# mc steps between binary coordinate saves
	
	svariables['delta_translation'] = ('0.02','float')			# translational move step size
	svariables['delta_rotation'] = ('0.0','float')				# molecular rotation move step size
	svariables['delta_volume'] = ('0.001','float')				# volume move step size (percent)

	#svariables['temperature'] = ('298.0','float')				# temperature (Kelvin)
	#svariables['temperature'] = ('1.5043','float')				# temperature (85 Kelvin in reduced units)
	svariables['temperature'] = ('0.72','float')				# temperature (85 Kelvin in reduced units)
	#svariables['pressure'] = ('1.0','float')				# pressure (bar)
	svariables['pressure'] = ('1.0','float')				# pressure (bar)
	svariables['initial_rho'] =  ('0.8442','float')				# initial packing fraction


	svariables['internal_dof'] = ('NO','string')				# allow internal dof ('YES' or 'NO')
	svariables['number_of_flexible_regions'] = ('2','int')			# number of flexible regions
	svariables['delta_theta'] = ('30.0, 30.0','float_array')		# maximum delta theta per region
	svariables['low_residues'] = ('123, 278','int_array')			# starting residue for each flexible region
	svariables['number_of_contiguous_residues'] = ('21, 5','int_array')	# number of contiguous residues for each region
	svariables['align_low'] = ('284','int')					# molecular alignment low residue
	svariables['align_high'] = ('350','int')				# molecular alignment high residue


	nsteps = locale.atoi(svariables['number_of_steps'][0])

	svariables['translation_ratio'] = (str(nsteps/100),'float')		# parameter to re-assess frequency of translation move step size
	svariables['rotation_ratio'] = (str(nsteps/100),'float')		# parameter to re-assess frequency of molecular rotation move step size 
	svariables['volume_ratio'] = (str(nsteps/100),'float')			# parameter to re-assess frequency of box-volume move step size 
	svariables['dtheta_ratio'] = (str(nsteps/100),'float')			# parameter to re-assess frequency of internal DOF move step size 

	return svariables


def unpack_variables(variables):

	runname					= variables['runname'][0]
	initial_pdb 				= variables['initial_pdb'][0]
	scattering_data_file 			= variables['scattering_data_file'][0]
	start_type				= variables['start_type'][0]
	number_of_molecules			= variables['number_of_molecules'][0]
	number_of_steps				= variables['number_of_steps'][0]
	number_of_cycles			= variables['number_of_cycles'][0]
	dcd_save_frequency			= variables['dcd_save_frequency'][0]
	delta_translation			= variables['delta_translation'][0]
	delta_rotation				= variables['delta_rotation'][0]
	delta_volume				= variables['delta_volume'][0]
	temperature				= variables['temperature'][0]
	pressure				= variables['pressure'][0]
	initial_rho				= variables['initial_rho'][0]
	internal_dof				= variables['internal_dof'][0]
	number_of_flexible_regions		= variables['number_of_flexible_regions'][0]
	delta_theta				= variables['delta_theta'][0]
	low_residues				= variables['low_residues'][0]
	number_of_contiguous_residues		= variables['number_of_contiguous_residues'][0]
	align_low				= variables['align_low'][0]
	align_high				= variables['align_high'][0]
	translation_ratio			= variables['translation_ratio'][0]
	rotation_ratio				= variables['rotation_ratio'][0]
	volume_ratio				= variables['volume_ratio'][0]
	dtheta_ratio				= variables['dtheta_ratio'][0]


	return runname,initial_pdb,scattering_data_file,start_type,number_of_molecules,number_of_steps,number_of_cycles,dcd_save_frequency,delta_translation,delta_rotation,delta_volume,temperature,pressure,initial_rho,internal_dof,number_of_flexible_regions,delta_theta,low_residues,number_of_contiguous_residues,align_low,align_high,translation_ratio,rotation_ratio,volume_ratio,dtheta_ratio
	

def overlap(allcoor,trialcoor,diameter):

        check=0
        tx=trialcoor[0] ; ty=trialcoor[1] ; tz=trialcoor[2]
        for coor in allcoor:
                x=coor[0] ; y=coor[1] ; z=coor[2]
                dist=math.sqrt((x-tx)*(x-tx)+(y-ty)*(y-ty)+(z-tz)*(z-tz))
                if dist < diameter:
                        check=1
        return check

def makefcc(diameter):

        pi=numpy.pi
        eta=0.74048

        matms=3
        matms=8
        natms=4*matms*matms*matms

        volp=natms*4.0*pi*((diameter/2.0)**3.0)/3.0
        volbox=volp/eta
        boxl=volbox**(1.0/3.0)

        xn = [-0.25,0.25,-0.25,0.25]
        yn = [-0.25,0.25,0.25,-0.25]
        zn = [-0.25,-0.25,0.25,0.25]

        print 'boxl = ',boxl

        sid=boxl/matms

        kk=0
        x=numpy.zeros(natms,float)
        y=numpy.zeros(natms,float)
        z=numpy.zeros(natms,float)
        for n in xrange(4):
                for i in xrange(matms):
                        for j in xrange(matms):
                                for k in xrange(matms):
                                        x[kk]=(float(k)-0.50+xn[n])*sid-0.5*boxl
                                        y[kk]=(float(j)-0.50+yn[n])*sid-0.5*boxl
                                        z[kk]=(float(i)-0.50+zn[n])*sid-0.5*boxl
                                        kk=kk+1

	coor=numpy.zeros((natms,3),numpy.float)
        coor[:,0]=x ; coor[:,1]=y ; coor[:,2]=z

        return coor,boxl

def setup(nmolecules,rho_start):
	
	radius = 0.5521 #  argon VDW radius / sigma argon

	molvolume = ((nmolecules*4.0*numpy.pi*(radius**3.0))/3.0)
	boxvolume = nmolecules/rho_start
	packing_fraction = molvolume/boxvolume
	print 'packing fraction = ',packing_fraction
	print 'boxvolume = ',boxvolume
	
	boxl=boxvolume**(1/3.)
	print 'boxl = ',boxl	
	print 'radius = ',radius
	diameter=2.0*radius

	allcoor=[] 
	for i in range(nmolecules):
		search=True		
		while(search):
			trialcoor=boxl*scipy.random.random(3)
			check=overlap(allcoor,trialcoor,diameter)
			if(check==0):
				allcoor.append(trialcoor)
				search=False

	#coor=numpy.array(allcoor)-(boxl/2.0) #	-L/2 to L/2 
	coor=numpy.array(allcoor)
	
	return coor,boxl,radius

def sumup(coor,natoms,rcut,rmin,boxl): #,v12,v16,w12,w6):

	xcoor=coor[:,0] ; ycoor=coor[:,1] ; zcoor=coor[:,2]
	rcutsq = rcut*rcut
	rminsq = rmin*rmin
	boxlinv = 1.0/boxl
	
	v12 = 0.0 ; v6 = 0.0 ; w12 = 0.0 ; w6 = 0.0

	for i in xrange(natoms-1):
		rxi=xcoor[i] ; ryi=ycoor[i] ; rzi=zcoor[i]
		for j in xrange(i+1,natoms):
			rxj=xcoor[j] ; ryj=ycoor[j] ; rzj=zcoor[j]
			rxij = rxi - rxj	
			ryij = ryi - ryj	
			rzij = rzi - rzj
	
			rxij = rxij-boxl*((rxij*boxlinv).round())
			ryij = ryij-boxl*((ryij*boxlinv).round())
			rzij = rzij-boxl*((rzij*boxlinv).round())

			rijsq = rxij*rxij + ryij*ryij + rzij*rzij
					
			if(rijsq < rcutsq):
				sr2 = 1.0/rijsq
				sr6 = sr2*sr2*sr2
				vij12 = sr6*sr6
				vij6 = -sr6
				v12 = v12 + vij12
				v6 = v6 + vij6
				w12 = w12 + vij12
				w6 = w6 + vij6*0.5	

	v12 = 4.0*v12
	v6 = 4.0*v6
	w12 = 48.0*w12/3.0
	w6 = 48.0*w6/3.0

	return v12,v6,w12,w6

def scale_coords_and_build_box(single_molecule,start_type):

	ljradius = 0.5521 #  argon VDW radius / sigma argon
	sigmaAr=3.405
	diameter_of_Ar = 3.76 # van der Waals diameter
	ljradius = (diameter_of_Ar/2.0) / sigmaAr


	single_molecule.center(0)
	mm1=single_molecule.calcminmax()

	dx=math.fabs(mm1[0][0]-mm1[1][0])
	dy=math.fabs(mm1[0][1]-mm1[1][1])
	dz=math.fabs(mm1[0][2]-mm1[1][2])

	if(dx == dy == dz == 0.0):
		dx = diameter_of_Ar
		dy = diameter_of_Ar
		dz = diameter_of_Ar

	maxdiamolecule = max(dx,dy,dz) ; print 'maximum molecular diameter = ',maxdiamolecule
	mindiamolecule = min(dx,dy,dz) ; print 'minimum molecular diameter = ',mindiamolecule

	print 'max diameter = ',maxdiamolecule
	print 'min diameter = ',mindiamolecule

	scale_coords = sigmaAr*maxdiamolecule/(ljradius*2.0*sigmaAr)

	print 'scale_coords = ',scale_coords		

	if(start_type==1):
	
		coor,boxl=makefcc(2.0*ljradius)
			
	elif(start_type==2):

		coor,boxl=readxyz(coorfilename)

	else:
		print 'incorrect start up type indicated ',start_type
		print '# 1--> make fcc coords, 2--> restart from PDB coords'
		print 'STOPPING NOW\nSTOPPING NOW\nSTOPPING NOW\n\n',3/0

	flag='new' ; filename='traj0.pdb'
#	write_frame(tempmol,coor,scale_lj_to_protein,flag,filename)

	return coor,ljradius,boxl,scale_coords


def duplicate_molecule_and_assign_coor(all_molecules,single_molecule,coor):

	initial_natoms = single_molecule.natoms()

	number_of_duplicates = len(coor)/initial_natoms
	print ' >> found ',number_of_duplicates,' molecular copies'
	frame = 0

	error = single_molecule.duplicate_molecule(all_molecules,number_of_duplicates,frame,coor)	

	all_molecules.coor()[0,:,:] = coor

	number_atoms = all_molecules.natoms()

	print all_molecules.index()
	
	print ' >> all natoms = ',number_atoms


	original_segname = all_molecules.segname() 

	segname = [] ;
	for i in xrange(2048):
		segname.append("ARG")

	all_molecules.setSegname(segname)

	all_molecules.write_pdb("dum_test.pdb",frame,'w')

	all_molecules.setSegname(original_segname)
	
	all_molecules.write_pdb("all_test.pdb",frame,'w')

	return


def read_restart_pdb(initial_pdb,all_molecules):

	scale_coords = 1.0
	
	ljradius = 0.5521 #  argon VDW radius / sigma argon
	sigmaAr=3.405
	diameter_of_Ar = 3.76 # van der Waals diameter
	ljradius = (diameter_of_Ar/2.0) / sigmaAr

	infile = open(initial_pdb,'r').readlines()
	lin = string.split(infile[0])
	boxl = locale.atof(lin[3])

	print '\n > boxl = ',boxl,'\n\n'
	
	all_molecules.read_pdb(initial_pdb)

	coor = all_molecules.coor()[0,:,:]

	return coor,ljradius,boxl,scale_coords



def nptmc(number_of_molecules,temperature,goalpressure,initial_rho,number_of_steps,nbins,delta_translation,delta_volume,runname,initial_pdb,translation_ratio,volume_ratio,specificmolvol,start_type,dcd_save_frequency):

	
	frame = 0 ; overall_frame = 0
	
	pi=numpy.pi

	all_molecules = sasmol.SasMol(1)
	
	#start_type: 1--> make fcc coords, 2--> restart from molecular coords(PDB)
	
	if(start_type==1):
	
		single_molecule = sasmol.SasMol(0)
		single_molecule.read_pdb(initial_pdb)
		mass = single_molecule.calcmass(); print '> single molecule mass = ',mass

		coor,radius,boxl,scale_coords = scale_coords_and_build_box(single_molecule,start_type)

		duplicate_molecule_and_assign_coor(all_molecules,single_molecule,coor)

		total_mass = all_molecules.calcmass() ; print ' > total mass = ',total_mass

	elif(start_type==2):

		coor,radius,boxl,scale_coords = read_restart_pdb(initial_pdb,all_molecules)
		total_mass = all_molecules.calcmass() ; print ' > total mass = ',total_mass
		all_molecules.write_pdb(runname+'_initial_coords.pdb',frame,'w')
	
	else:
		print 'incorrect start up type indicated ',start_type
		print '# 1--> make fcc coords, 2--> restart from pdb file'
		print 'STOPPING NOW\nSTOPPING NOW\nSTOPPING NOW\n\n',3/0

	deltag = boxl/(2.0*nbins)
	gr = numpy.zeros(nbins,numpy.float)

	dcdoutfile = all_molecules.open_dcd_write(runname+'.dcd')

	boxfile=open(runname+'_box.dat','w')
	pressfile=open(runname+'_press.dat','w')
	densityfile=open(runname+'_density.dat','w')

	natoms=number_of_molecules 	# for now

	volume = boxl**3.0
	boxlinv = 1.0/boxl
	density = natoms/volume

	rcut = 0.5*boxl

	dboxmx = boxl/40.0
	drmax = delta_translation
	rmin = 0.70
	beta = 1.0/temperature	# i.e. kb=1.0
	
	sr3 = (1.0/rcut)**3.0
	sr9 = sr3**3.0
	vlrc12 = 8.0*pi*density*natoms*sr9/9.0
	vlrc6 = -8.0*pi*density*natoms*sr3/3.0
	wlrc12 = 4.0*vlrc12
	wlrc6 =	2.0*vlrc6
	vlrc = vlrc12 + vlrc6
	wlrc = wlrc12 + wlrc6

	acm = 0.0 ; acatma = 0.0 ; acboxa = 0.0 ; acv = 0.0 ; acp = 0.0	; acd = 0.0
	acvsq = 0.0 ; acpsq = 0.0 ; acdsq = 0.0 ; flv = 0.0 ; flp = 0.0	; fld = 0.0

	v12,v6,w12,w6 = sumup(coor,natoms,rcut,rmin,boxl) #,v12,v6,w12,w6)
 
	vs = ( v12 + v6 + vlrc ) / natoms
	ws = ( w12 + w6 + wlrc ) / natoms
	ps = density * temperature + ( w12 + w6 + wlrc ) / volume

	v12 = v12 + vlrc12
	v6 = v6 + vlrc6
	w12 = w12 + wlrc12
	w6 = w6 + wlrc6

	print 'sr3 = ',sr3
	print 'sr9 = ',sr9
	print 'vlrc12 = ',vlrc12
	print 'vlrc6 = ',vlrc6
	print 'wlrc12 = ',wlrc12
	print 'wlrc6 = ',wlrc6
	print 'vlrc = ',vlrc
	print 'wlrc = ',wlrc
	print ''

	print 'v12 = ',v12
	print 'v6 = ',v6
	print 'w12 = ',w12
	print 'w6 = ',w6
	
	print 'density = ',density
	print 'initial energy/atom = ',vs
	print 'initial w/atom = ',ws
	print 'initial pressure = ',ps

	molvolume = ((natoms*4.0*numpy.pi*(radius**3.0))/3.0)

	m=0 ; v=0 ; tm=0 ; tv=0 
	sumeta=0.0

	fr=1

	print '\nSTEP\t\tVN\t%ACC_MOV\t%ACC_VOL\t ETA\t\t<ETA>\t\tDENSITY\t\t<DENSITY>\tPRESSURE\t<PRESSURE>\tBOXL\n'
		#print '%i\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%10.6f\t%8.6f\n' % (i+1,acatma*100.0/m,acboxa*100.0/v,density,acd/(m+v),pressure,acp/(m+v),boxl)




	for step in xrange(number_of_steps):
		if(step == 0): 
			frames_since_last_dcd_save = dcd_save_frequency
		else:
			frames_since_last_dcd_save += 1	
		
		for i in xrange(natoms):
			m=m+1 ; tm=tm+1		
			rxiold = coor[i,0]
			ryiold = coor[i,1]
			rziold = coor[i,2]
			thisi=i+1
			
			tgr = numpy.zeros(nbins,numpy.float)

			tgr,v12old,v6old,w12old,w6old=nenergy.lje(coor,rxiold,ryiold,rziold,thisi,boxl,deltag,nbins)

			rxinew = rxiold + ( 2.0 * random.random() - 1.0 ) * drmax
			ryinew = ryiold + ( 2.0 * random.random() - 1.0 ) * drmax
			rzinew = rziold + ( 2.0 * random.random() - 1.0 ) * drmax

			rxinew = rxinew-boxl*((rxinew*boxlinv).round())
			ryinew = ryinew-boxl*((ryinew*boxlinv).round())
			rzinew = rzinew-boxl*((rzinew*boxlinv).round())

			tgr = numpy.zeros(nbins,numpy.float)
			
			tgr,v12new,v6new,w12new,w6new=nenergy.lje(coor,rxinew,ryinew,rzinew,thisi,boxl,deltag,nbins)

			delv12 = v12new - v12old
			delv6 = v6new - v6old
			delw12 = w12new - w12old
			delw6 = w6new - w6old
			
			deltv = delv12 + delv6
			deltvb = beta * deltv

			if(deltvb < 75.0):
				if(deltv <= 0.0):
					v12=v12+delv12
					v6 =v6+delv6
					w12=w12+delw12
					w6 =w6+delw6
					coor[i,0]=rxinew
					coor[i,1]=ryinew
					coor[i,2]=rzinew
					acatma=acatma+1.0
					gr = gr + tgr
				elif(math.exp(-deltvb) > random.random()):
					v12=v12+delv12
					v6 =v6+delv6
					w12=w12+delw12
					w6 =w6+delw6
					coor[i,0]=rxinew
					coor[i,1]=ryinew
					coor[i,2]=rzinew
					acatma=acatma+1.0
					gr = gr + tgr

	
			#vlrc12 = 8.0*pi*density*natoms*sr9/9.0
			#vlrc6 = -8.0*pi*density*natoms*sr3/3.0
			#wlrc12 = 4.0*vlrc12
			#wlrc6 =	2.0*vlrc6

			vn = ( v12 + v6 ) / natoms
			#vn = ( v12 + vlrc12 + v6 + vlrc6 ) / natoms
			pressure = density * temperature + ( w12 + w6 ) / volume
			#pressure = density * temperature + ( w12 + wlrc12 + w6 + wlrc6 ) / volume

			acm = acm + 1.0
			acv = acv + vn
			acp = acp + pressure
			acd = acd + density
	
			acvsq = acvsq + vn**2.0
			acpsq = acpsq + pressure**2.0
			acdsq = acdsq + density**2.0
			
		#END of for i in xrange(natoms)
	
		v=v+1 ; tv=tv+1
		#boxlnew = -1.0
		#while(boxlnew < 0):
		#	boxlnew = boxl + ( 2.0 * random.random() - 1.0 ) * dboxmx
		
		boxlnew = boxl + ( 2.0 * random.random() - 1.0 ) * dboxmx
		ratbox = boxl / boxlnew
		rrbox = 1.0 / ratbox
		rcutn = rcut * rrbox

		rat6 = ratbox ** 6.0
		rat12 = rat6 * rat6

		v12new = v12 * rat12
		v6new = v6 * rat6
		w12new = w12 * rat12
		w6new = w6 * rat6

		deltv = v12new + v6new - v12 - v6
		dpv = goalpressure * (boxlnew**3.0 - volume)
		dvol= 3.0 * temperature * natoms * math.log(ratbox)
		delthb = beta * ( deltv + dpv + dvol )

			
		if(delthb < 75.0):
			if(delthb <= 0.0):
				v12 = v12new 
				v6  = v6new 
				w12 = w12new 
				w6  = w6new 

				coor=coor*rrbox

				boxl=boxlnew
				rcut=rcutn
				acboxa = acboxa + 1.0
			elif(math.exp(-delthb) > random.random()):
				v12 = v12new 
				v6  = v6new 
				w12 = w12new 
				w6  = w6new

				coor=coor*rrbox

				boxl=boxlnew
				rcut=rcutn
				acboxa = acboxa + 1.0
			

		boxlinv = 1.0/boxl
		volume = boxl**3.0
		density=natoms/volume

	
		#vlrc12 = 8.0*pi*density*natoms*sr9/9.0
		#vlrc6 = -8.0*pi*density*natoms*sr3/3.0
		#wlrc12 = 4.0*vlrc12
		#wlrc6 =	2.0*vlrc6
		vn = ( v12 + v6 ) /natoms
		#vn = ( v12 +vlrc12 + v6 + vlrc6) /natoms
		pressure = density * temperature + ( w12 + w6 ) / volume
		#pressure = density * temperature + ( w12 + wlrc12 + w6 + wlrc6) / volume

		eta=molvolume/boxl**3.0

		acm = acm + 1.0
		acv = acv + vn
		acp = acp + pressure
		acd = acd + density

		sumeta = sumeta + eta
	
		acvsq = acvsq + vn**2.0
		acpsq = acpsq + pressure**2.0
		acdsq = acdsq + density**2.0
	
		all_molecules.coor()[0,:,:] = coor
	
		if(frames_since_last_dcd_save == dcd_save_frequency):
	
			overall_frame += 1
			all_molecules.write_dcd_step(dcdoutfile,frame,overall_frame)
			frames_since_last_dcd_save = 0
		
		boxfile.write("%i\t%f\n" % (step,boxl))	
		pressfile.write("%i\t%f\t%f\n" % (step,pressure,acp/(tm+tv)))	
		densityfile.write("%i\t%f\t%f\n" % (step,density,acd/(tm+tv)))	
		boxfile.flush() ; pressfile.flush() ; densityfile.flush()
		if (step%(number_of_steps/100.))==0: 
			if (fr%20.)==0: 
				print '\nSTEP\t\tVN\t%ACC_MOV\t%ACC_VOL\t ETA\t\t<ETA>\t\tDENSITY\t\t<DENSITY>\tPRESSURE\t<PRESSURE>\tBOXL\n'
			print '%i\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%10.6f\t%8.6f\n' % (step,vn,acatma*100.0/m,acboxa*100.0/v,eta,sumeta/tv,density,acd/(tm+tv),pressure,acp/(tm+tv),boxl)

			fr+=1
		
		if (step%translation_ratio==0):
			ratio = acatma/(natoms*translation_ratio)
			if (ratio > 0.5):
				drmax=drmax*1.05
			else:
				drmax=drmax*0.95
			#print 'drmax = ',drmax
			acatma = 0.0
			m=0
		if (step%volume_ratio==0):
			bratio = acboxa/volume_ratio
			if (bratio > 0.5):
				dboxmx=dboxmx*1.05
			else:
				dboxmx=dboxmx*0.95
			#print 'dboxmx = ',dboxmx
			acboxa = 0.0
			v=0	

	print '\n\npercent mcmoves accepted ',acatma*100.0/tm
	print '\n\npercent mcvols accepted ',acboxa*100.0/tv

	print 'goal pressure = ',goalpressure

	#stop_here()	
	
	final_pdb = 'final_'+runname+'.pdb'
	outputpdbheader = open(final_pdb,'w')
	outputpdbheader.write('REMARK BOXL = '+str(boxl)+'\n')
	outputpdbheader.close()

	all_molecules.write_pdb(final_pdb,frame,'a')

	boxfile.close()
	pressfile.close()

        if (acboxa>0):
                #rho=1.0/((boxsum/accvol)**3.0)
		rho = acd/(tm+tv)
                print 'rho = ',rho
        else:
                rho=1.0/(boxl**3.0)

        goutfile=open(runname+'_gr.dat','w')
        for i in xrange(nbins):
                r=deltag*(i+0.5)
                vb=((i+1)**3-i**3)*deltag**3
                nid=(4./3)*numpy.pi*vb*rho
                #gr[i]=gr[i]/(accmove*natoms*nid)
                gr[i]=gr[i]/((acatma)*natoms*nid)
                #gr[i]=gr[i]/((tm)*natoms*nid)
                #goutfile.write('%i\t%f\t%f\n' % (i,r,gr[i]/nmolecules))
                goutfile.write('%i\t%f\t%f\n' % (i,r,gr[i]))

        goutfile.close()

	return


def print_input_parameters(variables):

	print '\n >> input variables :\n'
	for key in variables:
		print key,' = ',variables[key][0]

	print

	return

def main_program(variables):

	print '\n\n >> unpacking variables\n'

	runname,initial_pdb,scattering_data_file,start_type,number_of_molecules,number_of_steps,number_of_cycles,dcd_save_frequency,delta_translation,delta_rotation,delta_volume,temperature,pressure,initial_rho,internal_dof,number_of_flexible_regions,delta_theta,low_residues,number_of_contiguous_residues,align_low,align_high,translation_ratio,rotation_ratio,volume_ratio,dtheta_ratio = unpack_variables(variables)

	print ' >> calling nptmc'

	nbins = 1000
	specificmolvol = 9769

	nptmc(number_of_molecules,temperature,pressure,initial_rho,number_of_steps,nbins,delta_translation,delta_volume,runname,initial_pdb,translation_ratio,volume_ratio,specificmolvol,start_type,dcd_save_frequency)

	print '\n\n >> I am done'	

#	nptmc(nmolecules,temperature,goalpressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename,coorfilename,iratio,vratio,masterpdb,moleculepdb,specificmolvol,start_type)


if __name__=="__main__":


	svariables = get_parameters()

	error,variables=input_filter.type_check_and_convert(svariables)
	
	print_input_parameters(variables)

	main_program(variables)


	# for Argon, E/K_b = 119.8 K, sigma = 3.405E-10 m, M = 0.03994 kg/mol

	# length 	--> sigma
	# energy 	--> e
	# mass		--> m
	# time = sigma*sqrt(m/e)

	# potential energy == U* 	--> U/e
	# pressure == P* 		--> P * sigma^3 / e	
	# density == rho*		--> rho * sigma^3
	# temperature == T*		--> k_b * T / e

