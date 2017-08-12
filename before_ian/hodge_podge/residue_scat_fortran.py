import os,sys
import numpy
import math
import random
from matplotlib.pyplot import *

sys.path.append('/Users/curtisj/Desktop/smol')
import sasmol
import test

bcoh={'H':-3.739, 'D':6.671, 'C':6.646, 'N':9.29, 'O':5.803, 'S':2.847, 'P':5.13}

bcoh_res_h2o={'GLY':1.728, 'ALA':1.645, 'VAL':1.479, 'LEU':1.396, 'ILE':1.396,\
	'PHE':4.139, 'TYR':4.719, 'TRP':6.035, 'ASP':3.845, 'GLU':3.762,\
	'SER':2.225, 'THR':2.142, 'ASN':3.456, 'GLN': 3.373, 'LYS':1.586,\
	'ARG':3.466, 'HIS':4.959, 'MET':1.764, 'CYS':1.930, 'PRO':2.227, \
	'HSE':4.959}

bcoh_res_d2o={'GLY':2.769, 'ALA':2.686, 'VAL':2.520, 'LEU':2.437, 'ILE':2.437,\
	'PHE':5.180, 'TYR':6.802, 'TRP':8.118, 'ASP':4.886, 'GLU':4.803,\
	'SER':4.208, 'THR':4.224, 'ASN':6.580, 'GLN': 6.497, 'LYS':5.752,\
	'ARG':9.714, 'HIS':6.521, 'MET':2.805, 'CYS':4.013, 'PRO':2.227,\
	'HSE':4.959}

def random_vector_sphere():
	'''
	von Neumann method
	'''
	looking=True
	while(looking):
		r1=random.random()
		r2=random.random()
		r3=random.random()
		z1=1-2*r1 ; z2=1-2*r2 ; z3=1-2*r3
		sumz2=z1*z1+z2*z2+z3*z3
		if(sumz2 < 1):
			z=math.sqrt(sumz2)
			vector=[z1/z,z2/z,z3/z]
			looking=False
	return vector

def main():

	qlow=0.0 ; qhigh=0.50 ; deltaq=0.02 ; numqpoints=int(((qhigh-qlow)/deltaq)+1)
	nqvectors=2**14 #2**14 == 16384
	print 'numqpoints = ',numqpoints
	totnumq=nqvectors*numqpoints	
	print 'totnumq = ',totnumq 

	qx=numpy.zeros(totnumq,numpy.float)
	qy=numpy.zeros(totnumq,numpy.float)
	qz=numpy.zeros(totnumq,numpy.float)
	k=0
	for i in range(numqpoints):
		qmag=qlow+deltaq*i
		print 'qmag = ',qmag
		for j in range(nqvectors):
			vector=random_vector_sphere()
			qx[k]=qmag*vector[:][0]	
			qy[k]=qmag*vector[:][1]	
			qz[k]=qmag*vector[:][2]	
			k+=1
	
	#plot(qx,qy,label='qx,qy')
	#draw()
	#show()

	filename='lys1.pdb'

	mol1=sasmol.SasMol(0)
	mol1.readpdb(filename)
	name1=mol1.name()
	natoms=len(name1)
	resname=mol1.resname()

	ref_indexes=[]
	basis='CA'
	for i in xrange(natoms):
		if(basis==name1[i]): 
			ref_indexes.append(i)

	resmol=sasmol.SasMol(1)
	resmol.setCoor(mol1.coor()[ref_indexes])
	resmol.setNatoms(len(resmol.coor()))

###	OPEN	Hack with getting/setting resname

	nresname=[]
	for i in range(len(ref_indexes)):
		nresname.append(resname[ref_indexes[i]])	
	resmol.setResname(nresname)

	resnames=resmol.resname()
	nresidues=len(resnames)

	bresidue=[]
	for res in resnames:
		found=0
		for key,value in bcoh_res_d2o.iteritems():
			if key==res:
				bresidue.append(value)
				found=1
		if found==0:
			print 'no scattering length data for residue ',res
			print 'quitting program'
			sys.exit(1)

	print 'start loop over pairs'

	xcoor=resmol.coor()[:,0] ; ycoor=resmol.coor()[:,1] ; zcoor=resmol.coor()[:,2]	

	print 'xcoor = ',xcoor[-1]

	sq=numpy.zeros(numqpoints,numpy.float)
	outfile=open('sqfaad2o'+str(nqvectors)+'.dat','w')

	sq=test.pairsans(xcoor,ycoor,zcoor,qx,qy,qz,bresidue,nresidues,numqpoints,nqvectors)

	for q in range(numqpoints):
		print qlow+q*deltaq,sq[q]/sq[0]
		outfile.write('%f\t%f\n' % (qlow+q*deltaq,sq[q]/sq[0]))

	outfile.close()

if __name__=='__main__':
	
	main()

