import os,sys
import numpy
import math
import random
from matplotlib.pyplot import *

sys.path.append('/Users/curtisj/Desktop/smol')
import sasmol
import test2

bcoh={'H':-3.739, 'D':6.671, 'C':6.646, 'N':9.29, 'O':5.803, 'S':2.847, 'P':5.13}

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

	boxlength=19.6089213591
	qlow=0.0 ; qhigh=10.0 ; deltaq=0.1 ; numqpoints=int(((qhigh-qlow)/deltaq)+1)
	nqvectors=2048
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
	filename='frame_333.pdb'

	mol1=sasmol.SasMol(0)
	mol1.readpdb(filename)
	element=mol1.element()
	natoms=len(element)

	abatom=[]
	for atom in element:
		found=0
		for key,value in bcoh.iteritems():
			if key==atom:
				abatom.append(value)
				found=1
		if found==0:
			print 'no scattering length data for atom ',atom
			print 'quitting program'
			sys.exit(1)

	batom=numpy.array(abatom)

	print 'start loop over pairs'

	xcoor=mol1.coor()[:,0] ; ycoor=mol1.coor()[:,1] ; zcoor=mol1.coor()[:,2]	

	print 'xcoor = ',xcoor[-1]

	sq=numpy.zeros(numqpoints,numpy.float)
	outfile=open('sqf2048.dat','w')

	sq=test2.pairsans(xcoor,ycoor,zcoor,qx,qy,qz,batom,natoms,numqpoints,nqvectors,boxlength)

	for q in range(numqpoints):
		print qlow+q*deltaq,sq[q]
		outfile.write('%f\t%f\n' % (qlow+q*deltaq,sq[q]))

	outfile.close()

if __name__=='__main__':
	
	main()

