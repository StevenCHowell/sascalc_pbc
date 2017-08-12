import os,sys
import numpy
import math
import random
from matplotlib.pyplot import *

sys.path.append('/Users/curtisj/Desktop/smol')
import sasmol

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
	
	qlow=0.0 ; qhigh=0.30 ; deltaq=0.02 ; numqpoints=int(((qhigh-qlow)/deltaq)+1)
	nqvectors=4
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
	element=mol1.element()
	natoms=len(element)

	batom=[]
	for atom in element:
		found=0
		for key,value in bcoh.iteritems():
			if key==atom:
				batom.append(value)
				found=1
		if found==0:
			print 'no scattering length data for atom ',atom
			print 'quitting program'
			sys.exit(1)

	print 'start loop over pairs'

	xcoor=mol1.coor()[:,0] ; ycoor=mol1.coor()[:,1] ; zcoor=mol1.coor()[:,2]	

	print 'xcoor = ',xcoor[-1]

	sq=numpy.zeros(numqpoints,numpy.float)
	outfile=open('sq.dat','w')

	for i in range(natoms):
		b1=batom[i]
		xi=xcoor[i] ; yi=ycoor[i] ; zi=zcoor[i]
		for j in range(natoms):
			b2=batom[j]
			xj=xcoor[j] ; yj=ycoor[j] ; zj=zcoor[j]
			qs=0 ; qp=0
			for q in range(numqpoints):
				sum=0.0
				for qvec in range(nqvectors):
					tqx=qx[qs] ; tqy=qy[qs] ; tqz=qz[qs]
					u=xi*tqx+yi*tqy+zi*tqz
					v=xj*tqx+yj*tqy+zj*tqz
					cosu=math.cos(u) ; sinu=math.sin(u)
					cosv=math.cos(v) ; sinv=math.sin(v)
					sum=sum+b1*b2*(cosu*cosv+sinu*sinv)
					qs+=1
				sq[qp]+=sum
				qp+=1

	for qp in range(numqpoints):
		sq[qp]=sq[qp]/(nqvectors*natoms)
		outfile.write('%f\t%f\n' % (qlow+qp*deltaq,sq[qp]))

	outfile.close()

if __name__=='__main__':
	
	main()

