import string,locale
from scipy.integrate import simps
from scipy.integrate import trapz
import numpy


def integrand(i,q,r,gr):

	value = (gr[i] - 1.0) * r[i] * numpy.sin(q * r[i])
	
	return value

def read_gr(filename):

	data=open(filename,'r').readlines()
	nl=len(data)
	r=numpy.zeros(nl,numpy.float)
	gr=numpy.zeros(nl,numpy.float)


	sigma = 3.405

	for i in xrange(nl):
		lin=string.split(data[i])
		tr=locale.atof(lin[1])
		tgr=locale.atof(lin[2])
		r[i]=tr*3.405 ; gr[i]=tgr

	return r,gr

if __name__=='__main__': 

	pi = numpy.pi

	filename = 'argon_gr.txt'
	filename = 'gr.dat'
	ofilename = filename+'.int'

	r,gr = read_gr(filename)

	highq = 12.0
	nsteps = len(r)
	stepq = (highq)/nsteps

	lowq = stepq
	
	outfile=open(ofilename,'w')

	eta = 0.0318217 / 2
	no = eta
	fact = 4 * pi * no
	fact = eta

	rho = 864.0/(10.256618**3.0)
	rho = 1.0/(10.256618**3.0)

	for s in xrange(nsteps+1):
		this_q = lowq + stepq*s
		value = numpy.zeros(len(r),numpy.float)
		for i in xrange(len(r)):
			value[i] = integrand(i,this_q,r,gr)

		sq = (fact*simps(value)/this_q) + 1
		sq2 = (fact*trapz(value)/this_q) + 1

		outfile.write('%f\t%f\t%f\n' % (this_q,sq,sq2))

	outfile.close()

