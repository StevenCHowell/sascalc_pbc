import string,locale
from scipy import *
import numpy

data=open('gr.dat','r').readlines()
nl=len(data)

x=numpy.zeros(nl*2,numpy.float)
y=numpy.zeros(nl*2,numpy.float)

dr = 2.515615 - 2.513096

shift = 0.929761

outfile2=open('new_gr.dat','w')

for i in xrange(nl):
	lin=string.split(data[i])
	tx=locale.atof(lin[1])
	ty=locale.atof(lin[2])
	x[i]=tx ; y[i]=ty-shift
	lastx = tx ; lasty = ty-shift
	lasti = i
	outfile2.write('%f\t%f\n' % (tx,ty))

for i in xrange(nl):
	tx = lastx + dr*(i+1)
	ty = lasty
	x[lasti+i+1] = tx
	y[lasti+i+1] = ty
	outfile2.write('%f\t%f\n' % (tx,ty))

apod=[exp(-x/20.) for x in range(nl*2)]
napod=numpy.array(apod)
ay=napod*y



outfile2.close()

f=fft(y)
af=fft(ay)

outfile=open('gr.fft','w')
for i in range(len(f)):
	#outfile.write('%f\t%f\n' % (f[i].real(),f[i].imag()))
	outfile.write('%f\t%f\n' % (numpy.real(f[i]),numpy.real(af[i])))

outfile.close()

