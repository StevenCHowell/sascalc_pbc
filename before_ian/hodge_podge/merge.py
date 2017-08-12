import string,locale

infile1=open('ball.dat','r').readlines()
infile2=open('int.dat','r').readlines()
outfile=open('sum.dat','w')

for i in xrange(len(infile1)):
	lin1=string.split(infile1[i])
	lin2=string.split(infile2[i])
	ti1=locale.atof(lin1[1])
	ti2=locale.atof(lin2[1])	
	#outfile.write('%s\t%f\n' % (lin1[0],((ti1*108)+(ti1*ti2))))
	#outfile.write('%s\t%f\n' % (lin1[0], ti1*(108+ti2)))
	outfile.write('%s\t%f\n' % (lin1[0], ti1*(108-ti2)))

outfile.close()

