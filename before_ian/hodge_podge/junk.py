import sys,math
import numpy,scipy
import random
sys.path.append('/Users/curtisj/Desktop/smol/')
import sasmol
import mcenergy

def overlap(allcoor,trialcoor,diameter):

	check=0
	tx=trialcoor[0] ; ty=trialcoor[1] ; tz=trialcoor[2]	
	for coor in allcoor:
		x=coor[0] ; y=coor[1] ; z=coor[2]
		dist=math.sqrt((x-tx)*(x-tx)+(y-ty)*(y-ty)+(z-tz)*(z-tz))
		if dist < diameter:
			check=1
	return check

def setup(nmolecules,eta_goal):
	
#	radius == 1 == sigma
	radius=1.0 ; diameter=2*radius

	boxvolume=((nmolecules*4.0*numpy.pi*(radius**3.0))/3.0)/eta_goal
	boxl=boxvolume**(1/3.)
	print 'boxl = ',boxl	

	allcoor=[] #numpy.zeros([nmolecules,3],numpy.float)
	for i in range(nmolecules):
		search=True		
		while(search):
			trialcoor=boxl*scipy.random.random(3)
			check=overlap(allcoor,trialcoor,diameter)
			if(check==0):
				allcoor.append(trialcoor)
				search=False

	coor=numpy.array(allcoor)-(boxl/2.0) #	-L/2 to L/2 
	
	return coor,boxl

def writexyz(coor,filename,flag):
	
	natoms=len(coor)
	if(flag=='w'):
		outfile=open(filename,'w')
	elif(flag=='a'):
		outfile=open(filename,'a')
	outfile.write("%i\n\n" % natoms)
	for i in range(natoms):
		outfile.write('%s\t%f\t%f\t%f\n' % ('C',coor[i][0],coor[i][1],coor[i][2]))
	outfile.close()

def calcenergy(coor,natoms,boxl,deltag,nbins):
	xcoor=coor[:,0] ; ycoor=coor[:,1] ; zcoor=coor[:,2]
	tgr=numpy.zeros(nbins,numpy.int)
	tgr,energy=mcenergy.ljenergy(xcoor,ycoor,zcoor,natoms,boxl,deltag,nbins)
	return tgr,energy

def mcmove(coor,ran,natoms,beta,boxl,oldenergy,gr,deltag,nbins):
	result=0

	trialcoor=numpy.array(coor[ran]+(scipy.random.random(3)-0.5)) # *delx
	pbc=-boxl*((trialcoor/boxl).round())
	trialcoor+=pbc
	temp=coor.copy()
	coor[ran]=trialcoor
	tgr,newenergy=calcenergy(coor,natoms,boxl,deltag,nbins)
	if newenergy < oldenergy:
		result=1
		oldenergy = newenergy
		gr=gr+tgr
	elif random.random() < math.exp(-beta*(newenergy-oldenergy)):
		result=1
		oldenergy = newenergy
		gr=gr+tgr
		#print newenergy
		#print 'accepted an impotance sampling move'
	else:
		coor[ran]=temp[ran]
	del temp
	return result,oldenergy,gr

def mcvol():

	return

def calcsq():
	
	return

def nptmc(nmolecules,temperature,pressure,eta_goal,nsteps,nbins,runfilename):

	coor,boxl=setup(nmolecules,eta_goal)	
	writexyz(coor,runfilename,'w')
	natoms=nmolecules 	# for now
	deltag=boxl/(2.0*nbins)
	m=0 ; v=0 ; accmove=0 ; accvol=0

	gr=numpy.zeros(nbins,numpy.float)
	beta=1.0/temperature	# i.e. kb=1.0
	oldenergy=1E16		# initial energy is artifically high
	print 'starting MC run'
	for i in xrange(nsteps):
		
		if ((i+1)%(nsteps/10.))==0: 
			print 'i = ',i+1,' out of ',nsteps,' ',100.0*(i+1)/nsteps,' percent done'
			print '\n\npercent mcmoves accepted ',accmove*100.0/m
		ran=random.randint(0,nmolecules+1)
		if ran < nmolecules:
			result,oldenergy,gr=mcmove(coor,ran,natoms,beta,boxl,oldenergy,gr,deltag,nbins)
			if result==1:
				writexyz(coor,runfilename,'a')
				accmove+=1
			else:
				pass
			m+=1
		else:
			mcvol()
			accvol+=1
			v+=1

	goutfile=open('gr.dat','w')
	for i in xrange(nbins):
		r=deltag*(i+0.5)
		vb=((i+1)**3-i**3)*deltag**3
		nid=(4./3)*numpy.pi*vb*eta_goal # eta_goal should be rho
		gr[i]=gr[i]/(accmove*natoms*nid)
		goutfile.write('%i\t%f\n' % (i,gr[i]))

	goutfile.close()

	print '\n\npercent mcmoves accepted ',accmove*100.0/m
	print '\n\npercent mcvols accepted ',accvol*100.0/v

	print '\n\npercent mcmove =  ',100.0*m/nsteps
	print 'should be nmolecules/(nmolecules+1) ',100.0*(nmolecules/(nmolecules+1.))
	print 'percent mcvol =  ',100.0*v/nsteps
	print 'should be 1/(nmolecules+1) = ',100.0*1/(nmolecules+1)
	print '\n\n'

	return

if __name__=="__main__":


	real_temperature = 298.0		# Kelvin
	real_pressure = 1.0			# Atomsphere (bar)
	real_density = 200.0			# mg / mL

	real_molecular_weight = 14313.0758	# mass of lysozyme
	real_sigma = 14.10419			# Rg of lysozyme


	pressure = real_pressure*101325.0/(41.9E6) 	# reduced pressure P* = 413.52 atm
	temperature = real_temperature/119.8		# reduced temperature T* = 119.8 K
	#density = 	


	eta_goal=0.06

###	OPEN	Eventually eta_goal won't really be used

	runfilename='blah_test_'+str(eta_goal)+'.xyz'
	nsteps=1000000
	nmolecules=108
	nbins=1000 		# number of bins for g(r) calculation

	temperature=2.48  	# T*=1.0 	<---> T = 119.8 K (2.48 = 298 K)
	density=1.0		# rho*=1.0	<---> rho = 1680 kg/m^3
	dt=0.005		# dt*=0.005	<---> dt = 1.09E-14
	pressure=1.0		# P*=1.0	<---> P = 41.9 MPa			

	nptmc(nmolecules,temperature,pressure,eta_goal,nsteps,nbins,runfilename)



	# for Argon, E/K_b = 119.8 K, sigma = 3.405E-10 m, M = 0.03994 kg/mol

	# length 	--> sigma
	# energy 	--> e
	# mass		--> m
	# time = sigma*sqrt(m/e)

	# potential energy == U* 	--> U/e
	# pressure == P* 		--> P * sigma^3 / e	
	# density == rho*		--> rho * sigma^3
	# temperature == T*		--> k_b * T / e

