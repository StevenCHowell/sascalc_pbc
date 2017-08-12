import sys,math
import numpy,scipy
import random
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


def makefcc(max_diameter_molecule,number_of_molecules):

        pi=numpy.pi
        eta=0.74048

        matms=9
        natms=4*matms*matms*matms

          #natms = number_of_molecules

        volp=natms*4.0*pi*((max_diameter_molecule/2.0)**3.0)/3.0
        volbox=volp/eta
#        volbox=volp*natms

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

        com_coor=numpy.zeros((natms,3),numpy.float)
        com_coor[:,0]=x ; com_coor[:,1]=y ; com_coor[:,2]=z

        write_xyz('dum.xyz',com_coor,boxl)

        return com_coor,boxl




def setup(nmolecules,rho_start):
	
	radius = 0.208 #  argon tomic radius / sigma argon

	boxvolume=((nmolecules*4.0*numpy.pi*(radius**3.0))/3.0)/rho_start

	molvolume = ((nmolecules*4.0*numpy.pi*(radius**3.0))/3.0)
	boxvolume = nmolecules/rho_start

	packing_fraction = molvolume/boxvolume

	print 'packing fraction = ',packing_fraction
	print 'boxvolume = ',boxvolume
	
	boxl=boxvolume**(1/3.)
	print 'boxl = ',boxl	
	print 'radius = ',radius
	diameter=2.0*radius

	allcoor=[] #numpy.zeros([nmolecules,3],numpy.float)


	matms = 8
	natms = 4*matms*matms*matms

	if natms!=nmolecules:
		print 'error: natms = ',natms,' nmolecules = ',nmolecules
		sys.exit()

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

        allcoor=numpy.zeros((natms,3),numpy.float)
        allcoor[:,0]=x ; allcoor[:,1]=y ; allcoor[:,2]=z

#	for i in range(nmolecules):
#		search=True		
#		while(search):
#			trialcoor=boxl*scipy.random.random(3)
#			check=overlap(allcoor,trialcoor,diameter)
#			if(check==0):
#				allcoor.append(trialcoor)
#				search=False


	coor=numpy.array(allcoor)-(boxl/2.0) #	-L/2 to L/2 
	
	return coor,boxl,radius

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
	tgr,energy,virial=mcenergy.ljenergy(xcoor,ycoor,zcoor,natoms,boxl,deltag,nbins)
	return tgr,energy,virial

def mcmove(coor,ran,natoms,beta,boxl,oldenergy,gr,deltag,nbins):
	result=0

	trialcoor=numpy.array(coor[ran]+(scipy.random.random(3)-0.5)) # *delx
	pbc=-boxl*((trialcoor/boxl).round())
	trialcoor+=pbc
	temp=coor.copy()
	coor[ran]=trialcoor
	tgr,newenergy,virial=calcenergy(coor,natoms,boxl,deltag,nbins)
	volume=boxl**3.0
	rho=natoms/volume
	instpressure=0.0
	if newenergy < oldenergy:
		result=1
		oldenergy = newenergy
		gr=gr+tgr
		instpressure=(rho/beta)+virial/volume
	elif random.random() < math.exp(-beta*(newenergy-oldenergy)):
		result=1
		oldenergy = newenergy
		gr=gr+tgr
		instpressure=(rho/beta)+virial/volume
		#print newenergy
		#print 'accepted an impotance sampling move'
	else:
		coor[ran]=temp[ran]
	del temp


	return result,oldenergy,gr,instpressure

def nptmc(nmolecules,temperature,pressure,rho_start,nsteps,nbins,runfilename):

	coor,boxl,radius=setup(nmolecules,rho_start)	
	writexyz(coor,runfilename,'w')
	writexyz(coor/boxl,'test.xyz','w') 
	natoms=nmolecules 	# for now
	deltag=boxl/(2.0*nbins)
	m=0 ; v=0 ; accmove=0 ; accvol=0

	boxlmax=boxl*1.05
	vmax=boxlmax**3.0

	gr=numpy.zeros(nbins,numpy.float)
	beta=1.0/temperature	# i.e. kb=1.0
	oldenergy=1E16		# initial energy is artifically high
	print 'starting MC run'
	boxfile=open('box.dat','w')
	pressfile=open('press.dat','w')
	boxsum=0.0 ; presssum=0.0 ; numpress=0
	result=0 ; vresult=0
	for i in xrange(nsteps):
	
		if(i<100):
			print 'step = ',i,
			sys.stdout.flush()	
	
		if ((i+1)%(nsteps/10.))==0: 
			print 'i = ',i+1,' out of ',nsteps,' ',100.0*(i+1)/nsteps,' percent done'
			print '\n\npercent mcmoves accepted ',accmove*100.0/m
		ran=random.randint(0,nmolecules+1)
		if ran < nmolecules:
			result,oldenergy,gr,ip=mcmove(coor,ran,natoms,beta,boxl,oldenergy,gr,deltag,nbins)
			if result==1:
				writexyz(coor,runfilename,'a')
				accmove+=1
				if(ip<100):
					numpress+=1
					presssum+=ip
					avgp=presssum/numpress
					pressfile.write("%i\t%f\t%f\n" % (i,ip,avgp))
					pressfile.flush()
			else:
				pass
			m+=1
		else:
			pass
	boxfile.close()
	pressfile.close()
	rho=1.0/(boxl**3.0)
	goutfile=open('gr.dat','w')
	for i in xrange(nbins):
		r=deltag*(i+0.5)
		vb=((i+1)**3-i**3)*deltag**3
		nid=(4./3)*numpy.pi*vb*rho
		gr[i]=gr[i]/(accmove*natoms*nid)
		goutfile.write('%i\t%f\t%f\n' % (i,r,gr[i]/nmolecules))

	goutfile.close()

	print '\n\npercent mcmoves accepted ',accmove*100.0/m

	print '\n\npercent mcmove =  ',100.0*m/nsteps
	print 'should be nmolecules/(nmolecules+1) ',100.0*(nmolecules/(nmolecules+1.))
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

	rho_start=0.8442

	runfilename='blah_test_'+str(rho_start)+'.xyz'
	nsteps=100000
	nmolecules=2000
	nmolecules=2048


	nbins=1000 		# number of bins for g(r) calculation

	temperature=2.48  	# T*=1.0 	<---> T = 119.8 K (2.48 = 298 K)
	density=1.0		# rho*=1.0	<---> rho = 1680 kg/m^3
	dt=0.005		# dt*=0.005	<---> dt = 1.09E-14
	pressure=1.0		# P*=1.0	<---> P = 41.9 MPa			

	temperature=2.0
	temperature=1.5043

	rho_start = 0.4
	temperature = 1.51

	nptmc(nmolecules,temperature,pressure,rho_start,nsteps,nbins,runfilename)

	# for Argon, E/K_b = 119.8 K, sigma = 3.405E-10 m, M = 0.03994 kg/mol

	# length 	--> sigma
	# energy 	--> e
	# mass		--> m
	# time = sigma*sqrt(m/e)

	# potential energy == U* 	--> U/e
	# pressure == P* 		--> P * sigma^3 / e	
	# density == rho*		--> rho * sigma^3
	# temperature == T*		--> k_b * T / e

