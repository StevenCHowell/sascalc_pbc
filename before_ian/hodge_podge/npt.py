import sys,math,locale,string
import numpy,scipy
import random
import mcenergy
import mcvolenergy

def overlap(allcoor,trialcoor,diameter):

	check=0
	tx=trialcoor[0] ; ty=trialcoor[1] ; tz=trialcoor[2]	
	for coor in allcoor:
		x=coor[0] ; y=coor[1] ; z=coor[2]
		dist=math.sqrt((x-tx)*(x-tx)+(y-ty)*(y-ty)+(z-tz)*(z-tz))
		if dist < diameter:
			check=1
	return check

def readxyz(filename):
	
	print '\n>> reading in coordinates from: ',filename	
	radius = 0.208 #  argon atomic radius / sigma argon
	infile=open(filename,'r').readlines()
	
	allcoor=[]
	snatoms=string.split(infile[0]) 
	natoms=locale.atoi(snatoms[0])
	j=0 ; burp=0
	for i in xrange(len(infile)):
		tline=string.split(infile[i])
		if tline:
			if(tline[0] == 'C'):
				x=locale.atof(tline[1])			
				y=locale.atof(tline[2])			
				z=locale.atof(tline[3])			
				allcoor.append([x,y,z])	
				j+=1	
			if(tline[0] == 'B'):
				boxl=locale.atof(tline[1])			

	coor=numpy.array(allcoor)
	
	molvolume = ((natoms*4.0*numpy.pi*(radius**3.0))/3.0)
	boxvolume = boxl**3.0
	print '\nnatoms listed in header = ',natoms
	print 'read in = ',j,' atoms'
	print 'boxl = ',boxl,'\n'
	packing_fraction = molvolume/boxvolume
	print '\npacking fraction = ',packing_fraction
	print 'boxvolume = ',boxvolume
	print 'density = ',natoms/boxvolume

	return coor,radius,boxl

def setup(nmolecules,rho_start):
	
	radius = 0.208 #  argon atomic radius / sigma argon

	molvolume = ((nmolecules*4.0*numpy.pi*(radius**3.0))/3.0)
	boxvolume = nmolecules/rho_start
	packing_fraction = molvolume/boxvolume
	print 'packing fraction = ',packing_fraction
	print 'boxvolume = ',boxvolume
	print 'density = ',nmolecules/boxvolume
	
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

	coor=numpy.array(allcoor)-(boxl/2.0) #	-L/2 to L/2 
	
	return coor,boxl,radius

def writefxyz(coor,filename,boxl):
	
	natoms=len(coor)
	outfile=open(filename,'w')
	outfile.write("%i\n\n" % natoms)
	for i in range(natoms):
		outfile.write('%s\t%f\t%f\t%f\n' % ('C',coor[i][0],coor[i][1],coor[i][2]))
	outfile.write('%s\t%f\n' % ('B',boxl))
	outfile.close()

	return

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
	
	return

def calcenergy(coor,natoms,boxl,deltag,nbins):
	xcoor=coor[:,0] ; ycoor=coor[:,1] ; zcoor=coor[:,2]
	tgr=numpy.zeros(nbins,numpy.int)
	tgr,energy,virial=mcenergy.ljenergy(xcoor,ycoor,zcoor,natoms,boxl,deltag,nbins)
	return tgr,energy,virial

def calcenergy2(coor,natoms,boxl):
	#print 'starting a calcenergy2 query'
	xcoor=coor[:,0] ; ycoor=coor[:,1] ; zcoor=coor[:,2]
	#energy=mcvolenergy.ljvenergy(xcoor,ycoor,zcoor,natoms,boxl)
	pi=numpy.pi

	energy=0.0
	density=natoms/(boxl**3.0)
	corr=8.0*pi*density/3.0
	for i in xrange(natoms-1):
		xi=xcoor[i]
		yi=ycoor[i]
		zi=zcoor[i]
		for j in range(i+1,natoms):
			xj=xcoor[j]
			yj=ycoor[j]
			zj=zcoor[j]
                        rxij=xi-xj
                        ryij=yi-yj
                        rzij=zi-zj

			rxij=rxij-boxl*((rxij/boxl).round())
			ryij=ryij-boxl*((ryij/boxl).round())
			rzij=rzij-boxl*((rzij/boxl).round())
                        
			rij=math.sqrt(rxij*rxij+ryij*ryij+rzij*rzij)
                        rij2=rij*rij
                        if(rij < boxl/2.0):
                                irij2=1.0/rij2
                                irij6=irij2**3.0
                                irij12=irij6**2.0
                                val=4.0*(irij12-irij6) 
                        else:
                                irij=1.0/rij
                                irij3=irij**3.0
                                irij9=irij3**3.0
                                val=corr*((irij9/3.0)-irij3) 
                        
                        energy=energy+val
	#print 'done with a calcenergy2 query'
	return energy

def mcmove(coor,ran,natoms,beta,boxl,oldenergy,gr,deltag,nbins,delta_mov):
	result=0

	#trialcoor=numpy.array(coor[ran]+(scipy.random.random(3)-0.5)) # *delx
	trialcoor=numpy.array(coor[ran]+(scipy.random.random(3)*delta_mov)) # *delx
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
		instpressure=(rho/beta)+(virial/volume)
	elif random.random() < math.exp(-beta*(newenergy-oldenergy)):
		result=1
		oldenergy = newenergy
		gr=gr+tgr
		instpressure=(rho/beta)+(virial/volume)
		#print newenergy
		#print 'accepted an impotance sampling move'
	else:
		coor[ran]=temp[ran]
	del temp

	return result,oldenergy,gr,instpressure

def mcvol(coor,natoms,beta,boxl,oldenergy,vmax,pressure):

	delta_vol=1.0
	vmax=((boxl/40.0)**3.0)*delta_vol
	#vmax=((boxl/30.0)**3.0)*delta_vol

	result=0
	vo=boxl**3.0
	lnvn=math.log(vo)+(random.random()-0.5)*vmax
	vn=math.exp(lnvn)
	boxn=vn**(1.0/3.0)
	#print 'vo = ',vo,' vn = ',vn
	#print 'boxl = ',boxl,' boxn = ',boxn,' boxn/boxl = ',boxn/boxl,'\n\n'
	trialcoor=coor.copy()	
	trialcoor=trialcoor*(boxn/boxl)
	
	enn=calcenergy2(trialcoor,natoms,boxl)
	eno=oldenergy	

	#print 'eno = ',eno,' enn = ',enn
	
	#arg=-beta*((enn-eno)+pressure*(vn-vo) + -(natoms+1)*math.log(vn/vo)/beta)
	arg=((enn-eno)+pressure*(vn-vo)-(natoms+1)*math.log(vn/vo)/beta)

	if arg < 0.0:
		oldenergy=enn
		result=1
	elif abs(arg)<100:
		value=math.exp(-arg*beta)
		if(random.random()<value):
			oldenergy=enn
			result=1
		else:
			#coor=coor*(boxl/boxn)		# restore coordinates
			result=0
	else:
		#coor=coor*(boxl/boxn)		# restore coordinates
		result=0

	if result==0:
		rho=natoms/vo
		del trialcoor
		return result,oldenergy,boxl
	else:
		coor[:]=trialcoor[:]
		#coor=trialcoor.copy()
		rho=natoms/vn
		del trialcoor
		return result,enn,boxn
	
def calcsq():
	
	return

def nptmc(nmolecules,temperature,pressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename,coorfilename,restart_coor,npt):
	
	pi=numpy.pi
	if(restart_coor == 0):
		coor,boxl,radius=setup(nmolecules,rho_start)	
	else:
		coor,radius,boxl=readxyz(coorfilename)

	writexyz(coor,runfilename,'w')
	writexyz(coor/boxl,'test.xyz','w') 
	natoms=nmolecules 	# for now

	deltag=boxl/(2.0*nbins)
	m=0 ; v=0 ; accmove=0 ; accvol=0

	vmax=(boxl**3.0)*delta_vol/100.0
	gr=numpy.zeros(nbins,numpy.float)
	beta=1.0/temperature	# i.e. kb=1.0
	oldenergy=1E16		# initial energy is artifically high
	boxfile=open('box.dat','w')
	pressfile=open('press.dat','w')
	boxsum=0.0 ; presssum=0.0 ; numpress=0
	result=0 ; vresult=0  ; avgp=0.0

	if(npt==1):
		print 'starting NPT MC run'
	else:
		print 'starting NVT MC run'
		v=1
	
	print '\nSTEP\t%ACC_MOV\t%ACC_VOL\tDENSITY\t\t<PRESSURE>\tBOXL\n'
	for i in xrange(nsteps):
		
		if ((i+1)%(nsteps/10.))==0: 
			#print 'i = ',i+1,' out of ',nsteps,' ',100.0*(i+1)/nsteps,' percent done'
			#print 'i = ',i+1,'\t',accmove*100.0/m,'\t',accvol*100.0/v,'\t',natoms/boxl**3.0,'\t',avgp,'\n'
			print '%i\t%8.6f\t%8.6f\t%8.6f\t%10.6f\t%8.6f\n' % (i+1,accmove*100.0/m,accvol*100.0/v,natoms/boxl**3.0,avgp,boxl)
			#print '\n\npercent mcmoves accepted ',accmove*100.0/m
			#if(npt==1): print '\n\npercent mcvol accepted ',accvol*100.0/v
			#if(npt==1): print '\n\ndensity = ',natoms/(boxl**3.0)
		ran=random.randint(0,nmolecules+1)
		if ran < nmolecules:
			result,oldenergy,gr,ip=mcmove(coor,ran,natoms,beta,boxl,oldenergy,gr,deltag,nbins,delta_move)
			if result==1:
				writexyz(coor,runfilename,'a')
				accmove+=1
				irho=natoms/(boxl**3.0)
				rc=boxl/2.0
				irc3=1.0/(rc**3.0)
				irc9=irc3*irc3
				ptail=((16.0*pi*irho*irho)/3.0)*((2.0*irc9/3.0)-irc3 )			
				if(ip<1000000):
					numpress+=1
					presssum+=ip+ptail
					avgp=presssum/numpress
					pressfile.write("%i\t%f\t%f\t%f\n" % (i,ip+ptail,avgp,ptail))
					pressfile.flush()
			else:
				pass
			m+=1
		elif(npt==1):
			#blah1=coor[0,0]
			vresult,oldenergy,boxl=mcvol(coor,natoms,beta,boxl,oldenergy,vmax,pressure)
				
			if vresult==1:
			#	print 'blah1[0,0] = ',blah1
			#	print 'coor[0,0] = ',coor[0,0],'\n'
				writexyz(coor,runfilename,'a')
				accvol+=1
				#print 'new boxl = ',boxl
				boxsum=boxsum+boxl
				boxfile.write("%i\t%f\t%f\n" % (i,boxl,boxsum/accvol))
				boxfile.flush()
				#print 'accepted a mcvol move',i
			else:
				#print 'failed a mcvol move',i
				pass
			v+=1

	boxfile.close()
	pressfile.close()
	if (accvol>0 and boxsum>0):
		rho=1.0/((boxsum/accvol)**3.0)
		print 'accvol = ',accvol
		print 'boxsum = ',boxsum
		print 'rho = ',rho
	else:
		rho=1.0/(boxl**3.0)
	goutfile=open('gr.dat','w')
	for i in xrange(nbins):
		r=deltag*(i+0.5)
		vb=((i+1)**3-i**3)*deltag**3
		nid=(4./3)*numpy.pi*vb*rho
		gr[i]=gr[i]/(accmove*natoms*nid)
		goutfile.write('%i\t%f\t%f\n' % (i,r,gr[i]/nmolecules))

	goutfile.close()

#	print '\n\npercent mcmoves accepted ',accmove*100.0/m
#	if(npt==1):	
#		print '\npercent mcvols accepted ',accvol*100.0/v
#		print 'percent mcvol =  ',100.0*v/nsteps,'\n'
#		print 'should be 1/(nmolecules+1) = ',100.0*1/(nmolecules+1)
#
#	print '\n\npercent mcmove =  ',100.0*m/nsteps
#	print 'should be nmolecules/(nmolecules+1) ',100.0*(nmolecules/(nmolecules+1.))
#	print '\n\n'

	writefxyz(coor,'final_'+runfilename,boxl)

	return

if __name__=="__main__":

	real_temperature = 298.0		# Kelvin
	real_pressure = 1.0			# Atomsphere (bar)
	real_density = 200.0			# mg / mL

	real_molecular_weight = 14313.0758	# mass of lysozyme
	real_sigma = 14.10419			# Rg of lysozyme

	pressure = real_pressure*101325.0/(41.9E6) 	# reduced pressure P* = 413.52 atm
	temperature = real_temperature/119.8		# reduced temperature T* = 119.8 K
	dt=0.005				# dt*=0.005	<---> dt = 1.09E-14
	temperature=2.48  	# T*=1.0 	<---> T = 119.8 K (2.48 = 298 K)
	density=1.0		# rho*=1.0	<---> rho = 1680 kg/m^3

	rho_start=0.8442
	#rho_start=0.8

	delta_move = 0.15	# translational move step size (percent)
	delta_vol = 0.001	# volume move step size (percent)

	runfilename='blah_test_'+str(rho_start)+'.xyz'
	coorfilename='final_'+runfilename

	nsteps=1000000
	nmolecules=108
	nbins=1000 		# number of bins for g(r) calculation

	pressure=1.0		# P*=1.0	<---> P = 41.9 MPa			
	#pressure=6.0		# P*=1.0	<---> P = 41.9 MPa			

	temperature=1.5043
	#temperature=2.0
	
	restart_coor=0		# 0==no, 1==yes
	npt=0		# npt=1 or npt=0==nvt

	nptmc(nmolecules,temperature,pressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename,coorfilename,restart_coor,npt)

	# for Argon, E/K_b = 119.8 K, sigma = 3.405E-10 m, M = 0.03994 kg/mol

	# length 	--> sigma
	# energy 	--> e
	# mass		--> m
	# time = sigma*sqrt(m/e)

	# potential energy == U* 	--> U/e
	# pressure == P* 		--> P * sigma^3 / e	
	# density == rho*		--> rho * sigma^3
	# temperature == T*		--> k_b * T / e

