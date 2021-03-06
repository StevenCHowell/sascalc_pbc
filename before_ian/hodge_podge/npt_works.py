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

def calcenergy2(coor,natoms,boxl):
	xcoor=coor[:,0] ; ycoor=coor[:,1] ; zcoor=coor[:,2]
	nnbins=1000 
	tgr=numpy.zeros(nnbins,numpy.int)
	ddeltag=boxl/(2.0*nnbins)
	tgr,energy,virial=mcenergy.ljenergy(xcoor,ycoor,zcoor,natoms,boxl,ddeltag,nnbins)
	del tgr
	
	return energy

def mcmove(coor,ran,natoms,beta,boxl,oldenergy,gr,deltag,nbins,delta_mov):
	result=0

	#trialcoor=numpy.array(coor[ran]+(scipy.random.random(3)-0.5)) # *delx
	trialcoor=numpy.array(coor[ran]+(scipy.random.random(3)-delta_mov)) # *delx
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

def mcvol2(coor,natoms,beta,boxl,oldenergy,gr,deltag,nbins,vmax,pressure):

	result=0
	vo=boxl**3.0
	logvmax=math.log(vmax)
	lnvn=math.log(vo)+(random.random()-0.5)*logvmax
	vn=math.exp(lnvn)

	boxn=vn**(1.0/3.0)
	trialcoor=coor.copy()	
	trialcoor=trialcoor*(boxn/boxl) # *delx
	
	tgr,enn,virial=calcenergy(trialcoor,natoms,boxl,deltag,nbins)
	eno=oldenergy	
	
	arg=-beta*((enn-eno)+pressure*(vn-vo) + -(natoms+1)*math.log(vn/vo)/beta)

	if abs(arg)<10:
		value=math.exp(arg)
		#print 'value = ',value
		#print 'eno = ',eno
		#print 'enn = ',enn
		#print 'vn = ',vn
		#print 'vo = ',vo
		#print 'arg = ',arg
		#print 'beta = ',beta
		if(random.random()>value):
			#coor=coor*(boxl/boxn)		# restore coordinates
			result=0
		else:
			#print 'accepted a mcvol move in mcvol'
			#boxl=boxn
			oldenergy=enn
			#print 'boxl = ',boxl
			result=1
	else:
		#coor=coor*(boxl/boxn)		# restore coordinates
		result=0

	if result==0:
		rho=natoms/vo
		instpressure=(rho/beta)+virial/vo
		return result,oldenergy,gr,instpressure,boxl
	else:
		coor[:]=trialcoor[:]
		rho=natoms/vn
		instpressure=(rho/beta)+virial/vn
		return result,enn,gr,instpressure,boxn
	

def mcvol(coor,natoms,beta,boxl,oldenergy,vmax,pressure):

	vo=boxl**3.0
	if(vmax>0.0):
		logvmax=math.log(vmax)
	else:
		print 'vmax = ',vmax
		print 'STOPPING HERE', 3/0

	result=0
	if(vo>0.0):
		lnvn=math.log(vo)+(random.random()-0.5)*logvmax
	else:
		print 'vo = ',vo
		print 'STOPPING HERE', 3/0
	if(abs(lnvn)<100.0):
		vn=math.exp(lnvn)
	else:
		print 'lnvn = ',lnvn
		print 'STOPPING HERE', 3/0

	boxn=vn**(1.0/3.0)

	trialcoor=coor.copy()	
	trialcoor=trialcoor*(boxn/boxl) # *delx

#	pbc=-boxn*((tempcoor/boxn).round())
#	tempcoor+=pbc

	enn=calcenergy2(trialcoor,natoms,boxn)
	#enn=oldenergy	
	eno=oldenergy	

	#print 'eno = ',eno
	#print 'enn = ',enn
	#print 'vn = ',vn
	#print 'vo = ',vo
	#print 'vn/vo = ',vn/vo

	if((vn/vo)>0.0):
		arg=-beta*((enn-eno)+pressure*(vn-vo) + -(natoms+1)*math.log(vn/vo)/beta)
	else:
		arg=10000;

	if abs(arg)<10:
		value=math.exp(arg)
		#value=0.8
		#print 'value = ',value
		#print 'eno = ',eno
		#print 'enn = ',enn
		#print 'vn = ',vn
		#print 'vo = ',vo
		#print 'arg = ',arg
		#print 'beta = ',beta
		if(random.random()>value):
			#coor=coor*(boxl/boxn)		# restore coordinates
			result=0
		else:
			#print 'accepted a mcvol move in mcvol'
			boxl=boxn
			oldenergy=enn
			#print 'boxl = ',boxl
			result=1

	else:
		#coor=coor*(boxl/boxn)		# restore coordinates
		result=0
	if result==0:
		return result,oldenergy,boxl
	else:
		coor[:]=trialcoor[:]
		return result,oldenergy,boxl
		
def calcsq():
	
	return

def nptmc(nmolecules,temperature,pressure,rho_start,nsteps,nbins,delta_mov,delta_vol,runfilename):

	coor,boxl,radius=setup(nmolecules,rho_start)	
	writexyz(coor,runfilename,'w')
	writexyz(coor/boxl,'test.xyz','w') 
	natoms=nmolecules 	# for now
	deltag=boxl/(2.0*nbins)
	m=0 ; v=0 ; accmove=0 ; accvol=0

	vmax=(boxl**3.0)*delta_vol/100.0

	gr=numpy.zeros(nbins,numpy.float)
	beta=1.0/temperature	# i.e. kb=1.0
	oldenergy=1E16		# initial energy is artifically high
	print 'starting MC run'
	boxfile=open('box.dat','w')
	pressfile=open('press.dat','w')
	boxsum=0.0 ; presssum=0.0 ; numpress=0
	result=0 ; vresult=0 
	npt=1
	for i in xrange(nsteps):
		
		if ((i+1)%(nsteps/10.))==0: 
			print 'i = ',i+1,' out of ',nsteps,' ',100.0*(i+1)/nsteps,' percent done'
			print '\n\npercent mcmoves accepted ',accmove*100.0/m
			if(npt==1): print '\n\npercent mcvol accepted ',accvol*100.0/v
		ran=random.randint(0,nmolecules+1)
		if ran < nmolecules:
			result,oldenergy,gr,ip=mcmove(coor,ran,natoms,beta,boxl,oldenergy,gr,deltag,nbins,delta_mov)
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
		elif(npt==1):
			#blah1=coor[0,0]
			vresult,oldenergy,gr,ip,boxl=mcvol2(coor,natoms,beta,boxl,oldenergy,gr,deltag,nbins,vmax,pressure)
			#vresult,oldenergy,boxl=mcvol(coor,natoms,beta,boxl,oldenergy,vmax,pressure)
				
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

	print '\n\npercent mcmoves accepted ',accmove*100.0/m
	if(npt==1):	
		print '\n\npercent mcvols accepted ',accvol*100.0/v
		print 'percent mcvol =  ',100.0*v/nsteps
		print 'should be 1/(nmolecules+1) = ',100.0*1/(nmolecules+1)

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

	delta_move = 0.5	# translational move step size (angstrom)
	delta_vol = 0.0001	# volume move step size (percent)

	runfilename='blah_test_'+str(rho_start)+'.xyz'
	nsteps=100000
	nmolecules=108
	nmolecules=864
	nbins=1000 		# number of bins for g(r) calculation

	temperature=2.48  	# T*=1.0 	<---> T = 119.8 K (2.48 = 298 K)
	density=1.0		# rho*=1.0	<---> rho = 1680 kg/m^3
	dt=0.005		# dt*=0.005	<---> dt = 1.09E-14
	pressure=1.0		# P*=1.0	<---> P = 41.9 MPa			
	pressure=2.0		# P*=1.0	<---> P = 41.9 MPa			

	temperature=2.0
	temperature=1.5043

	temperature = 85.0/119.8
	
	density = 0.797381  # (1339.6 kg/m^3)
	rho_start = density


	nptmc(nmolecules,temperature,pressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename)



	# for Argon, E/K_b = 119.8 K, sigma = 3.405E-10 m, M = 0.03994 kg/mol

	# length 	--> sigma
	# energy 	--> e
	# mass		--> m
	# time = sigma*sqrt(m/e)

	# potential energy == U* 	--> U/e
	# pressure == P* 		--> P * sigma^3 / e	
	# density == rho*		--> rho * sigma^3
	# temperature == T*		--> k_b * T / e

