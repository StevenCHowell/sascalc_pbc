import sys,math,locale,string
import numpy,scipy
import random
import newmcenergy

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

	return coor,radius,boxl

def setup(nmolecules,rho_start):
	
	radius = 0.208 #  argon atomic radius / sigma argon

	molvolume = ((nmolecules*4.0*numpy.pi*(radius**3.0))/3.0)
	boxvolume = nmolecules/rho_start
	packing_fraction = molvolume/boxvolume
	print 'packing fraction = ',packing_fraction
	print 'boxvolume = ',boxvolume
	
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

	#coor=numpy.array(allcoor)-(boxl/2.0) #	-L/2 to L/2 
	coor=numpy.array(allcoor)
	
	return coor,boxl,radius

def writefxyz(coor,filename,boxl):
	
	natoms=len(coor)
	outfile=open(filename,'w')
	outfile.write("%i\n" % natoms)
	outfile.write('%s\t%f\n' % ('B',boxl))
	for i in range(natoms):
		outfile.write('%s\t%f\t%f\t%f\n' % ('C',coor[i][0],coor[i][1],coor[i][2]))
	outfile.close()

	return

def writexyz(coor,boxl,filename,flag):
	
	natoms=len(coor)
	if(flag=='w'):
		outfile=open(filename,'w')
	elif(flag=='a'):
		outfile=open(filename,'a')
	outfile.write("%i\n" % natoms)
	outfile.write('%s\t%f\n' % ('B',boxl))
	for i in range(natoms):
		outfile.write('%s\t%f\t%f\t%f\n' % ('C',coor[i][0],coor[i][1],coor[i][2]))
	outfile.close()
	
	return

def calcsq():
	
	return


def sumup(coor,natoms,rcut,rmin,boxl): #,v12,v16,w12,w6):

	xcoor=coor[:,0] ; ycoor=coor[:,1] ; zcoor=coor[:,2]
	overlap = 0.0
	rcutsq = rcut*rcut
	rminsq = rmin*rmin
	boxlinv = 1.0/boxl
	
	v12 = 0.0 ; v6 = 0.0 ; w12 = 0.0 ; w6 = 0.0

	for i in xrange(natoms-1):
		rxi=xcoor[i] ; ryi=ycoor[i] ; rzi=zcoor[i]
		for j in xrange(i+1,natoms):
			rxj=xcoor[j] ; ryj=ycoor[j] ; rzj=zcoor[j]
			rxij = rxi - rxj	
			ryij = ryi - ryj	
			rzij = rzi - rzj
	
			rxij = rxij-boxl*((rxij*boxlinv).round())
			ryij = ryij-boxl*((ryij*boxlinv).round())
			rzij = rzij-boxl*((rzij*boxlinv).round())

			rijsq = rxij*rxij + ryij*ryij + rzij*rzij
					
			#if(rijsq < rminsq):
			#	overlap=1.0
			if(rijsq < rcutsq):
			#elif(rijsq < rcutsq):
				sr2 = 1.0/rijsq
				sr6 = sr2*sr2*sr2
				vij12 = sr6*sr6
				vij6 = -sr6
				v12 = v12 + vij12
				v6 = v6 + vij6
				w12 = w12 + vij12
				w6 = w6 + vij6*0.5	

		#print 'rxi = ',rxi

	v12 = 4.0*v12
	v6 = 4.0*v6
	w12 = 48.0*w12/3.0
	w6 = 48.0*w6/3.0

	return overlap,v12,v6,w12,w6

def energy(coor,natoms,rxi,ryi,rzi,i,rcut,boxl):

	xcoor=coor[:,0] ; ycoor=coor[:,1] ; zcoor=coor[:,2]
	rcutsq = rcut*rcut
	boxlinv = 1.0/boxl
	
	v12 = 0.0 ; v6 = 0.0 ; w12 = 0.0 ; w6 = 0.0

	for j in xrange(natoms):
		if(i != j):
			
			rxij = rxi-xcoor[j]		
			ryij = ryi-ycoor[j]		
			rzij = rzi-zcoor[j]		

			rxij = rxij-boxl*((rxij*boxlinv).round())
			ryij = ryij-boxl*((ryij*boxlinv).round())
			rzij = rzij-boxl*((rzij*boxlinv).round())

			rijsq = rxij*rxij + ryij*ryij + rzij*rzij

			if(rijsq < rcutsq):
				sr2 = 1.0/rijsq
				sr6 = sr2*sr2*sr2
				vij12 = sr6*sr6
				vij6 = -sr6
				v12 = v12 + vij12
				v6 = v6 + vij6
				w12 = w12 + vij12
				w6 = w6 + vij6*0.5	
				
	v12 = 4.0*v12
	v6 = 4.0*v6
	w12 = 48.0*w12/3.0
	w6 = 48.0*w6/3.0

	return v12,v6,w12,w6

def nptmc(nmolecules,temperature,goalpressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename,coorfilename,restart_coor,iratio,vratio):


	boxfile=open('box.dat','w')
	pressfile=open('press.dat','w')
	densityfile=open('density.dat','w')
	
	pi=numpy.pi
	if(restart_coor == 0):
		coor,boxl,radius=setup(nmolecules,rho_start)	
	else:
		coor,radius,boxl=readxyz(coorfilename)

	writexyz(coor,boxl,runfilename,'w')
	writexyz(coor/boxl,boxl,'test.xyz','w') 
	natoms=nmolecules 	# for now

	volume = boxl**3.0
	boxlinv = 1.0/boxl
	density = natoms/volume

	rcut = 0.5*boxl

	dboxmx = boxl/40.0
	drmax = delta_move
	rmin = 0.70
	beta = 1.0/temperature	# i.e. kb=1.0
	
	sr3 = (1.0/rcut)**3.0
	sr9 = sr3**3.0
	vlrc12 = 8.0*pi*density*natoms*sr9/9.0
	vlrc6 = -8.0*pi*density*natoms*sr3/3.0
	wlrc12 = 4.0*vlrc12
	wlrc6 =	2.0*vlrc6
	vlrc = vlrc12 + vlrc6
	wlrc = wlrc12 + wlrc6

	acm = 0.0 ; acatma = 0.0 ; acboxa = 0.0 ; acv = 0.0 ; acp = 0.0	; acd = 0.0
	acvsq = 0.0 ; acpsq = 0.0 ; acdsq = 0.0 ; flv = 0.0 ; flp = 0.0	; fld = 0.0

	overlap,v12,v6,w12,w6 = sumup(coor,natoms,rcut,rmin,boxl) #,v12,v6,w12,w6)
 
	vs = ( v12 + v6 + vlrc ) / natoms
	ws = ( w12 + w6 + wlrc ) / natoms
	ps = density * temperature + ( w12 + w6 + wlrc ) / volume

	v12 = v12 + vlrc12
	v6 = v6 + vlrc6
	w12 = w12 + wlrc12
	w6 = w6 + wlrc6

	print 'sr3 = ',sr3
	print 'sr9 = ',sr9
	print 'vlrc12 = ',vlrc12
	print 'vlrc6 = ',vlrc6
	print 'wlrc12 = ',wlrc12
	print 'wlrc6 = ',wlrc6
	print 'vlrc = ',vlrc
	print 'wlrc = ',wlrc
	print ''

	print 'v12 = ',v12
	print 'v6 = ',v6
	print 'w12 = ',w12
	print 'w6 = ',w6
	
	print 'density = ',density
	print 'initial energy/atom = ',vs
	print 'initial w/atom = ',ws
	print 'initial pressure = ',ps

	molvolume = ((natoms*4.0*numpy.pi*(radius**3.0))/3.0)

	m=0 ; v=0 ; tm=0 ; tv=0 
	sumeta=0.0

	print '\nSTEP\t\tVN\t%ACC_MOV\t%ACC_VOL\t ETA\t\t<ETA>\t\tDENSITY\t\t<DENSITY>\tPRESSURE\t<PRESSURE>\tBOXL\n'
		#print '%i\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%10.6f\t%8.6f\n' % (i+1,acatma*100.0/m,acboxa*100.0/v,density,acd/(m+v),pressure,acp/(m+v),boxl)

	xcoor=coor[:,0] ; ycoor=coor[:,1] ; zcoor=coor[:,2]
	for step in xrange(nsteps):
	
		for i in xrange(natoms):
			m=m+1 ; tm=tm+1		
			rxiold = coor[i,0]
			ryiold = coor[i,1]
			rziold = coor[i,2]
			v12old=0.0;v6old=0.0;w12old=0.0;w6old=0.0	
			#print 'rxiold = ',rxiold
			#print 'ryiold = ',ryiold
			#print 'rziold = ',rziold
			thisi=i+1	
			#v12old,v6old,w12old,w6old=energy(coor,natoms,rxiold,ryiold,rziold,i,rcut,boxl)
			#print 'v12old = ',v12old
			#print 'v6old = ',v6old
			#print 'w12old = ',w12old
			#print 'w6old = ',w6old
			
			v12old,v6old,w12old,w6old=energy(coor,natoms,rxiold,ryiold,rziold,i,rcut,boxl)
			#v12,v6,w12,w6=newmcenergy.ljenergy(xcoor,ycoor,zcoor,rxiold,ryiold,rziold,natoms,thisi,boxl)
			#v12old=v12 ; v6old=v6 ; w12old=w12 ; w6old=w6	
			#print 'v12old = ',v12old
			#print 'v6old = ',v6old
			#print 'w12old = ',w12old
			#print 'w6old = ',w6old

			rxinew = rxiold + ( 2.0 * random.random() - 1.0 ) * drmax
			ryinew = ryiold + ( 2.0 * random.random() - 1.0 ) * drmax
			rzinew = rziold + ( 2.0 * random.random() - 1.0 ) * drmax

			#rxinew = rxiold - 23.1
			#ryinew = ryiold - 23.1
			#rzinew = rziold - 23.1

			rxinew = rxinew-boxl*((rxinew*boxlinv).round())
			ryinew = ryinew-boxl*((ryinew*boxlinv).round())
			rzinew = rzinew-boxl*((rzinew*boxlinv).round())
		
			#v12,v6,w12,w6=newmcenergy.ljenergy(xcoor,ycoor,zcoor,rxinew,ryinew,rziold,natoms,ti,boxl)
			v12new,v6new,w12new,w6new=energy(coor,natoms,rxinew,ryinew,rzinew,i,rcut,boxl)
			#v12new=v12 ; v6new=v6 ; w12new=w12 ; w6new=w6	
			#print 'v12new = ',v12new
			#print 'v6new = ',v6new
			#print 'w12new = ',w12new
			#print 'w6new = ',w6new
	

			delv12 = v12new - v12old
			if(delv12 == 0.0):
				print 'v12new = ',v12new
				print 'v12old = ',v12old
				print 3/0

			delv6 = v6new - v6old
			delw12 = w12new - w12old
			delw6 = w6new - w6old
			deltv = delv12 + delv6
			deltvb = beta * deltv

			#print 'deltv = ',deltv
			#print 'deltvb = ',deltvb

			if(deltvb < 75.0):
				if(deltv < 0.0):
					v12=v12+delv12
					v6 =v6+delv6
					w12=w12+delw12
					w6 =w6+delw6
					coor[i,0]=rxinew
					coor[i,1]=ryinew
					coor[i,2]=rzinew
					acatma=acatma+1.0
				elif(math.exp(-deltvb) > random.random()):
				#elif(math.exp(-deltvb) > 0.5):
					v12=v12+delv12
					v6 =v6+delv6
					w12=w12+delw12
					w6 =w6+delw6
					coor[i,0]=rxinew
					coor[i,1]=ryinew
					coor[i,2]=rzinew
					acatma=acatma+1.0

			vn = ( v12 + v6 ) / natoms
			pressure = density * temperature + ( w12 + w6 ) / volume

			#print 'vn = ',vn
			#print 'w12 = ',w12
			#print 'w6 = ',w6
			#print 'volume = ',volume
#
#			print 'pressure = ',pressure
#			print 3/0

			acm = acm + 1.0
			acv = acv + vn
			acp = acp + pressure
			acd = acd + density
	
			acvsq = acvsq + vn**2.0
			acpsq = acpsq + pressure**2.0
			acdsq = acdsq + density**2.0

		#for i in xrange(natoms):
		v=v+1 ; tv=tv+1
		boxlnew = boxl + ( 2.0 * random.random() - 1.0 ) * dboxmx
		ratbox = boxl / boxlnew
		rrbox = 1.0 / ratbox
		rcutn = rcut * rrbox

		rat6 = ratbox ** 6.0
		rat12 = rat6 * rat6

		v12new = v12 * rat12
		v6new = v6 * rat6
		w12new = w12 * rat12
		w6new = w6 * rat6

		deltv = v12new + v6new - v12 - v6
		dpv = goalpressure * (boxlnew**3.0 - volume)
		dvol= 3.0 * temperature * natoms * math.log(ratbox)
		delthb = beta * ( deltv + dpv + dvol )


		if(delthb < 75.0):
			if(delthb < 0.0):
				v12 = v12new
				v6  = v6new
				w12 = w12new
				w6  = w6new

				coor=coor*rrbox

				boxl=boxlnew
				rcut=rcutn
				acboxa = acboxa + 1.0
			elif(math.exp(-delthb) > random.random()):
				v12 = v12new
				v6  = v6new
				w12 = w12new
				w6  = w6new

				coor=coor*rrbox

				boxl=boxlnew
				rcut=rcutn
				acboxa = acboxa + 1.0
			

		boxlinv = 1.0/boxl
		volume = boxl**3.0
		density=natoms/volume

		vn = ( v12 + v6 ) /natoms
		pressure = density * temperature + ( w12 + w6 ) / volume

		eta=molvolume/boxl**3.0

		acm = acm + 1.0
		acv = acv + vn
		acp = acp + pressure
		acd = acd + density

		sumeta = sumeta + eta
	
		acvsq = acvsq + vn**2.0
		acpsq = acpsq + pressure**2.0
		acdsq = acdsq + density**2.0
		
		writexyz(coor,boxl,runfilename,'a')

		boxfile.write("%i\t%f\n" % (step,boxl))	
		pressfile.write("%i\t%f\t%f\n" % (step,pressure,acp/(tm+tv)))	
		densityfile.write("%i\t%f\t%f\n" % (step,density,acd/(tm+tv)))	
		boxfile.flush() ; pressfile.flush() ; densityfile.flush()
		if (step%(nsteps/100.))==0: 
			print '%i\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%10.6f\t%8.6f\n' % (step,vn,acatma*100.0/m,acboxa*100.0/v,eta,sumeta/tv,density,acd/(tm+tv),pressure,acp/(tm+tv),boxl)

		if (step%iratio==0):
			ratio = acatma/(natoms*iratio)
			if (ratio > 0.5):
				drmax=drmax*1.05
			else:
				drmax=drmax*0.95
			#print 'drmax = ',drmax
			acatma = 0.0
			m=0
		if (step%vratio==0):
			bratio = acboxa/vratio
			if (bratio > 0.5):
				dboxmx=dboxmx*1.05
			else:
				dboxmx=dboxmx*0.95
			#print 'dboxmx = ',dboxmx
			acboxa = 0.0
			v=0	

	print '\n\npercent mcmoves accepted ',acatma*100.0/tm
	print '\n\npercent mcvols accepted ',acboxa*100.0/tv

	writefxyz(coor,'final_'+runfilename,boxl)

	boxfile.close()
	pressfile.close()

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
	rho_start=1.5917

	delta_move = 0.02	# translational move step size (angstrom)
	delta_vol = 0.001	# volume move step size (percent)

	runfilename='blah_test_'+str(rho_start)+'.xyz'
	coorfilename='final_'+runfilename
	restart_coor=1		# 0==no, 1==yes

	nsteps=500
	nmolecules=108
	nbins=1000 		# number of bins for g(r) calculation

	temperature=2.48  	# T*=1.0 	<---> T = 119.8 K (2.49 = 298 K)
	density=1.0		# rho*=1.0	<---> rho = 1680 kg/m^3
	dt=0.005		# dt*=0.005	<---> dt = 1.09E-14
	pressure=50.0		# P*=1.0	<---> P = 41.9 MPa			
	#pressure=5.2		# P*=1.0	<---> P = 41.9 MPa			

	temperature=1.5043
	temperature=2.0
	temperature=2.49	# T=298K
	temperature=1.62	# T=193K (-80 C)

	Tcentigrade = 25.0
	Tkelvin = Tcentigrade + 273.15 
	Treduced = Tkelvin/119.8
	
	temperature = Treduced
	
	print 'T(C) = ',Tcentigrade
	print 'T(K) = ',Tkelvin
	print 'T*   = ',temperature

	#Pa_to_atm = 1.0/101325.0  # 1 atm = 101325 Pa
	
	print 'P(MPa) = ',pressure*41.9
	print 'P(atm) = ',pressure*41.9*9.869233

	iratio=nsteps/100
	vratio=nsteps/100

	nptmc(nmolecules,temperature,pressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename,coorfilename,restart_coor,iratio,vratio)

	# for Argon, E/K_b = 119.8 K, sigma = 3.405E-10 m, M = 0.03994 kg/mol

	# length 	--> sigma
	# energy 	--> e
	# mass		--> m
	# time = sigma*sqrt(m/e)

	# potential energy == U* 	--> U/e
	# pressure == P* 		--> P * sigma^3 / e	
	# density == rho*		--> rho * sigma^3
	# temperature == T*		--> k_b * T / e

