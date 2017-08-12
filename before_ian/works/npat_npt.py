import sys,os,math,locale,string
import numpy,scipy
import random
import nenergy

from sassie.sasmol import sasmol as sasmol

def parameters():

	nsteps=1000
	nmolecules=108
	nmolecules=864
	
	start_type=3		# 1--> make fcc coords, 2--> restart from molecular coords, 3--> restart from LJ coords
				# NOTE: start_type 2 is not implemented!

	rho_start=1.5917
	rho_start=0.797381
	rho_start=0.8442
	
	runfilename='blah_test_'+str(rho_start)+'.xyz'
	coorfilename='final_'+runfilename
	coorfilename="final_blah_test_0.797381.xyz"
	coorfilename="final_0.8442.xyz"
	coorfilename="final_blah_test_0.8442.xyz"
	
	real_temperature = 298.0		# Kelvin
	real_pressure = 1.0			# Atomsphere (bar)
	real_density = 200.0			# mg / mL

	real_molecular_weight = 14313.0758	# mass of lysozyme
	real_sigma = 14.10419			# Rg of lysozyme
	
	masterpdb='ca_fcc_lys1.pdb' 		# for start_type 2 (fcc pdb)
	moleculepdb='ca_lys1.pdb' 		# single protein pdb
	specificmolvol = 9769.0  		# cm^3/mol for lysozyme

	pressure = real_pressure*101325.0/(41.9E6) 	# reduced pressure P* = 413.52 atm
	temperature = real_temperature/119.8		# reduced temperature T* = 119.8 K

	nbins=1000 		# number of bins for g(r) calculation

	density=1.0		# rho*=1.0	<---> rho = 1680 kg/m^3  !!! not used ... here for reference
	goalpressure=1.0	# P*=1.0	<---> P = 41.9 MPa			

	Tcentigrade = 25.0
	Tkelvin = Tcentigrade + 273.15 
	Tkelvin = 85.0
	Treduced = Tkelvin/119.8
	
	temperature = Treduced
	
	print 'T(C) = ',Tcentigrade
	print 'T(K) = ',Tkelvin
	print 'T*   = ',temperature

	print 'P(MPa) = ',goalpressure*41.9
	print 'P(atm) = ',goalpressure*41.9*9.869233

	delta_move = 0.02	# translational move step size (angstrom)
	delta_vol = 0.001	# volume move step size (percent)
	
	iratio=nsteps/100
	vratio=nsteps/100

	return nmolecules,temperature,goalpressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename,coorfilename,iratio,vratio,masterpdb,moleculepdb,specificmolvol,start_type

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
	radius = 0.5521 #  argon VDW radius / sigma argon
	
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

	return coor,boxl

def makefcc(diameter):

        pi=numpy.pi
        eta=0.74048

        matms=3
        natms=4*matms*matms*matms

        volp=natms*4.0*pi*((diameter/2.0)**3.0)/3.0
        volbox=volp/eta
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

	coor=numpy.zeros((natms,3),numpy.float)
        coor[:,0]=x ; coor[:,1]=y ; coor[:,2]=z

        return coor,boxl

def setup(nmolecules,rho_start):
	
	radius = 0.5521 #  argon VDW radius / sigma argon

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
	sigmaAr=3.405	
	natoms=len(coor)
	if(flag=='w'):
		outfile=open(filename,'w')
	elif(flag=='a'):
		outfile=open(filename,'a')
	outfile.write("%i\n" % natoms)
	outfile.write('%s\t%f\n' % ('B',boxl))
	for i in range(natoms):
		#outfile.write('%s\t%f\t%f\t%f\n' % ('C',coor[i][0],coor[i][1],coor[i][2]))
		outfile.write('%s\t%f\t%f\t%f\n' % ('Ar',coor[i][0]*sigmaAr,coor[i][1]*sigmaAr,coor[i][2]*sigmaAr))
	outfile.close()
	
	return

def sumup(coor,natoms,rcut,rmin,boxl): #,v12,v16,w12,w6):

	xcoor=coor[:,0] ; ycoor=coor[:,1] ; zcoor=coor[:,2]
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

def write_frame(tempmol,coor,scale_lj_to_protein,flag,filename):

	pi=numpy.pi
	
	npatoms=tempmol.natoms()
	comprotein=coor*scale_lj_to_protein
	k=1
	outfile1='temp.pdb' ; outfile2='traj.pdb'
	for i in xrange(len(coor)):
		tempmol.center(0)
        	ang=random.random()*1.0*pi; tempmol.rotate(0,'x',ang) ; tempmol.center(0)
        	ang=random.random()*1.0*pi; tempmol.rotate(0,'y',ang) ; tempmol.center(0)
        	ang=random.random()*1.0*pi; tempmol.rotate(0,'z',ang) ; tempmol.center(0)

        	tempmol.moveto(0,[comprotein[i][0],comprotein[i][1],comprotein[i][2]])
        	if(i<10):       segname='000'+str(i+1)
        	elif(i<100):    segname='00'+str(i+1)
        	elif(i<1000):   segname='0'+str(i+1)
        	elif(i<10000):   segname=str(i+1)
        	else: print 'bug in code NATOMS > 9999',3/0
        	segarray=[]  ; indexarray=[]
        	for s in xrange(npatoms):
                	segarray.append(segname)
                	indexarray.append(k)
                	k+=1
        	tempmol.setSegname(segarray)

        	if(i==0 and flag=='new'):
                	tempmol.writepdb(outfile1,0,'w')
        	else:
                	tempmol.writepdb(outfile1,0,'a')

	ost='cat temp.pdb | grep -v END > barf.pdb'
	os.system(ost)
	outend=open('end.txt','w')
	outend.write('END\n')
	outend.close()
	if(flag=='new'):
		ost='cat barf.pdb end.txt > '+filename
		os.system(ost)
		ost= 'rm -f barf.pdb temp.pdb end.txt'	
		os.system(ost)
	else:
		ost='cat '+outfile2+' barf.pdb end.txt > ttemp.pdb'
		os.system(ost)
		ost='mv ttemp.pdb '+outfile2
		os.system(ost)
		ost= 'rm -f barf.pdb temp.pdb end.txt'	
		os.system(ost)
	
	return 


def prot_to_lj():

	## this is not implemented yet !!!

	molecules=[] 
	for i in xrange(nmolecules):
		tmol=sasmol.SasMol(i)
		molecules.append(tmol)	
	
	return

def lj_to_prot(moleculepdb,start_type):

	ljradius = 0.5521 #  argon VDW radius / sigma argon
	sigmaAr=3.405

	tempmol=sasmol.SasMol(0)
	tempmol.readpdb(moleculepdb)

	tempmol.center(0)
	mm1=tempmol.calcminmax()

	dx=math.fabs(mm1[0][0]-mm1[1][0])
	dy=math.fabs(mm1[0][1]-mm1[1][1])
	dz=math.fabs(mm1[0][2]-mm1[1][2])

	maxdiamolecule = max(dx,dy,dz) ; print 'maximum molecular diameter = ',maxdiamolecule
	mindiamolecule = min(dx,dy,dz) ; print 'minimum molecular diameter = ',mindiamolecule

	print 'max diameter = ',maxdiamolecule
	print 'min diameter = ',mindiamolecule

	scale_lj_to_protein = sigmaAr*maxdiamolecule/(ljradius*2.0*sigmaAr)

	print 'scale_lj_to_protein = ',scale_lj_to_protein		

	if(start_type==1):
	
		coor,boxl=makefcc(2.0*ljradius)
			
	elif(start_type==3):

		coor,boxl=readxyz(coorfilename)

	else:
		print 'incorrect start up type indicated ',start_type
		print '# 1--> make fcc coords, 3--> restart from LJ coords'
		print 'STOPPING NOW\nSTOPPING NOW\nSTOPPING NOW\n\n',3/0

	flag='new' ; filename='traj0.pdb'
	write_frame(tempmol,coor,scale_lj_to_protein,flag,filename)
	del tempmol

	return coor,ljradius,boxl,scale_lj_to_protein


def nptmc(nmolecules,temperature,goalpressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename,coorfilename,iratio,vratio,masterpdb,moleculepdb,specificmolvol,start_type):

	
	pi=numpy.pi

	lj=sasmol.SasAtm(0)

	msystem=sasmol.SasAss(0)
	msystem.readpdb(masterpdb)

	molecule=sasmol.SasMol(0)
	molecule.readpdb(moleculepdb)
	
	#start_type: 1--> make fcc coords, 2--> restart from molecular coords(PDB), 3--> restart from LJ coords(XYZ)
	
	if(start_type==1):

		coor,radius,boxl,scale_lj_to_protein=lj_to_prot(moleculepdb,start_type)

	elif(start_type==2):

		print 'incorrect start up type indicated ',start_type
		print '# 1--> make fcc coords, 3--> restart from LJ coords'
		print 'STOPPING NOW\nSTOPPING NOW\nSTOPPING NOW\n\n',3/0

		prot_to_lj(msystem)
	
	elif(start_type==3):
	
		coor,radius,boxl,scale_lj_to_protein=lj_to_prot(moleculepdb,start_type)

	else:
		print 'incorrect start up type indicated ',start_type
		print '# 1--> make fcc coords, 3--> restart from LJ coords'
		print 'STOPPING NOW\nSTOPPING NOW\nSTOPPING NOW\n\n',3/0

	#coor,boxl,radius=setup(nmolecules,rho_start)	
	#coor,radius,boxl=readxyz(coorfilename)

	deltag = boxl/(2.0*nbins)
	gr = numpy.zeros(nbins,numpy.float)

	writexyz(coor,boxl,runfilename,'w')
	writexyz(coor/boxl,boxl,'test.xyz','w') 
	
	boxfile=open('box.dat','w')
	pressfile=open('press.dat','w')
	densityfile=open('density.dat','w')

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

	v12,v6,w12,w6 = sumup(coor,natoms,rcut,rmin,boxl) #,v12,v6,w12,w6)
 
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

	fr=1
	pdbfiles=[]

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
			thisi=i+1
			
			tgr = numpy.zeros(nbins,numpy.float)

			tgr,v12old,v6old,w12old,w6old=nenergy.lje(coor,rxiold,ryiold,rziold,thisi,boxl,deltag,nbins)

			rxinew = rxiold + ( 2.0 * random.random() - 1.0 ) * drmax
			ryinew = ryiold + ( 2.0 * random.random() - 1.0 ) * drmax
			rzinew = rziold + ( 2.0 * random.random() - 1.0 ) * drmax

			rxinew = rxinew-boxl*((rxinew*boxlinv).round())
			ryinew = ryinew-boxl*((ryinew*boxlinv).round())
			rzinew = rzinew-boxl*((rzinew*boxlinv).round())

			tgr = numpy.zeros(nbins,numpy.float)
			
			tgr,v12new,v6new,w12new,w6new=nenergy.lje(coor,rxinew,ryinew,rzinew,thisi,boxl,deltag,nbins)

			delv12 = v12new - v12old
			delv6 = v6new - v6old
			delw12 = w12new - w12old
			delw6 = w6new - w6old
			
			deltv = delv12 + delv6
			deltvb = beta * deltv

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
					gr = gr + tgr
				elif(math.exp(-deltvb) > random.random()):
					v12=v12+delv12
					v6 =v6+delv6
					w12=w12+delw12
					w6 =w6+delw6
					coor[i,0]=rxinew
					coor[i,1]=ryinew
					coor[i,2]=rzinew
					acatma=acatma+1.0
					gr = gr + tgr

	
			vlrc12 = 8.0*pi*density*natoms*sr9/9.0
			vlrc6 = -8.0*pi*density*natoms*sr3/3.0
			wlrc12 = 4.0*vlrc12
			wlrc6 =	2.0*vlrc6

			#vn = ( v12 + v6 ) / natoms
			vn = ( v12 + vlrc12 + v6 + vlrc6 ) / natoms
			#pressure = density * temperature + ( w12 + w6 ) / volume
			pressure = density * temperature + ( w12 + wlrc12 + w6 + wlrc6 ) / volume

			acm = acm + 1.0
			acv = acv + vn
			acp = acp + pressure
			acd = acd + density
	
			acvsq = acvsq + vn**2.0
			acpsq = acpsq + pressure**2.0
			acdsq = acdsq + density**2.0
			
		#END of for i in xrange(natoms)
		
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

		vlrc12 = 8.0*pi*density*natoms*sr9/9.0
		vlrc6 = -8.0*pi*density*natoms*sr3/3.0
		wlrc12 = 4.0*vlrc12
		wlrc6 =	2.0*vlrc6
		#vn = ( v12 + v6 ) /natoms
		vn = ( v12 +vlrc12 + v6 + vlrc6) /natoms
		#pressure = density * temperature + ( w12 + w6 ) / volume
		pressure = density * temperature + ( w12 + wlrc12 + w6 + wlrc6) / volume

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

			flag='new'
			filename='traj_'+str(fr)+'.pdb'
			pdbfiles.append(filename)
#			write_frame(molecule,coor,scale_lj_to_protein,flag,filename)
			fr+=1
		
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

	ost='catdcd -o traj'+str(goalpressure)+'.dcd -otype dcd -stype pdb -s traj0.pdb '
	st=' ' ; rmst='rm -f '	
	for f in xrange(len(pdbfiles)):
		st=st+' -pdb '+pdbfiles[f]+' '
		rmst=rmst+' '+pdbfiles[f]+' '
#	os.system(ost+st)
#	os.system(rmst)

	print '\n\npercent mcmoves accepted ',acatma*100.0/tm
	print '\n\npercent mcvols accepted ',acboxa*100.0/tv

	print 'goal pressure = ',goalpressure

	writefxyz(coor,'final_'+runfilename,boxl)

	boxfile.close()
	pressfile.close()

        if (acboxa>0):
                #rho=1.0/((boxsum/accvol)**3.0)
		rho = acd/(tm+tv)
                print 'rho = ',rho
        else:
                rho=1.0/(boxl**3.0)

        goutfile=open('gr.dat','w')
        for i in xrange(nbins):
                r=deltag*(i+0.5)
                vb=((i+1)**3-i**3)*deltag**3
                nid=(4./3)*numpy.pi*vb*rho
                #gr[i]=gr[i]/(accmove*natoms*nid)
                gr[i]=gr[i]/((acatma)*natoms*nid)
                #gr[i]=gr[i]/((tm)*natoms*nid)
                #goutfile.write('%i\t%f\t%f\n' % (i,r,gr[i]/nmolecules))
                goutfile.write('%i\t%f\t%f\n' % (i,r,gr[i]))

        goutfile.close()

	return

if __name__=="__main__":

	nmolecules,temperature,goalpressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename,coorfilename,iratio,vratio,masterpdb,moleculepdb,specificmolvol,start_type = parameters()

	nptmc(nmolecules,temperature,goalpressure,rho_start,nsteps,nbins,delta_move,delta_vol,runfilename,coorfilename,iratio,vratio,masterpdb,moleculepdb,specificmolvol,start_type)

	# for Argon, E/K_b = 119.8 K, sigma = 3.405E-10 m, M = 0.03994 kg/mol

	# length 	--> sigma
	# energy 	--> e
	# mass		--> m
	# time = sigma*sqrt(m/e)

	# potential energy == U* 	--> U/e
	# pressure == P* 		--> P * sigma^3 / e	
	# density == rho*		--> rho * sigma^3
	# temperature == T*		--> k_b * T / e

