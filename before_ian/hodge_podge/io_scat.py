import os,sys
import numpy
from matplotlib.pyplot import *

sys.path.append('/Users/curtisj/Desktop/smol')
import sasmol


bcoh={'H':-0.374, 'D':0.667, 'C':0.665, 'N':0.936, 'O':0.580, 'S':0.285}

binc={'H':2.5274, 'D':0.404, 'C':0, 'N':0, 'O':0, 'S':0}

sigcoh={'H':1.76, 'D':5.59, 'C':5.55, 'N':11.01, 'O':4.23, 'S':1.02}
siginc={'H':79.9, 'D':2.04, 'C':0.001, 'N':0.49, 'O':0.00, 'S':0.007}

#sigc=sigcoh_i/4PI
#sigi=siginc_i/4PI

filename='lys1.pdb'

mol1=sasmol.SasMol(0)

mol1.readpdb(filename)

com=mol1.calccom()

print 'com = ',com

element=mol1.element()

coh=0.0 ; inc=0.0
for atom in element:
	found=0
	for key,value in sigcoh.iteritems():
		if key==atom:
			coh=coh+value
			found=1
	for key,value in siginc.iteritems():
		if key==atom:
			inc=inc+value
			found=1
	if found==0:
		print 'need scattering data for ',atom

natoms=len(element)
natoms2=natoms*natoms

coherent = coh*natoms2/(4*numpy.pi)
incoherent = inc*natoms/(4*numpy.pi)

intensity = coherent+incoherent

print 'fraction coherent = ',coherent/intensity, ' percent = ',100.0*coherent/intensity
print 'fraction incoherent = ',incoherent/intensity , ' percent = ',100.0*incoherent/intensity

water_coherent=1.76*2+4.23 ; d2o_coherent=5.59*2+4.23
water_incoherent=79.9*2 ; d2o_incoherent=2.04*2

water_density=0.997 	# g/cm^3 @ 293K
d2o_density=1.104 	# g/cm^3 @ 293K

h2o=sasmol.SasMol(1)
h2o.readpdb('h2o.pdb')
mass_h2o=h2o.calcmass('element')

d2o=sasmol.SasMol(2)
d2o.readpdb('d2o.pdb')
mass_d2o=d2o.calcmass('element')

print 'mass h2o = ',mass_h2o
print 'mass d2o = ',mass_d2o

Na=6.02214E23

nh2o_1ml = water_density*Na/mass_h2o
nd2o_1ml = d2o_density*Na/mass_d2o

mass_mol1=mol1.calcmass('element')

#conc =[ x for x in xrange(1,5*150,1)]

nmol1=[] ; conc=[]
for x in xrange(1,5*150,1):
	val=x*0.2
	conc.append(val)
	nmol1.append(val*Na/(1000*mass_mol1))

na_nmol1=numpy.array(nmol1)

water_scat_coherent=nh2o_1ml*water_coherent
water_scat_incoherent=nh2o_1ml*water_incoherent
d2o_scat_coherent=nd2o_1ml*d2o_coherent
d2o_scat_incoherent=nd2o_1ml*d2o_incoherent

mol1_scat_coherent=na_nmol1*coherent
mol1_scat_incoherent=na_nmol1*incoherent

total_coherent_in_h2o=water_scat_coherent+mol1_scat_coherent
total_incoherent_in_h2o=water_scat_incoherent+mol1_scat_incoherent

total_coherent_in_d2o=d2o_scat_coherent+mol1_scat_coherent
total_incoherent_in_d2o=d2o_scat_incoherent+mol1_scat_incoherent

#fig=figure()
#ax = fig.add_subplot(111)

#ax.plot(conc,total_coherent_in_h2o/total_incoherent_in_h2o,label='h2o')
#ax.plot(conc,total_coherent_in_d2o/total_incoherent_in_d2o,label='d2o')

plot(conc,total_coherent_in_h2o/total_incoherent_in_h2o,label='protein in h2o')
plot(conc,total_coherent_in_d2o/total_incoherent_in_d2o,label='protein in d2o')

#ax.legend()
#legend()
ylabel('Coherent Scattering / Incoherent Scattering: "Signal to Noise"')
xlabel('Protein Concentration (mg/mL)')
title('Lysozyme Scattering at I(q=0)')

treh=sasmol.SasMol(3)
treh.readpdb('treh.pdb')
mass_treh=treh.calcmass('element')
treh_element=treh.element()

na_treh_40percent=0.4*Na/mass_treh

tcoh=0.0 ; tinc=0.0
for atom in treh_element:
	found=0
	for key,value in sigcoh.iteritems():
		if key==atom:
			tcoh=tcoh+value
			found=1
	for key,value in siginc.iteritems():
		if key==atom:
			tinc=tinc+value
			found=1
	if found==0:
		print 'need scattering data for ',atom

treh_scat_coherent=na_treh_40percent*tcoh ; 
treh_scat_incoherent=na_treh_40percent*tinc ; 

total_coherent_in_h2o=water_scat_coherent+mol1_scat_coherent+treh_scat_coherent
total_incoherent_in_h2o=water_scat_incoherent+mol1_scat_incoherent+treh_scat_incoherent

total_coherent_in_d2o=d2o_scat_coherent+mol1_scat_coherent+treh_scat_coherent
total_incoherent_in_d2o=d2o_scat_incoherent+mol1_scat_incoherent+treh_scat_incoherent

plot(conc,total_coherent_in_h2o/total_incoherent_in_h2o,label='protein, 40% trehalose in h2o')
plot(conc,total_coherent_in_d2o/total_incoherent_in_d2o,label='protein, 40% trehalose in d2o')

legend(loc='lower right')

draw()

show()


