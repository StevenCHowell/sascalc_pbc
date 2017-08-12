'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os,sys,locale,string
import sassie.sasmol.sasmol as sasmol
import sassie.interface.input_filter as input_filter


### OPEN:	HAS NOT BEEN IMPLEMENTED OR TESTED !!!!
### OPEN:	HAS NOT BEEN IMPLEMENTED OR TESTED !!!!
### OPEN:	HAS NOT BEEN IMPLEMENTED OR TESTED !!!!

def check_input(variables):
      
	runname					= variables['runname'][0]
	initial_pdb 				= variables['initial_pdb'][0]
	scattering_data_file 			= variables['scattering_data_file'][0]
	start_type				= variables['start_type'][0]
	number_of_molecules			= variables['number_of_molecules'][0]
	number_of_steps				= variables['number_of_steps'][0]
	number_of_cycles			= variables['number_of_cycles'][0]
	delta_translation			= variables['delta_translation'][0]
	delta_rotation				= variables['delta_rotation'][0]
	delta_volume				= variables['delta_volume'][0]
	temperature				= variables['temperature'][0]
	pressure				= variables['pressure'][0]
	initial_rho				= variables['initial_rho'][0]
	internal_dof				= variables['internal_dof'][0]
	number_of_flexible_regions		= variables['number_of_flexible_regions'][0]
	delta_theta				= variables['delta_theta'][0]
	low_residues				= variables['low_residues'][0]
	number_of_contiguous_residues		= variables['number_of_contiguous_resiudes'][0]
	align_low				= variables['align_low'][0]
	align_high				= variables['align_high'][0]
	translation_ratio			= variables['translation_ratio'][0]
	rotation_ratio				= variables['rotation_ratio'][0]
	volume_ratio				= variables['volume_ratio'][0]
	dtheta_ratio				= variables['dtheta_ratio'][0]


###	OPEN:	temporary hardwire for path

	path = './'

	error=[]
        ev,rv,wv=input_filter.check_permissions(path)
        if(not ev or not rv or not wv):
		error.append('permission error in input file path '+path+'  [code = '+str(ev)+str(rv)+str(wv)+']')
		if(ev==False):
			error.append('path does not exist')
		elif(rv==False):
			error.append('read permission not allowed')
		elif(wv==False):
			error.append('write permission not allowed')
		return error

	pdbfile=path+'/'+pdbfile
     	ev,value=input_filter.check_pdb_dcd(initial_pdb,'pdb')

	if(ev == 0):
        	error.append('input pdb file, '+pdbfile[3:]+', does not exist')
		return error
	if(value == 0):
        	error.append( 'input pdb file, '+pdbfile[3:]+', is not a valid pdb file')
		return error

	if(number_of_steps < 1):
        	error.append( 'trials = '+str(trials)+'?')
		return error
	elif(temperature < 0):
        	error.append( 'use a positive temperature, temperature = '+str(temp))
		return error
#	elif(moltype != 'protein' and moltype != 'rna'):
#		error.append('only protein and rna backbone dihedral move sets are defined, you entered : '+str(moltype))
#		return error	
#	elif(cutoff < 1.0):
#        	error.append( 'use a larger cutoff value, cutoff = '+str(cutoff))
#		return error
		return error
	elif(len(delta_theta) != number_of_flexible_regions):
        	error.append( 'the number of dtheta values does not match the number of ranges, dtheta = '+str(delta_theta)+' numranges = '+str(number_of_flexible_regions))
		return error
	elif(len(low_residues) != number_of_flexible_regions):
        	error.append( 'the number of low residue values does not match the number of ranges, lowres = '+str(low_residues)+' numranges = '+str(number_of_flexible_regions))
		return error
	elif(len(number_of_contiguous_residues) != number_of_flexible_regions):
        	error.append( 'the number of contiguous residues does not match the number of ranges, lowres = '+str(number_of_contiguous_residues)+' numranges = '+str(number_of_flexible_regions))
		return error

	for th in delta_theta:
		if(th > 180.0):
			delta_theta[th]=180.0
		elif(th < 0.0):
			delta_theta[th]=0.0

	locvariables=['resid']
	value,result=input_filter.get_pdb_stats(pdbfile,locvariables)
	resid=map(int,result[0])

#	if(resid[0] != 1):
#		error.append('amino acid residues in starting pdbfile '+pdbfile+' must start at resid = 1 : '+str(resid[0]))
#		return error

	number_aa = resid[-1] - resid[0]+1

	for i in xrange(resid[0],number_aa):
		ti=i+1
		if ti not in resid:
			error.append('amino acid '+str(ti)+' is missing from pdbfile'+pdbfile)
			print 'amino acid '+str(ti)+' is missing from pdbfile'+pdbfile
			return error

	for i in xrange(number_of_flexible_regions-1):
		if(low_residues[i]< resid[0]):
			error.append('low residue is lower than the n-terminal amino acid number, reslow = '+str(low_residues))
			return error	
		elif(low_residues[i]>low_residues[i+1]):
			error.append('low residue values must increase from low to high, reslow = '+str(low_residues))
			return error	
		elif(low_residues[i]+number_of_flexible_regions[i] > low_residues[i+1]):
			error.append('low residue values plus number contiguous overlap, reslow = '+str(low_residues)+' numcont = '+str(number_of_flexible_regions))
			return error	
		elif(low_residues[i]+number_of_flexible_regions[i] > number_aa+1):
			error.append('your low residue plus number contiguous exceed the number of amino acids ('+str(number_aa)+'), reslow = '+str(low_residues)+' numcont = '+str(number_of_flexible_regions))
			return error	
	if(reslow[-1]+number_of_flexible_regions[-1] > number_aa+1):
			error.append('your low residue plus number contiguous exceed the number of amino acids ('+str(number_aa)+'), reslow = '+str(reslow)+' numcont = '+str(number_of_flexible_regions))
			return error	

        if(align_low not in resid):
                error.append('input pdb file, '+str(pdbfile)+' does not have low residue amino acid, '+str(align_low)+', range = '+str(resid[0])+' : '+str(resid[-1]))
                return error
        elif(align_high not in resid):
                error.append('input pdb file, '+str(pdbfile)+' does not have high residue amino acid, '+str(align_high)+', range = '+str(resid[0])+' : '+str(resid[-1]))
                return error
        elif(align_high - align_low < 3):
                error.append('alignment basis is too small (less than 3 points) or low residue > high residue')
                return error

	# check residue topology and atom ordering
		
	m1=sasmol.SasMol(0)
	#error = m1.read_pdb(pdbfile,fastread=True,saspdbrx_topology=True)
	if(len(error)>0):
		error.append(result)
		return error

	return error



