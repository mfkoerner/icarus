#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 19:14:26 2017

@author: mitchell
"""
#######################################
# Requires .cif's in own folder within main one
# Requires baseDir and a list of compound names
# This is not refined code!!!! Just shows the concepts
#######################################


# modules
import os
import pymatgen as mg




# *************** #
# INPUT VARIABLES #
# *************** #


baseDir = 'DIRECTORY_WITH_.cif_FOLDER_HERE'
compoundlist = os.listdir('.cif')
compoundlist = set(compoundlist)
compoundlist = compoundlist - {'.DS_Store'}

# ********************** #
# END OF INPUT VARIABLES #
# ********************** #


os.chdir(baseDir)

# starting variables
compounds = {}

# grabbing cifs
for i in compoundlist:
	compounds[i] = mg.Structure.from_file(i)


os.chdir('..')

for i in compoundlist:
	dirname = i.split('_')[1].rstrip('.cif')
	compname = i.split('_')[0]
	os.mkdir(dirname)
	os.chdir(dirname)
	compounds[i].to(filename="POSCAR")
	with open(compname, 'w') as f:
		f.write('')
	os.chdir('..')

