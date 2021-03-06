#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 23:57:28 2017

@author: mitchell
"""


jobText = """#!/bin/sh
############################################################################## 
# First we have some directives to the queueing system. Must be in 
# the first block of comment lines. 
#
##PBS -q largemem
#PBS -l nodes={0}:ppn={1}
#PBS -l walltime={2}:00:00:00
# 
# Make sure that we are in the same subdirectory as where the qsub command 
# is issued. 
# 
cd $PBS_O_WORKDIR 
#
#  Determine the nodes, num process, etc.
#  cat $PBS_NODEFILE > nodes
#  oddly,, this version puts all nodes on one line...
#  mpirun wants separate lines though.
while read machine
do
echo $machine
done < $PBS_NODEFILE > nodes
# Get number of nodes allocated
NO_OF_NODES=`cat $PBS_NODEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{{print $1}}'`
NODE_LIST=`cat $PBS_NODEFILE`
NUM_CORES=`cat $PBS_NODEFILE | wc -w`
#
# Our list of nodes...                      
echo $NODE_LIST
# 
# 
# Run the executable. *DO NOT PUT* a '&' at the end - it will not 
# work. 
#
ulimit -s unlimited 
#VASP.5.4.4 non-collinear build (has SOC)
/sw/mrl/mvapich/bin/mpirun -n $NUM_CORES -machinefile nodes /home/burak/vasp.5.4.4/bin/vasp_ncl > out                                         
#
# End of script-file. 
# 
############################################################################### 
"""

jobText533 = """#!/bin/sh
############################################################################## 
# First we have some directives to the queueing system. Must be in 
# the first block of comment lines. 
#
##PBS -q largemem
#PBS -l nodes={0}:ppn={1}
#PBS -l walltime={2}:00:00:00
# 
# Make sure that we are in the same subdirectory as where the qsub command 
# is issued. 
# 
cd $PBS_O_WORKDIR 
#
#  Determine the nodes, num process, etc.
#  cat $PBS_NODEFILE > nodes
#  oddly,, this version puts all nodes on one line...
#  mpirun wants separate lines though.
while read machine
do
echo $machine
done < $PBS_NODEFILE > nodes
# Get number of nodes allocated
NO_OF_NODES=`cat $PBS_NODEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{{print $1}}'`
NODE_LIST=`cat $PBS_NODEFILE`
NUM_CORES=`cat $PBS_NODEFILE | wc -w`
#
# Our list of nodes...                      
echo $NODE_LIST
# 
# 
# Run the executable. *DO NOT PUT* a '&' at the end - it will not 
# work. 
#
ulimit -s unlimited 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/openmpi/lib
#VASP.5.3.3 Jakoah Build (can handle multiple nodes)
/opt/openmpi/bin/mpirun -np $NUM_CORES -machinefile nodes /home/jdb/bin/vasp533_openmpi-jakoah > out                                         
#
# End of script-file. 
# 
############################################################################### 
"""


import icarus.config as config
import os
import shutil
from subprocess import call
import numpy as np
from collections import OrderedDict
import sys




notes = '''run from code most recently updated on Friday Feb. 23, 2018 (autoVASP.py in https://github.com/mfkoerner/icarus)
The cutoff energies as of now are decided to be 700 eV if contains C N O F period 2 anions, and 500 eV otherwise
All K-POINT meshes are automatic mode density of 30 except absorption which is deisity of 60 and line mode for bands of course'''
a700 = {'C', 'N', 'O', 'F'}             #If we have any of the following elements, encut is automatically 700 instead of 500

recommended_PAW={ 'H':'', 'He':'', 'Li':'_sv', 'Be':'', 'B':'', 'C':'', 'N':'',
                 'O':'', 'F':'', 'Ne':'', 'Na':'_pv', 'Mg':'', 'Al':'', 'Si':''
                 , 'P':'', 'S':'', 'Cl':'', 'Ar':'', 'K':'_sv', 'Ca':'_sv', 
                 'Sc':'_sv', 'Ti':'_sv', 'V':'_sv', 'Cr':'_pv', 'Mn':'_pv', 
                 'Fe':'', 'Co':'', 'Ni':'', 'Cu':'', 'Zn':'', 'Ga':'_d', 'Ge'
                 :'_d', 'As':'', 'Se':'', 'Br':'', 'Kr':'', 'Rb':'_sv', 'Sr'
                 :'_sv', 'Y':'_sv', 'Zr':'_sv', 'Nb':'_sv', 'Mo':'_sv', 'Tc'
                 :'_pv', 'Ru':'_pv', 'Rh':'_pv', 'Pd':'', 'Ag':'', 'Cd':'', 
                 'In':'_d', 'Sn':'_d', 'Sb':'', 'Te':'', 'I':'', 'Xe':'', 'Cs'
                 :'_sv', 'Ba':'_sv', 'La':'', 'Ce':'', 'Pr':'_3', 'Nd':'_3', 
                 'Pm':'_3', 'Sm':'_3', 'Eu':'_2', 'Gd':'_3', 'Tb':'_3', 'Dy'
                 :'_3', 'Ho':'_3', 'Er':'_3', 'Tm':'_3', 'Yb':'_2', 'Lu':'_3', 
                 'Hf':'_pv', 'Ta':'_pv', 'W':'_sv', 'Re':'', 'Os':'', 'Ir':'', 
                 'Pt':'', 'Au':'', 'Hg':'', 'Tl':'_d', 'Pb':'_d', 'Bi':'_d', 
                 'Po':'_d', 'At':'_d', 'Rn':'', 'Fr':'_sv', 'Ra':'_sv', 'Ac':''
                 , 'Th':'', 'Pa':'', 'U':'', 'Np':'', 'Pu':'', 'Am':'', 'Cm':''}

filestart = """automatic generation - MP scheme
 0
Monkhorst-Pack
"""

needforjob = {'INCAR','KPOINTS','POSCAR','POTCAR'}

def write_KPOINTS(mesh, mode='auto'):
    if mode[0]=='a' or mode[0]=='A':
        assert(isinstance(mesh, int))
        outlines = ['fully automatic', ' 0', 'Auto', '     {}'.format(mesh)]
    if mode[0]=='M' or mode[0]=='m':
        assert(len(mesh) == 3)
        outlines = ['grid automatic', ' 0', 'Monkhorst-Pack', '  {} {} {}'.format(*mesh)]
    outlines = [i + '\n' for i in outlines]
    with open('KPOINTS', 'w') as f:
        f.writelines(outlines)
#INCAR files
def write_INCAR(indict, filepath = 'INCAR'):
    outlines = []
    if 'System' in indict.keys():
        outlines.append('System = {}\n'.format(indict.pop('System')))
    else:
        outlines.append('System = unknown\n')
    for key in indict.keys():
        if isinstance(indict[key], tuple) or isinstance(indict[key], list):
            outlines.append('{:10}= {:15}# {}\n'.format(key, str(indict[key][0]), indict[key][1]))
        else:
            outlines.append('{:10}= {:15}\n'.format(key, str(indict[key])))
    with open(filepath, 'w') as f:
        f.writelines(outlines)
def apply_all(indict, compound, symprec = config.VASP_symprec):
    indict['GGA'] = ('PE', 'Perdew-Burke-Ernzerhof')
    indict['LREAL'] = ('.FALSE.', 'Reciprocal space projection')
    indict['ISMEAR'] = (-5, '-5 -> tet, 0 -> gaussian')
    indict['EDIFF'] = (0.00001, 'stopping criterion for electronic self-consistency')
    indict['NSW'] = (0, 'no strucural relaxation')
    indict['SYMPREC'] = (symprec, 'symmetry requirements')
def parallel(indict):
    indict['NCORE'] = 4
    indict['NSIM'] = 4
def start_incar(compound, section, params_path = '../params', Apply_all = True, Parallel = False,
    Band = False, Soc = False, Continue_run = False, Hse = False, Absorb = None, DOS = False):
    indict = OrderedDict()
    indict['System'] = '{} {}'.format(compound, section)
    if Apply_all:
        apply_all(indict, compound)
    if Band:
        band(indict)
    if Soc:
        soc(indict)
    if Hse:
        hse(indict)
    if Continue_run:
        continue_run(indict)
    if Parallel:
        parallel(indict)
    if not Absorb is None:
        assert(isinstance(Absorb, int)), "Absorb needs to be and integer for start_incar"
        absorb(indict, Absorb)
    if DOS:
        dos(indict)
    permaread(indict, path = params_path)
    return(indict)
def dos(indict):
    indict['NEDOS'] = (4001, 'number of energy levels calculated')
    indict['EMIN'] =  (-25, 'minimum energy for DOS calcs')
    indict['EMAX'] =  (15, 'max energy for DOS calcs')
    indict['SIGMA'] = (0.001, 'smearing of electrons')
    indict['EDIFF'] = (1e-7, 'energy convergence criterion')
    indict['ISIF'] =  (0, 'symmetry off')
    indict['ICHARG']= (11, 'non-self-consistent from CHGCAR')
    indict['LORBIT']= (11, 'write DOSCAR and PROCAR with all quantum numbers')
    indict['ISMEAR']= (0, '-5 -> tet, 0 -> gaussian')
def band(indict):
    indict['ICHARG'] = (11, 'non-self-consistent from CHGCAR')
    indict['LORBIT'] = (11, 'Write DOSCAR and PROCAR, quick scheme, "l" only')
    indict['ISMEAR'] = (0, '-5 -> tet, 0 -> gaussian')
    indict['ISYM']   = (0, "don't enforce symmetry")
def soc(indict):
    indict['LSORBIT'] = ('.TRUE.', 'turn on spin-orbit coupling')
    indict['LNONCOLLINEAR'] = '.TRUE.'
def continue_run(indict):
    indict['ICHARG'] = (1, 'start from previous CHGCAR')
def hse(indict):
    indict['ISMEAR'] = (0, '-5 -> tet, 0 -> gaussian')
    indict['LHFCALC'] = ('.TRUE.', 'Hartree-Fock exchange')
    indict['HFSCREEN'] = (0.2, 'screening distance for HF exchange')
    indict['ALGO'] = ('All', 'iterative solver parameter')
    indict['TIME'] = (0.4, 'iterative solver parameter')
    indict['PRECFOCK'] = ('Fast', 'Normal gives better quality for high cost')
def absorb(indict, NBANDS):
    """ If you can get away with it, use tetrahedral ISMEAR """
    indict['NBANDS'] = (NBANDS, 'Number of bands calculated')
    indict['NEDOS'] = (2000, 'number of DOS energies calculated (big -> fine features)')
    indict['LOPTICS'] = ('.TRUE.', 'Turns on optical absorbtion calculation')
    indict['ICHARG'] = (11, 'non-self-consistent from CHGCAR')
    indict['LWAVE'] = ('.FALSE.', 'saving space by not writing WAVECAR')
    indict['SIGMA'] = (0.002, 'electron smearing: default is 0.2 eV')


def write_job(NODES = config.VASP_nodes, PPN = config.VASP_ppn, WALLTIME = config.VASP_walltime, version = 'auto'):
    """write the job text

        Inputs:
            NODES:      how many nodes to run on
            PPN:        Processers per node. Braid can do 16, 20, 24 or less
            WALLTIME:   Requested time for job
            version:    auto bases off of versionXXX file in parent directory
                        '5.3.3' gives 5.3.3 and '5.4.4' gives 5.4.4
                        no other options at present
    """
    if version == 'auto':
        if 'version533' in os.listdir('../'):
            version = '5.3.3'
        else:
            version = '5.4.4'
    with open('job', 'w') as f:
        if version == '5.4.4':
            f.write(jobText.format(NODES,PPN,WALLTIME))
        elif version == '5.3.3':
            f.write(jobText533.format(NODES,PPN,WALLTIME))
        else:
            raise ValueError('invalid version for write_job')
def subjob():
    assert(needforjob.issubset(set(os.listdir())))
    call( ['qsub', 'job'] )
    with open('autonotes', 'w') as f:
        f.write(notes)
def getNBANDS(path = 'OUTCAR'):
    with open(path, 'r') as f:
        inlines = f.readlines()
    implines = [i for i in inlines if 'NBANDS' in i]
    assert(len(implines) == 1)
    line = implines[0]
    return(int(line[line.index('NBANDS')+7:].strip(' ')))
def permawrite(indict, path='params'):
    outlines = []
    for key in indict.keys():
        if isinstance(indict[key], tuple) or isinstance(indict[key], list):
            outlines.append('{:10}= {:15}# {}\n'.format(key, str(indict[key][0]), indict[key][1]))
        else:
            outlines.append('{:10}= {:15}\n'.format(key, str(indict[key])))
    print('writing the following to params:\n' + ''.join(outlines))
    with open(path, 'a') as f:
        f.writelines(outlines)
def permaread(indict, path='params'):
    with open(path, 'r') as f:
        inlines = [i.rstrip('\n') for i in f.readlines()]
    for line in inlines:
        key = line.split('=')[0].strip(' ')
        value = line.split('=')[1].strip(' ')
        if '#' in value:
            value = value[:value.index('#')].rstrip(' ')
        try:
            value = eval(value)
        except:
            pass
        indict[key] = value
def autoencut(path = 'POSCAR'):
    with open(path, 'r') as f:
        inlines = [i.rstrip('\n') for i in f.readlines()]
    atoms = set(inlines[5].split())
    if len(atoms.intersection(a700)) > 0:
        return(700)
    else:
        return(500)
def autocompound():
    return(os.getcwd().split('/')[-1])






###############################################################################
###############################################################################
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#--------------------------------Running Stuff--------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
###############################################################################
###############################################################################


def runDIRECTORY(dirs = config.VASP_directories):
    """ directory setup """
    for i in dirs:
        if i not in os.listdir():
            os.mkdir(i)
    if 'params' not in os.listdir():
        with open('params', 'w') as f:
            f.write('')
def runPOTCAR():
    """
    To begin you need POSCAR, job, hybridjob (put this in folder that holds the compounds)
    Also will need KPOINTS for bandstructure
    Make sure job has 4 nodes with 4 ppn or change around INCAR above
    """
    with open('POSCAR', 'r') as f:
        inlines = [i.rstrip('\n') for i in f.readlines()]
    matline = inlines[5].split()
    POSnames = [i + recommended_PAW[i] for i in matline]
    with open('POTCAR', 'w') as f:
        call(['cat'] + 
             ['/home/jdb/vasp_support/potentials/potpaw_PBE.54/{}/POTCAR'.format(i)
             for i in POSnames], stdout = f)
    write_job()


class Run(object):
    """Run class to initiate vasp runs

    Inputs:
        compound:   name of material (defaults to directory name)
        NODES:      number of nodes to run on
        PPN:        processers per node
        WALLTIME:   time that the run goes for
    """
    def __init__(self, compound = None, NODES = None, PPN = None, WALLTIME = None, version = 'auto', assume_memerror = True):
        if compound is None:
            compound = autocompound()
        self.compound = compound

        if NODES is None:
            self.NODES = config.VASP_nodes
        else:
            self.NODES = NODES

        if PPN is None:
            self.PPN = config.VASP_ppn
        else:
            self.PPN = PPN

        if WALLTIME is None:
            self.WALLTIME = config.VASP_walltime
        else:
            self.WALLTIME = WALLTIME

        if version == 'auto':
            if 'version533' in os.listdir('.'):
                self.version = '5.3.3'
            else:
                self.version = '5.4.4'
        else:
            self.version = version

        if assume_memerror and self.version == '5.3.3':
            if NODES is None:
                self.NODES = 4
            if WALLTIME is None:
                WALLTIME = 2 * config.VASP_walltime

    def STATIC(self, kpoints = 30, encut = 'auto', submit = True):
        """Runs a static job

        Inputs
            submit:
                True is default and whole function runs
                False runs everything but submit it
                'only' assumes it has already been run with submit = False
        """
        os.chdir('static')
        if not submit == 'only':
            shutil.copy('../POSCAR', '.')
            shutil.copy('../POTCAR', '.')
            write_KPOINTS(kpoints)
            tester = {}
            permaread(tester, '../params')
            if not 'ENCUT' in tester.keys():
                if encut == 'auto':
                    permawrite({'ENCUT': autoencut(path='POSCAR')}, path = '../params')
                else:
                    permawrite({'ENCUT': encut}, path = '../params')
            indict = start_incar(self.compound, 'static', Parallel = False)
            write_INCAR(indict)
            write_job(NODES = self.NODES, PPN = self.PPN, WALLTIME = self.WALLTIME)
        if submit:
            subjob()
        os.chdir('../')

    def BAND(self):
        """
        band portion (REMEMBER TO PUT KPOINTS IN BAND FOLDER)

        check to make sure static runs went alright
        remember to use wavecheck program to check runs
        
        make sure bands KPOINTS file is already in the band directory
        """        
        os.chdir('band')
        indict = start_incar(self.compound, 'band', Parallel = False, Band = True)
        write_INCAR(indict)
        write_job(NODES = self.NODES, PPN = self.PPN, WALLTIME = self.WALLTIME)
        call( ['cp', '../static/POTCAR', '../static/POSCAR',
               '../static/CHGCAR', '.'] )
        subjob()
        os.chdir('../')
    def SOC(self):
        """ SOC portion (Run at same time as bands) """
        os.chdir('pbesoc')
        write_job(NODES = self.NODES, PPN = self.PPN, WALLTIME = self.WALLTIME)
        call( ['cp', '../static/POTCAR', '../static/POSCAR',
              '../static/CHGCAR', '../static/KPOINTS', '.'] )
        indict = start_incar(self.compound, 'pbesoc calculation', Parallel = False, Soc = True, Continue_run = True)
        write_INCAR(indict)
        subjob()
        os.chdir('../')
    def SOCBAND(self):
        """
        SOC band portion 

        Make sure that pbesoc has finished running. It's ok if band did not finish
        """
        os.chdir('SOCband')
        write_job(NODES = self.NODES, PPN = self.PPN, WALLTIME = self.WALLTIME)
        call( ['cp', '../static/POTCAR', '../pbesoc/POSCAR',
               '../pbesoc/CHGCAR', '../band/KPOINTS', '.'] )
        indict = start_incar(self.compound, 'Spin Orbit Coupling Band Structure', Parallel = False, Soc = True, Band = True)
        write_INCAR(indict)
        subjob()
        os.chdir('../')
    def HSESOC(self):
        """ HSESOC section (run at same time as bands)"""
        os.chdir('HSESOC')
        call( ['cp', '../static/POTCAR', '../static/POSCAR', '../static/CHGCAR',
               '../static/KPOINTS', '.'] )
        write_job(NODES = self.NODES, PPN = self.PPN, WALLTIME = self.WALLTIME)
        indict = start_incar(self.compound, 'HSESOC static', Parallel = False, Soc = True, Continue_run = True, Hse = True)
        write_INCAR(indict)
        subjob()
        os.chdir('../')
    def SHBAND(self):
        """ SHband presetup section
        (SHband is pretty manual for KPOINTS and qsub)
        make sure that HSESOC is finished to run this """
        os.chdir('SHband')
        call( ['cp', '../POTCAR', '../HSESOC/POSCAR', '../HSESOC/CHGCAR', '.'] )
        write_job(NODES = self.NODES, PPN = self.PPN, WALLTIME = self.WALLTIME)
        indict = start_incar(self.compound, 'SOC + HSE06 band structure calculation', Parallel = False, Soc = True, Hse = True, Band = True)
        print('Remember to prep KPOINTS based on graphs!')
        os.chdir('../')
    def ABSORB(self, KPOINTS = 'auto', NBANDS = 'auto'):
        """ Absorbtion spectrum calculation """
        os.chdir('absorb')
        write_job(NODES = self.NODES, PPN = self.PPN, WALLTIME = self.WALLTIME)
        call( ['cp', '../static/POTCAR', '../static/POSCAR', '../static/CHGCAR', '.'] )
        if KPOINTS == 'auto':
            with open('../static/KPOINTS', 'r') as f:
                inlines = f.readlines()
            density = int(inlines[3].strip(' '))
            assert ( inlines[2][0] in {'A','a'} ),  "Needs to be automatic scheme in kpoints"
            KPOINTS = 3*density
        write_KPOINTS(KPOINTS)
        if NBANDS == 'auto':
            NBANDS = 2*getNBANDS('../static/OUTCAR')
        indict = start_incar(self.compound, 'Absorption run', Parallel = False, Absorb = NBANDS)
        write_INCAR(indict)
        subjob()
        os.chdir('../')
    def DOS(self):
        """ DOS Calculation """
        os.chdir('DOS')
        write_job(NODES = self.NODES, PPN = self.PPN, WALLTIME = self.WALLTIME)
        call( [ 'cp', '../static/POTCAR', '../static/POSCAR', '../static/CHGCAR', '.'])
        with open('../static/KPOINTS', 'r') as f:
            inlines = f.readlines()
        assert( inlines[2][0] in {'A','a'} ) #Needs to be automatic scheme in kpoints
        density = int(inlines[3].strip(' '))
        write_KPOINTS(3*density)
        indict = start_incar(self.compound, 'DOS run', Parallel = False, DOS = True)
        write_INCAR(indict)
        subjob()
        os.chdir('../')


###############################
# Interacting with status.txt #
###############################

class Status(object):
    """Status is used to interact with the status.txt file
    It needs a bigger docstring"""
    def __init__(self, filepath = config.statuspath):
        self.filepath = filepath
        self.dictform = self.statusread()
        self.all_mpids = set(self.dictform.keys())
        self.subspace = self.get_mpids('comments', 'subspacematrix')
        self.killed = self.get_mpids('comments', 'Killed (not job.e)')
        self.metallic = self.get_mpids('comments', 'metallic')
        self.timeout = self.get_mpids('comments', 'timeout')
        self.mem = self.get_mpids('comments', 'Insufficient Memory')

    def statuswrite(self, filepath = None):
        """Writes a status dictionary to the official status file

        Inputs
            status: dictionary of {mpid1: {trait: value}, mpid2: {trait: value}}
            filepath: str filepath to status.txt file

        Outputs
            overwrites the file status.txt with the new status data

        Notes on static and absorb codes
            Both static and absorb have the same code
                0  : not done
                1  : done
                2  : in progress
                3  : investigate
                -1 : abort (permanently)
        """
        if filepath is None:
            filepath = self.filepath
        outlines = ['{:<12}{:<7}{:<7}{:<12}{}'.format(i['mpid'], i['static'], i['absorb'], i['origin'], i['comments'])
                    for i in self.dictform.values()]
        header = '{:<12}{:<7}{:<7}{:<12}{}\n'.format('mpid', 'static', 'absorb', 'origin', 'comments')
        outlines = [i + '\n' for i in outlines]
        outlines.sort()
        outlines = [header] + outlines
        with open(filepath, 'w') as f:
            f.writelines(outlines)

    def statusread(self, filepath = None):
        """Reads the status.txt file to a dictionary

        Inputs
            filepath: str filepath to status.txt file

        Outputs
            status: dictionary of {mpid1: {trait: value}, mpid2: {trait: value}}
        """
        if filepath is None:
            filepath = self.filepath
        slens = [12, 19, 26, 38]    #needs to match the sizes in outlines and header for statuswrite
        with open(filepath, 'r') as f:
            inlines = [i.rstrip('\n') for i in f.readlines()]
        header = inlines.pop(0)
        headerlist = header.split()
        inlist = [[str(i[:slens[0]]).strip(), int(i[slens[0]:slens[1]]), 
                   int(i[slens[1]:slens[2]]), str(i[slens[2]:slens[3]]).strip(),
                   str(i[slens[3]:]).strip()] for i in inlines]
        status = {i[0]: {headerlist[j]: i[j] for j in range(5)} for i in inlist}
        return(status)

    def get_mpids(self, attribute, value, invert = False, source = None):
        """gets a list of mpids with a given attribute and status

        Inputs:
            status:     dictionary of {mpid1: {trait: value}, mpid2: {trait: value}}
            attribute:  one of "absorb", "static", "mpid", "origin", "comments"
            value:      what you want your attribute to be equal to
                        May input a list to check for multiple values
            invert:     default False, if True inverts selection
            source:     dictionary to search from, defaults to self.dictform

        Outputs:
            mpids:      list of ids which match the above criteria
        """
        if source == None:
            source = self.dictform
        all_mpids = set(source.keys())
        if type(value) in {list, set, tuple}:
            values = value
        else:
            values = [value]
        mpidslist = []
        for val in values:
            mpidslist.append({i for i in all_mpids if source[i][attribute] == val})
        mpids = mpidslist[0]
        if len(mpidslist) > 1:
            for i in range(len(mpidslist) - 1):
                mpids = mpids.union(mpidslist[i + 1])
        if invert:
            mpids = all_mpids - mpids
        return(mpids)

    def set_from_file(self, filepath):
        """takes lines from filepath as input and returns a set of those strings"""
        with open(filepath, 'r') as f:
            values = {i.rstrip('\n') for i in f.readlines()}
        return(values)

    def set_to_file(self, yourset, filepath):
        """writes a set to a file separated by lines"""
        outlines = [i + '\n' for i in yourset]
        with open(filepath, 'w') as f:
            f.writelines(outlines)

    def attribute_to_value(self, mpids, attribute, value):
        """sets value of attribute to value for all mpids given

        Inputs:
            mpids:      iterable of mpid strings to set
            attribute:  one of "absorb", "static", "mpid", "origin", "comments"
            value:      what you want your attribute to be set to for all of mpids

        Outputs:
            modifies self.dictform (base data) to include this update
        """
        if type(mpids) == str:          #handles single mpid passed in as a string
            mpids = [mpids]
        for mpid in mpids:
            self.dictform[mpid][attribute] = value
    def toss(self, mpids, comment):
        """sets static and absorb to -1 and then sets the comment on all mpids given

        Inputs:
            mpids:      iterable of mpid strings to set
            comment:    the comment you wnat to associate with this dead compound

        Outputs:
            modifies self.dictform (base data) to include this update
            specifically, sets static and absorb to -1
            sets comments to comment
        """
        if type(mpids) == str:          #handles single mpid passed in as a string
            mpids = [mpids]
        for mpid in mpids:
            self.dictform[mpid]['static'] = -1
            self.dictform[mpid]['absorb'] = -1
            self.dictform[mpid]['comments'] = comment
    def update_from_playground(self, attribute):
        """Very dangerous function to be used after finishing a batch of calculations

        Inputs:
            attribute: either "absorb" or "static"

        Filesystem inputs:
            done:       currentdirectory/done needs to be the mpids that have finished recently
            bad:        currentdirectory/bad  needs to be the mpids which failed recently

        Outputs:
            done:       set of stuff from done file
            bad:        set of stuff from bad  file
            modifies self.dictform base data to reflect the new changes
        """
        assert (attribute in ['absorb', 'static']), 'attribute needs to be one of 2 approved values'
        done = self.set_from_file('done')
        bad = self.set_from_file('bad')
        self.attribute_to_value(done, attribute, 1)
        self.attribute_to_value(bad,  attribute, 3)
        print('remember to set comments for bad')
        return(done, bad)
    def check(self, mpids, do_print=True, do_return=False):
        """Checks the currently stored values for an iterable of ids (or a single string)

        Inputs:
            mpids:      iterable of mpids to check
            do_print:   default True, prints output in easy to read format
            do_return:  default False, returns dictionary of entries

        Outputs:
            prints the status of those ids if do_print=True
            returns a smaller version of self.dictform with only the requested entries if do_return=True
        """
        if type(mpids) == str:
            mpids = [mpids]
        out = {i: self.dictform[i] for i in mpids}
        if do_print:
            outlines = ['{:<12}{:<7}{:<7}{:<12}{}'.format(self.dictform[i]['mpid'], self.dictform[i]['static'], 
                self.dictform[i]['absorb'], self.dictform[i]['origin'], self.dictform[i]['comments']) for i in mpids]
            header = '{:<12}{:<7}{:<7}{:<12}{}\n'.format('mpid', 'static', 'absorb', 'origin', 'comments')
            print(header)
            for item in outlines:
                print(item)
        if do_return:
            return(out)

    def add_origin(self, mpids, origin, resource_warnings = True, mpids_from_text = False, comments = ''):
        """Adds mpids of the new origin to the database

        Inputs:
            mpids:              set of mpids in the new list (filepath to text file if mpids_from_text = True)
            origin:             date of new set in 'Mmm-DD-YYYY' format
            resource_warnings:  default True, print warnings of bad runs that are no longer needed
            mpids_from_text:    default False, read mpids from a text file instead of input set

        Outputs:
            Modifies self.dictform (underlying data) to include new mpids with status 0 0 if they are new
            Anything already in dictionary is ignored
            if resource_warnings is on, it warns about mpids with status not 1 1 which are not in the new set
            new:                set of mpids that were added
            old:                set of mpids that already existed but were not in new one
            overlap:            set of mpids that already existed and were in new one
            give_up:            set of mpids that should be abandoned (returned regardless)
        """
        # Error handling and making sure mpids is a set of mp-###### strings
        if mpids_from_text:
            assert (type(mpids) == str), "using mpids_from_text requires filepath string for mpids input"
            mpids = self.set_from_file(mpids)
        if type(mpids) in {'list', 'tuple'}:
            mpids == set(mpids)
        assert (type(mpids) == set), "mpids should be input as a set"
        # Relevant subsets
        new = mpids - self.all_mpids
        old = self.all_mpids - mpids
        overlap = self.all_mpids.intersection(mpids)
        olddict = self.check(old, do_print = False, do_return = True)
        give_up = self.get_mpids('absorb', [1, -1], invert = True, source = olddict)
        # Modify base data with new data
        for mpid in new:
            self.dictform[mpid] = { 'mpid': mpid,
                                    'static': 0,
                                    'absorb': 0,
                                    'origin': origin,
                                    'comments': comments
                                    }
        # Print warnings
        if resource_warnings:
            print('Compounds to give up on:')
            self.check(give_up)
        # return results
        return(new, old, overlap, give_up)







###############################
#     Old high level funcs    #
###############################

def POTwrite(): 
    with open('POSCAR', 'r') as f:
        inlines = [i.rstrip('\n') for i in f.readlines()]
    matline = inlines[5].split()
    POSnames = [i + recommended_PAW[i] for i in matline]
    with open('POTCAR', 'w') as f:
        call(['cat'] + 
             ['/home/jdb/vasp_support/potentials/potpaw_PBE.54/{}/POTCAR'.format(i)
             for i in POSnames], stdout = f)
    return
    
    """ directory setup """
def directory_setup():
    for i in ['static', 'band', 'pbesoc', 'SOCband',
              'HSESOC', 'SHband']:
        if i not in os.listdir():
            os.mkdir(i)
        shutil.copy('../job', i)
    return
# def main(option):
#     if option == 'makePotcar':
#         print('Making POTCAR')
#         POTwrite()
#         return
#     elif option == 'plotBand':
#         stuf
#     elif option == 'plotOrbital':
#         stuf
# if __name__ == '__main__':
#     main(sys.argv[1])
def diffFile(file1, file2, file3):
    """Takes line contents of file1 - line contents of file2 (as sets) and writes it to file3
    May input a list of filenames to file2 to remove multiple files"""
    with open(file1, 'r') as f:
        start = set([i.rstrip('\n') for i in f.readlines()])
    if isinstance(file2, str):
        with open(file2, 'r') as f:
            remove = set([i.rstrip('\n') for i in f.readlines()])
        end = start - remove
    else:
        end = start
        for file in file2:
            with open(file, 'r') as f: #This could fail if elements of file2 are not strings
                remove = set([i.rstrip('\n') for i in f.readlines()])
            end = end - remove
    with open(file3, 'w') as f:
        f.writelines([i + '\n' for i in list(end)])


# Old Functions #
    # def KTEST(self, RUN = True):
    #     """ ktest portion """
    #     os.chdir('ktest')
    #     indict = start_incar(self.compound, 'ktest')
    #     indict['ENCUT'] = 1200
    #     write_INCAR(indict)
    #     for i in range(20,61,10):
    #         if '{:02d}'.format(i) not in os.listdir():
    #             os.mkdir('{:02d}'.format(i))
    #     for i in range(20,61,10):
    #         os.chdir('{:02d}'.format(i))
    #         write_KPOINTS(i)
    #         write_job(NODES = self.NODES, PPN = self.PPN, WALLTIME = self.WALLTIME)
    #         shutil.copy('../../POSCAR', '.')
    #         shutil.copy('../../POTCAR', '.')
    #         shutil.copy('../INCAR', '.')        
    #         os.chdir('../')
    #     if RUN:
    #         for i in range(20,61,10):
    #             os.chdir('{:02d}'.format(i))
    #             subjob()
    #             os.chdir('../')
    #         os.chdir ('../')
    # def ENTEST(self, RUN = True):
    #     """
    #     ENTEST PORTION

    #     Remember to check successful runs with wavecheck function (already in path)
    #     Simple solution:
    #     for i in {03..15..2}; do echo $i; tail -2 ${i}/out; done
        
    #     Better one:
    #     for i in {03..15..2}; do
    #     echo $i
    #     cat $i/OUTCAR | grep TOTEN | tail -1 | sed 's/free  energy   TOTEN  =      //g'
    #     echo -e '\n'
    #     done
        
    #     Remember to copy good kpoints file to base directory for each
    #     (Rb and K directories)
    #     """
    #     os.chdir('entest')
    #     shutil.copy('../ktest/KPOINTS', '../')
    #     for i in range(4, 12):
    #         if '{:02d}'.format(i) not in os.listdir():
    #             os.mkdir('{:02d}'.format(i))
    #     for i in range(4, 12):
    #         os.chdir('{:02d}'.format(i))
    #         write_job()
    #         shutil.copy('../../ktest/KPOINTS', '.')
    #         shutil.copy('../../ktest/KPOINTS', '../../')
    #         shutil.copy('../../POSCAR', '.')
    #         shutil.copy('../../POTCAR', '.')
    #         indict = start_incar(self.compound, 'entest')
    #         indict['ENCUT'] = 100*i
    #         write_INCAR(indict)
    #         os.chdir('../')
    #     if RUN:
    #         for i in range(4, 12):
    #             os.chdir('{:02d}'.format(i))
    #             subjob()
    #             os.chdir('../')
    #         os.chdir('../')


