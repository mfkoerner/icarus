#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 20:22:25 2017

@author: mitchell
"""


##########################################################
#  Module for crystal analysis
##########################################################

import icarus.config as config
import numpy as np, os
from pymatgen import MPRester
import pickle
import pymatgen.io.vasp as io
import pymatgen.electronic_structure.plotter as plotter
import collections
import math
import json
from zipfile import ZipFile


API_KEY = config.matprojapi



# Dictionary of effective ionic radii
#--------#--------#
'''
 Requires text file in format of:
     Element_1 ionic_radius_1
     Element_2 ionic_radius_2
     C 0.16
'''
#--------#--------#
def relpath(filelocation):
    return(os.path.join(os.path.dirname(__file__), filelocation))
def radii(filelocation):
    d = {}
    with open(filelocation) as f:
        for line in f:
            (key, val) = line.split()
            d[key] = float(val)
    return(d)

def reciprocal(Lattice):
    """Returns the reciprocal of a lattice
    Lattice: input as [[a1],[a2],[a3]]
    Output: [[b1],[b2],[b3]]"""
    tsprod = np.dot(Lattice[0], np.cross(Lattice[1], Lattice[2]))    #Triple scalar product
    rLattice = 2*np.pi * np.array([np.cross(Lattice[1], Lattice[2]), 
                               np.cross(Lattice[2], Lattice[0]), 
                               np.cross(Lattice[0], Lattice[1])]) / tsprod
    return(rLattice)
def rec_abc(a, b, c, Ain, Bin, Cin):
    """Returns the reciprocal of a set of lattice vectors and angles (a, b, c, alpha, beta, gamma)
    data: input as (a, b, c, alpha, beta, gamma)
    output: (a*, b*, c*, alpha*, beta*, gamma*, reciprocal volume, real volume)"""
    # a, b, c, A, B, C = data
    A = Ain * np.pi / 180 ; B = Bin * np.pi / 180 ; C = Cin * np.pi / 180
    V = np.sqrt(1 - np.cos(A)**2 - np.cos(B)**2 - np.cos(C)**2 + 2*np.cos(A)*np.cos(B)*np.cos(C)) * a * b * c
    d = 2 * np.pi * b * c * np.sin(A) / V
    e = 2 * np.pi * c * a * np.sin(B) / V
    f = 2 * np.pi * a * b * np.sin(C) / V
    D = np.arccos((np.cos(B)*np.cos(C) - np.cos(A)) / (np.sin(B) * np.sin(C)))
    E = np.arccos((np.cos(C)*np.cos(A) - np.cos(B)) / (np.sin(C) * np.sin(A)))
    F = np.arccos((np.cos(A)*np.cos(B) - np.cos(C)) / (np.sin(A) * np.sin(B)))
    # RV = np.sqrt(1 - np.cos(D)**2 - np.cos(E)**2 - np.cos(F)**2 + 2*np.cos(D)*np.cos(E)*np.cos(F)) * d * e * f   # Calculating reciprocal volume from vectors
    RV = 8 * np.pi**3 / V                   # Reciprocal volume as related to real volume
    D *= 180 / np.pi ; E *= 180 / np.pi ; F *= 180 / np.pi
    # print(d, e, f, D, E, F, RV, V)
    return(d, e, f, D, E, F, RV, V)

def separate3(data):
    """Separates data into 3 separate instances"""
    a = np.array([i[0] for i in data])
    b = np.array([i[1] for i in data])
    c = np.array([i[2] for i in data])
    return(a,b,c)

def vol_from_abc(a, b, c, Ain, Bin, Cin):
    '''returns volume from vectors and angles'''
    A = Ain * np.pi / 180 ; B = Bin * np.pi / 180 ; C = Cin * np.pi / 180
    V = np.sqrt(1 - np.cos(A)**2 - np.cos(B)**2 - np.cos(C)**2 + 2*np.cos(A)*np.cos(B)*np.cos(C)) * a * b * c
    return(V)

def download_cifs(mpids, zipfilename='cifs.zip', style='cifs.computed'):
    """deposits cif files from materials project into a zipped directory

    Inputs
        mpids: iterable of materials id strings such as 'mp-19' to get POSCARS for
        zipfilename: filename of zip folder created
        style: choose which cif to get from materialsproject

    Outputs
        places cif files into filepath/zipfilename
    """
    with MPRester(config.matprojapi) as mpr:
        docs = mpr.query({'material_id': {'$in': mpids}}, ['material_id', 'pretty_formula','cifs.computed'])
    with ZipFile(zipfilename, 'w') as f:
        for d in docs:
            f.writestr('{}.cif'.format(d['material_id']),d[style])


# def fullprint(*args, **kwargs):
#     from pprint import pprint
#     import numpy
#     opt = numpy.get_printoptions()
#     numpy.set_printoptions(threshold='nan')
#     print(*args, **kwargs)
#     numpy.set_printoptions(**opt)


# Flatten out a list
#--------#--------#
'''
Takes a series of nested lists and flattens it out
Eg. [1,[2,3,[4]],[5,6]] -> [1,2,3,4,5,6]
Credit: http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
'''
#--------#--------#
def flatten(l):
  out = []
  for item in l:
    if isinstance(item, (list, tuple)):
      out.extend(flatten(item))
    else:
      out.append(item)
  return out

def lengthsappend(lengths, label1, label2, specificDistance):
    try:
        lengths[label1].append([label2, specificDistance])
    except KeyError:
        lengths[label1] = [[label2, specificDistance]]
    try:
        lengths[label2].append([label1, specificDistance])
    except KeyError:
        lengths[label2] = [[label1, specificDistance]]

def is_odd(num):
    return num & 0x1


def swaprows(matrix, index1, index2):
    '''
    Only works with symmetric square matrices
    Swaps desired rows and columns
    '''
    out = np.zeros([len(matrix), len(matrix)])
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if j == i:
                if j == index1:
                    out[i,j] = matrix[index2,index2]
                if j == index2:
                    out[i,j] = matrix[index1,index1]
            if j > i:
                if i == index1:
                    k = index2
                elif i == index2:
                    k = index1
                else:
                    k = i
                if j == index1:
                    l = index2
                elif j==index2:
                    l = index1
                else:
                    l = j
                out[i, j] = matrix[k, l]
                out[j, i] = matrix[k, l]
    return out

def getHeader(filepath):
    '''Returns header of a file without the '# '
    or the newline character
    '''
    with open(filepath, 'r') as f:
        return(f.readline()[2:]).rstrip('\n')

def makeTextFile(compound, filename, lines, baseDir = '.'):
    '''Writes list of strings (lines) to a text file
    file path is baseDir/compound/filename
    '''
    with open('{}/{}/{}'.format(baseDir, compound, filename), 'w') as f:
        f.writelines([i+'\n' for i in lines])
def get_data(ids):
    with MPRester(API_KEY) as mpr:
        data = {i: mpr.get_data(i) for i in ids}
        data = {i: data[i][0] for i in data}
    return(data)
def nearly_equal(value1, value2, tolerance):
    return abs(value1-value2) <= tolerance
def header_from_structure(structure):
    elements = [str(i.specie) for i in structure]
    header = ', '.join(elements)
    return(header)
def pickleload(fp):
    with open(fp, 'rb') as f:
        return(pickle.load(f))
def get_BSobject(vasprun_fp = 'vasprun.xml'):
    run = io.outputs.BSVasprun(vasprun_fp)
    bandstructure = run.get_band_structure()
    return(bandstructure)
def quickplot(vasprun_fp = 'vasprun.xml'):
    bandstructure = get_BSobject(vasprun_fp)
    myplot = plotter.BSPlotter(bandstructure)
    myplot.show()
def importband(mpid):
    with MPRester(API_KEY) as mpr:
        return(mpr.get_bandstructure_by_material_id(mpid))
def get_vbm_loose(BS, cutoff = 0.026):
    """
    Returns data about the VBM.

    Returns:
        dict as {"band_index","kpoint_index","kpoint","energy"}
        - "band_index": A dict with spin keys pointing to a list of the
        indices of the band containing the VBM (please note that you
        can have several bands sharing the VBM) {Spin.up:[],
        Spin.down:[]}
        - "kpoint_index": The list of indices in BS.kpoints for the
        kpoint vbm. Please note that there can be several
        kpoint_indices relating to the same kpoint (e.g., Gamma can
        occur at different spots in the band structure line plot)
        - "kpoint": The kpoint (as a kpoint object)
        - "energy": The energy of the VBM
        - "projections": The projections along sites and orbitals of the
        VBM if any projection data is available (else it is an empty
        dictionnary). The format is similar to the projections field in
        BandStructure: {spin:{'Orbital': [proj]}} where the array
        [proj] is ordered according to the sites in structure
"""
    if BS.is_metal():
        return {"band_index": [], "kpoint_index": [],
                "kpoint": [], "energy": None, "projections": {}}
    max_tmp = -float("inf")
    index = None
    kpointvbm = None
    for spin, v in BS.bands.items():
        for i, j in zip(*np.where(v < BS.efermi)):
            if v[i, j] > max_tmp:
                max_tmp = float(v[i, j])
                index = j
                kpointvbm = BS.kpoints[j]

    list_ind_kpts = []
    if kpointvbm.label is not None:
        for i in range(len(BS.kpoints)):
            if BS.kpoints[i].label == kpointvbm.label:
                list_ind_kpts.append(i)
    else:
        list_ind_kpts.append(index)
    # get all other bands sharing the vbm
    list_ind_band = collections.defaultdict(list)
    for spin in BS.bands:
        for i in range(BS.nb_bands):
            if math.fabs(BS.bands[spin][i][index] - max_tmp) < 0.001:
                list_ind_band[spin].append(i)
            elif (max_tmp - BS.bands[spin][i][index] < cutoff) and (max_tmp - BS.bands[spin][i][index]) > 0 :
                list_ind_band[spin].append(i)
                print('!!!!!!! Extra VBM\n\n')
                print(BS.bands[spin][i][index] - max_tmp)
                print(i)
    proj = {}
    for spin, v in BS.projections.items():
        if len(list_ind_band[spin]) == 0:
            continue
        proj[spin] = v[list_ind_band[spin][0]][list_ind_kpts[0]]
    return {'band_index': list_ind_band,
            'kpoint_index': list_ind_kpts,
            'kpoint': kpointvbm, 'energy': max_tmp,
            'projections': proj}

def get_cbm_loose(BS, cutoff = 0.026):
    """
    Returns data about the CBM.

    Returns:
        {"band_index","kpoint_index","kpoint","energy"}
        - "band_index": A dict with spin keys pointing to a list of the
        indices of the band containing the VBM (please note that you
        can have several bands sharing the VBM) {Spin.up:[],
        Spin.down:[]}
        - "kpoint_index": The list of indices in BS.kpoints for the
        kpoint vbm. Please note that there can be several
        kpoint_indices relating to the same kpoint (e.g., Gamma can
        occur at different spots in the band structure line plot)
        - "kpoint": The kpoint (as a kpoint object)
        - "energy": The energy of the VBM
        - "projections": The projections along sites and orbitals of the
        VBM if any projection data is available (else it is an empty
        dictionnary). The format is similar to the projections field in
        BandStructure: {spin:{'Orbital': [proj]}} where the array
        [proj] is ordered according to the sites in structure
    """
    if BS.is_metal():
        return {"band_index": [], "kpoint_index": [],
                "kpoint": [], "energy": None, "projections": {}}
    max_tmp = float("inf")

    index = None
    kpointcbm = None
    for spin, v in BS.bands.items():
        for i, j in zip(*np.where(v > BS.efermi)):
            if v[i, j] < max_tmp:
                max_tmp = float(v[i, j])
                index = j
                kpointcbm = BS.kpoints[j]

    list_index_kpoints = []
    if kpointcbm.label is not None:
        for i in range(len(BS.kpoints)):
            if BS.kpoints[i].label == kpointcbm.label:
                list_index_kpoints.append(i)
    else:
        list_index_kpoints.append(index)

    # get all other bands sharing the cbm
    list_index_band = collections.defaultdict(list)
    for spin in BS.bands:
        for i in range(BS.nb_bands):
            if math.fabs(BS.bands[spin][i][index] - max_tmp) < 0.001:
                list_index_band[spin].append(i)
            elif (BS.bands[spin][i][index] - max_tmp < cutoff) and (BS.bands[spin][i][index] - max_tmp) > 0:
                list_index_band[spin].append(i)
                print('!!!!!!! Extra CBM\n\n')
                print(BS.bands[spin][i][index] - max_tmp)
                print(i)
    proj = {}
    for spin, v in BS.projections.items():
        if len(list_index_band[spin]) == 0:
            continue
        proj[spin] = v[list_index_band[spin][0]][list_index_kpoints[0]]

    return {'band_index': list_index_band,
            'kpoint_index': list_index_kpoints,
            'kpoint': kpointcbm, 'energy': max_tmp,
            'projections': proj}

###########################
# Section for rarity data #
###########################

# Grab data itself
with open(relpath('data/rarity/HHI_Production.json'), 'r') as f:
    HHI_P = json.load(f)
with open(relpath('data/rarity/HHI_Reserves.json'), 'r') as f:
    HHI_R = json.load(f)
with open(relpath('data/rarity/surface_abundance.json'), 'r') as f:
    abundance = json.load(f)
with open(relpath('data/rarity/weights.csv'), 'r') as f:
    indata = [line.rstrip('\n').split(',') for line in f.readlines()]
weights = {i[1]: float(i[3].strip(']').strip('[')) for i in indata[1:110]}
rarity_allowed_atoms = {'H','Li','Be','B','C','N','O','F','Na','Mg','Al','Si','P','S','Cl','K','Ca','Sc','Ti','V','Cr',
'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Rb','Sr','Y','Zr','Nb','Mo','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',
'I','Cs','Ba','La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi'}

# manipulate data to more convenient format
data = {}
for element in HHI_P:
    AN = element['x']
    form = element['formula']
    kzt = element['kZT']
    HHI_link = element['url']
    prod = element['y']
    reserve = 'none'
    for el2 in HHI_R:
        if el2['x'] == AN:
            reserve = el2['y']
    abu = 'none'
    abu_link = 'none'
    for el3 in abundance:
        if el3['x'] == AN:
            abu = el3['y']
            abu_link = el3['url']
    weight = weights[form]
    data[AN] = {
        'atomic_number': AN,
        'formula': form,
        'kZT': kzt,
        'HHI_Production': prod,
        'HHI_Reserve': reserve,
        'abundance': abu,
        'HHI_url': HHI_link,
        'abundance_url': abu_link,
        'weight': weight
    }
data_by_form = {i['formula']: i for i in data.values()}

# useful funcitons
def abundance_and_HHI(formula, strform = False):
    """Calculates crustal abundance and Herfindahi-Hirschman index for a given formula

    Inputs:     formula as a dictionary of {'element1': number1, 'element2': number2, ...}
                strform would allow one to input elements as they would in dict_formula_from_str defaults to false
                strform input looks like MgH2SO5 (exceptions in dict_formula_from_str documentation)

    Outputs:    tuple of (abundance, HHI_Production, HHI_Reserve)
                with abundance units of Parts Per Million (ppm)
                HHI can range from 0 (spread evenly over infinite countries) to 10,000 (all from one country)

    exceptions: Cannot handle Actinides
    """
    if strform:
        formula = dict_formula_from_str(formula)
    components = []
    for atom in formula.keys():
        # print(atom)
        data = data_by_form[atom]
        scarcity = 1/data['abundance']
        number = formula[atom]
        weight = data['weight']
        hhip = data['HHI_Production']
        hhir = data['HHI_Reserve']
        components.append([number * weight, scarcity, hhip, hhir])
    components = np.array(components).T       #makes it of format [[scarcity1, scarcity2], [totweight1, totweight2]]
    avgScarcity = np.sum(components[0] * components[1]) / np.sum(components[0])
    avgHHIP = np.sum(components[0] * components[2]) / np.sum(components[0])
    avgHHIR = np.sum(components[0] * components[3]) / np.sum(components[0])    
    return((1/avgScarcity, avgHHIP, avgHHIR))

def dict_formula_from_str(strform):
    """Creates a dictionary format of a compounds from a string format

    Inputs:     string of the standard form like "CsCl", "C6H12O", "CO2", "MgH2SO5"
                capitalization matters!
                cannot handle "()" so Na3Ca(BO2)5 would not work
                can handle double counting an element so CH3NH3PbI3 returns {'C': 1, 'I': 3, 'Pb': 1, 'H': 6, 'N': 1}

    Outputs:    formula as a dictionary of {'element1': number1, 'element2': number2, ...}

    Exceptions: Cannot handle pharentesis
    """
    UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    LOWER = 'abcdefghijklmnopqrstuvwxyz'
    NUMBER = '0123456789'
    upper = np.where([x in UPPER for x in strform])[0]
    lower = np.where([x in LOWER for x in strform])[0]
    number = np.where([x in NUMBER for x in strform])[0]
    out = {}
    for i in range(len(upper)):
        start = upper[i]
        # Finding element name and number
        try:
            end = upper[i+1]
        except IndexError:
            end = len(strform)
        try:
            numstart = min(set(range(start,end)).intersection(set(number)))
            nel = int(strform[numstart:end])
        except ValueError:
            nel = 1
            numstart = end
        element = strform[start:numstart]
        # end of element name and number
        if element in set(out.keys()):          #Checking for double counting
            out[element] += nel
        else:                                   #Non double counted (normal)
            out[element] = nel
    return(out)

