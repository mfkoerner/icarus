#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 18:07:28 2017

@author: mitchell

!!!!!!!
!!!!!!!
This code exists for legacy purposes. It is likely not functional. Please contact Mitchell Koerner (GitHub mfkoerner) for more details.

Requires special version of pymatgen/electronic_structure/bandstructure.py to work (has get_{v,c}bm_round)
I don't think this is relevant anymore, but it could be an issue if you are having issues
baseDir needs to be updated to a directory with bandsXXX.p pickled lists of bandstructure objects
!!!!!!!
!!!!!!!
"""

"""
Needed updates: (old)
    needs to be run
    # need formula (not pretty formula) as well
    26mev degeneracy allowance
"""
"""        outlist = [mpid,form,SGnum,PG,CS,OT,spol,EaH,Eg,Edg,dK,mh,me,V,*angles,
                   *vectors,colstring(Kcbmc),colstring(Kvbmc),colstring(Kcbmf),
                   colstring(Kvmbf), colstring(icsd_ids), colstring(rangles), colstring(rabc)]
"""
header = 'mpid,formula,sgno,pgsymbol,xtalsys,oxtype,spol,ehull,eg,edg,dk,mh,me,V,a,b,c,alpha,beta,gamma,formula_dict,k_cbm_abs,k_vbm_abs,k_cbm_frac,k_vmb_frac,icsd,rangles,rabc'

blacklist = ['mp-562720']

import config
import pickle
import time
import numpy as np
import os
from pymatgen import MPRester
from bisect import bisect
import scipy.optimize
import mymatgen as mmg
import importlib as imp
imp.reload(mmg)

minpoints = 2
maxpoints = 2
melect = 9.109e-31
hbar = 1.054e-34
aunits = 1.602e-19 * 1e-20

#import scipy.optimize
#from bisect import bisect

baseDir = 'INPUT_DIRECTORY_HERE'
os.chdir(baseDir)

mpr = MPRester(config.matprojapi)
results = mpr.query(criteria={"has_bandstructure" : True},properties=[
        "material_id", "spacegroup", "e_above_hull", "oxide_type", "pretty_formula",
        "volume", "cif", "icsd_ids", "formula"])
ids = [i['material_id'] for i in results]
idreferenced = {i['material_id']: i for i in results}

def colstring(inlist):
    return(':'.join([str(x) for x in inlist]))
def Func(x, a):
    return(a*(x)**2)

def Curve(xs, ys, xf, yf):
    x = [np.linalg.norm(i - xf) for i in xs]
    y = [i - yf for i in ys]
    return(scipy.optimize.curve_fit(Func, x, y)[0][0])

def effmass(BS, extrema):
    if len(extrema['band_index']) == 0:
        raise AttributeError('Band Structure has no band gap')
    out = []
    sindex = [i['start_index'] for i in BS.branches]
    print('sindex', sindex)
    eindex = [i['end_index'] for i in BS.branches]
    print('eindex', eindex)
    for key in extrema['band_index'].keys():
        bandids = extrema['band_index'][key]
        for bandid in bandids:
            band = BS.bands[key][bandid]
            print('band', band)
            for kindex in extrema['kpoint_index']:
                print('working on kindex', kindex)
                if kindex in sindex:
                    branchid = sindex.index(kindex)
                if kindex in eindex:
                    branchid = eindex.index(kindex)
                else:
                    branchid = bisect(eindex, kindex)
                branch = list(range(sindex[branchid], eindex[branchid] + 1))
                print('branch', branch)
                nleft = branch.index(kindex)
                nright = len(branch) - nleft - 1
                if nleft > minpoints:
                    if nleft < maxpoints:
                        myslice = slice(0, nleft + 1)
                    else:
                        myslice = slice(nleft - maxpoints, nleft + 1)
                    curve = Curve([BS.kpoints[i].cart_coords for i in branch[myslice]],
                                  band[branch[myslice]], extrema['kpoint'].cart_coords,
                                  band[kindex])
                    print([BS.kpoints[i].cart_coords for i in branch[myslice]],
                                  band[branch[myslice]], extrema['kpoint'].cart_coords,
                                  band[kindex])
                    out.append(curve)
                if nright > minpoints:
                    if nright < maxpoints:
                        myslice = slice(nleft, len(branch) + 1)
                    else:
                        myslice = slice(nleft, nleft + maxpoints + 1)
                    curve = Curve([BS.kpoints[i].cart_coords for i in branch[myslice]],
                                  band[branch[myslice]], extrema['kpoint'].cart_coords,
                                  band[kindex])
                    print([BS.kpoints[i].cart_coords for i in branch[myslice]],
                                  band[branch[myslice]], extrema['kpoint'].cart_coords,
                                  band[kindex])
                    out.append(curve)
    out = np.array(out)
    return(out)



def main(n = 1):
    os.chdir(baseDir)
    t1 = time.perf_counter()
    with open('bands{:03d}.p'.format(n), 'rb') as f:
        bands = pickle.load(f)
    print('bands: {}s'.format(time.perf_counter()-t1))
    partialids = list(bands.keys())
    outlines = []
    for mpid in partialids:
        if mpid not in blacklist:
            print('working on', mpid)
            myband = bands[mpid]
            try:
                Kcbmc = myband.get_cbm()['kpoint'].cart_coords
                Kvbmc = myband.get_vbm()['kpoint'].cart_coords
                Kcbmf = myband.get_cbm()['kpoint'].frac_coords
                Kvmbf = myband.get_vbm()['kpoint'].frac_coords            
                dK = np.linalg.norm(Kcbmc-Kvbmc)
                myvbm = mmg.get_vbm_loose(myband)
                mycbm = mmg.get_cbm_loose(myband)
                vcurve = effmass(myband, myvbm) * -1
                ccurve = effmass(myband, mycbm)
                if len(vcurve) != 0:
                    mhlist = [abs(hbar**2 / (2 * melect * i * aunits)) for i in vcurve]
                    mh = ':'.join([str(i) for i in mhlist])
                else:
                    mh = np.nan
                if len(ccurve) != 0:
                    melist = [abs(hbar**2 / (2 * melect * i * aunits)) for i in ccurve]
                    me = ':'.join([str(i) for i in melist])
                else:
                    me = np.nan
            except AttributeError:
                Kcbmc, Kvbmc, Kcbmf, Kvmbf, dK, mh, me = ([np.nan]*3, [np.nan]*3,
                            [np.nan]*3, [np.nan]*3, np.nan, np.nan, np.nan)
            Eg = myband.get_band_gap()['energy']
            Edg = myband.get_direct_band_gap()
            try:
                data = idreferenced[mpid]
            except KeyError:
                data = mpr.get_data(mpid)[0]
            form = data['pretty_formula']
            SGnum = data['spacegroup']['number']
            PG = data['spacegroup']['point_group']
            CS = data['spacegroup']['crystal_system']
            OT = data['oxide_type']
            EaH = data['e_above_hull']
            V = data['volume']
            cif = [i.split() for i in data['cif'].split('\n')]
            lengths = np.array(cif[3:9])
            angles = lengths[:3,1]
            vectors = lengths[3:,1]
            spol = myband.is_spin_polarized
            icsd_ids = data['icsd_ids']
            rangles = myband.lattice_rec.angles
            rabc = myband.lattice_rec.abc
            try:
                formula = data['formula']
            except:
                formula = np.nan
            outlist = [mpid,form,SGnum,PG,CS,OT,spol,EaH,Eg,Edg,dK,mh,me,V,*vectors,
                       *angles,formula,colstring(Kcbmc),colstring(Kvbmc),colstring(Kcbmf),
                       colstring(Kvmbf), colstring(icsd_ids), colstring(rangles), colstring(rabc)]
            strlist = [str(i) for i in outlist]
            outstring = ','.join(strlist)
            outlines.append(outstring)
    return(outlines)
'''
if __name__ == '__main__':
    os.chdir(baseDir)
    for n in range(213, 266):       #266 for real run and change baseDir back
        print('working on', n)
        with open('spec_data{:03d}.txt'.format(n), 'w') as f:
            f.writelines([i + '\n' for i in main(n)])
    with open('header.txt', 'w') as f:
        f.write(header)'''
#if __name__ == '__main__':
#    myband = bands['mp-532575']
#    emass = effmass(myband, myband.get_vbm())
            
'''
for n in range(9, 76):
    with open('spec_data{:03d}.txt'.format(n), 'w') as f:
        f.writelines([i + '\n' for i in main(n)])
'''
'''
mpid,formula,SGn,PG,CS,oxide,EaH,E_gap,E_digap,dK,cv*,cc*,VolumeofUnitCell,alpha,beta,gamma,a,b,c,K_cbm(cart),K_vbm(cart),K_cbm(frac),K_vmb(frac) frac and abs k_cbm
'''
