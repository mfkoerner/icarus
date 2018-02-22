#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 14:18:39 2017

@author: mitchell
"""

import config
import PyQt5
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.patches as patches

plt.style.use(config.plot_style)

#-------------------------------------------------------#
# SET OUTLIERS, REFERENCE COMPOUNDS, AESTHETIC SETTINGS #
#-------------------------------------------------------#

outliers = ['mp-14134', 'mp-769844', 'mp-11585', 'mp-758317', 'mp-20526']
ref = {
    5238:'CuGaS$_2$',
    406:'CdTe',
    149:'Si',
    2534:'GaAs',
    804:'GaN',
    2133:'ZnO',
    2624:'AlSb'
}
scaleDot = 1e2
c1, c2, c3, c4, c5 = ('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd')
fs = 24 # annotation fontsize
ms = 50

def gtrans(txt):
    if txt != '\\Gamma':
        return(txt)
    else:
        return(r'$\Gamma$')

def sizeplot(xdata, ydata, sizedata, scale = 1):
    f, ax = plt.subplots()
    return(f, ax, ax.scatter(xdata, ydata, s = np.array(sizedata) * scale, 
               facecolors = 'none', edgecolors = 'b'))

def nan2str(arg):
    if arg is np.nan:
        return('nan')
    else:
        return(arg)

def nan2nan(arg):
    if arg == 'nan' or arg == 'none' or arg == 'none5':
        return(np.nan)
    else:
        return(arg)

def scatter(xdata, ydata, c = 'k', ax = None, linewidth = None, alpha = None,
            marker = 'o'):
    new = ax is None
    if new:
        f, ax = plt.subplots()
    if new:
        return f, ax, ax.scatter(xdata, ydata, marker = marker, facecolors = 'none', edgecolors = c,
               linewidth = linewidth, alpha = alpha)
    else:
        return ax.scatter(xdata, ydata, marker = marker, facecolors = 'none', edgecolors = c,
               linewidth = linewidth, alpha = alpha)
def zoomscatter(xdata, ydata, zoom = 2, limits = None, zlimits = None, marker = 'o', loc = 1, size = None, scale = 1):
    if scale != 1:
        size = np.array([i * scale for i in size])
    if limits is None:
        limits = [min(xdata), max(xdata), min(ydata), max(ydata)]
    f, ax = plt.subplots()
    ax.scatter(xdata, ydata, s = size)
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    axins = zoomed_inset_axes(ax, zoom, loc = loc)
    if zlimits is None:
        xdis = limits[1] - limits[0]
        ydis = limits[3] - limits[2]
        zlimits = [limits[0] + .4*xdis, limits[0] + .6*xdis, 
        limits[2] + .4*ydis, limits[2] + .6*ydis]
    axins.set_xlim(zlimits[0], zlimits[1])
    axins.set_ylim(zlimits[2], zlimits[3])
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    axins.tick_params(axis = 'both', which = 'both', left = 'off', bottom = 'off')
    axins.scatter(xdata, ydata, marker = marker, s = size)
    mark_inset(ax, axins, loc1=2, loc2=4, fc = 'none', ec='0.5')
    return(f, ax, axins)
def zoomscatterbicolor(xdata, ydata, select, zoom = 2, limits = None, zlimits = None, marker = 'o', loc = 1, size = None, scale = 1, c = 'b', cs = 'r'):   #WIP
    if scale != 1:
        size = np.array([i * scale for i in size])
    if limits is None:
        limits = [min(xdata), max(xdata), min(ydata), max(ydata)]
    f, ax = plt.subplots()
    sel = ax.scatter(xdata[select], ydata[select], s = size, c = cs)
    unsel = ax.scatter(xdata[~select], ydata[~select], s = size, c = c)
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    axins = zoomed_inset_axes(ax, zoom, loc = loc)
    if zlimits is None:
        xdis = limits[1] - limits[0]
        ydis = limits[3] - limits[2]
        zlimits = [limits[0] + .4*xdis, limits[0] + .6*xdis, 
        limits[2] + .4*ydis, limits[2] + .6*ydis]
    axins.set_xlim(zlimits[0], zlimits[1])
    axins.set_ylim(zlimits[2], zlimits[3])
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    axins.tick_params(axis = 'both', which = 'both', left = 'off', bottom = 'off')
    axins.scatter(xdata, ydata, marker = marker, s = size)
    mark_inset(ax, axins, loc1=2, loc2=4, fc = 'none', ec='0.5')
    return(f, ax, axins)
def bandplot(kdist, energy, mpid, formula, dirChange, pathBreak, xrange = None, yrange = None):
    if len(energy) == 1:
        energyfull = energy[1]
    elif len(energy) == 2:
        print(energy)
        energyfull = np.array([energy[1], energy[-1]])
    if xrange is None:
        xrange = (np.min(kdist), np.max(kdist))
    if yrange is None:
        yrange = (np.min(energyfull), np.max(energyfull))
    ymin, ymax = yrange
    xmin, xmax = xrange
    f, ax = plt.subplots()
    ax.plot(kdist, energy[1], '-' , c=( 0.121569 , 0.466667 , 0.705882 ) , lw=3)
    if len(energy) == 2:
        ax.plot(kdist, energy[-1], 'r--', lw=3)
    for k in dirChange:
        ax.plot( (k,k) , (ymin,ymax) , 'k-' , lw=1 )
    for k in pathBreak:
        ax.plot( (k,k) , (ymin,ymax) , 'k-' , lw=2 )
    ax.plot( (xmin,xmax) , (0,0) , 'k--' , lw=1 )     #fermi level
    ax.set_xlim( xmin , xmax )
    if yrange is None:
        ax.set_ylim( ymin , ymax )
    else:
        ax.set_ylim(yrange)
    ax.set_ylabel(r'$E-E_f$ (eV)')
    ax.tick_params( axis='x' , which='both' , bottom='off' , top='off' , labelbottom='off' )
    ax.set_title('{}   {}'.format(mpid, formula))
    return(f, ax)
def get_kdist(BS):
    k = BS.kpoints
    cart = [i.cart_coords for i in k]
    kdist = [0.]
    dirChange, pathBreak = [], []
    for i in range(1, len(cart)):
        kstep = np.linalg.norm(cart[i] - cart[i-1])
        if k[i].label is not None and k[i-1].label is not None and i != 1:
            if k[i].label != k[i-1].label:
                pathBreak.append( kdist[-1] )
                kdist.append( kdist[-1] )
            else:
                dirChange.append( kdist[-1] )
                kdist.append( kdist[-1] )
        else:
            kdist.append( kdist[-1] + kstep )
    pathBreak = np.array(pathBreak)
    pathBreak = pathBreak[np.where(pathBreak>0)]
    pathBreak = np.array(list(set(pathBreak)))
    kdist = np.array(kdist)
    dirChange = np.array(list(set(dirChange)))
    if not BS.is_spin_polarized:
        assert(len(BS.bands)) == 1
        key = list(BS.bands.keys())[0]
        energy = BS.bands[key]
        energy += -BS.efermi
        energy = {1: energy.T}
    else:
        assert(len(BS.bands)) == 2
        keys = list(BS.bands.keys())
        energy = {1: BS.bands[keys[0]], -1 : BS.bands[keys[1]]}
        energy = {i: energy[i].T -BS.efermi for i in energy}
    return(kdist, energy, dirChange, pathBreak)
def labelBZ( ax , kdist , dirChange , pathBreak , labels, fdy = 0.03):
    yl = ax.get_ylim()
    joined = np.sort( np.concatenate((
        np.array([np.min(kdist)]),
        dirChange,pathBreak,
        np.array([np.max(kdist)]))) )
    for j , v in enumerate(joined):
        if v in pathBreak:
            s = '|'.join((labels.pop(0),labels.pop(0)))
        else:
            s = labels.pop(0)
        ax.text( v , yl[0]-(yl[1]-yl[0])*fdy , s , horizontalalignment='center' )
def getlabels(BS):
    out = []
    kold = ''
    for k in BS.kpoints:
        if k.label is not None:
            if k.label is not kold:
                out.append(k.label)
                kold = k.label
    out = [gtrans(i) for i in out]
    return(out)
def print_full_wide(x):
    pd.set_option('display.max_columns', x.shape[1])
    print(x)
    pd.reset_option('display.max_columns')
def colon_sep_nums(x):
    return([float(i) for i in str(x).split(':')])
def round2(x):
    return(np.round(x, 2))
def round3(x):
    return(np.round(x, 3))
def round4(x):
    return(np.round(x, 4))
def plot_eg_dk_m(df):
    '''plot dk vs. E_g with 1/mean(m_e) and 1/mean(m_h) dot size'''
    fig, ax = plt.subplots()
    ax.scatter(df['eg'], df['dk'], s=1./df['memean']*scaleDot, c=c1, label='', lw=1, alpha=0.2)
    ax.scatter(df['eg'] ,df['dk'] ,s=1./df['mhmean']*scaleDot, c=c2, label='', lw=1, alpha=0.1)
    ax.scatter(-10, -10, c=c1, s=100, label='electrons', lw=1, alpha=1)
    ax.scatter(-10, -10, c=c2, s=100, label='holes', lw=1, alpha=1)
    ax.set_xlabel('$E_g$ (eV)')
    ax.set_ylabel('$\Delta k$ ($\AA^{-1}$)')
    ax.set_ylim(-0.2,2)
    ax.set_xlim(-0.5,5)
    ax.set_title('marker size inversely proportional to $m^*$')
    ax.legend()
    return fig , ax

def plot_mh_me(DF):
    '''plot mean(m_h) vs. mean(m_e) with E_g dot size'''
    df = DF.ix[DF.eg>0] # remove metals
    fig, ax = plt.subplots()
    ax.loglog((1e-2,1e5), (1e-2,1e5), 'k--')
    sh = ax.scatter(df['memean'], df['mhmean'], s=df['eg']*scaleDot, c=c1, label='', lw=1, alpha=0.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('avg($m_e^*$) $m_e^{-1}$')
    ax.set_ylabel('avg($m_h^*$) $m_e^{-1}$')
    ax.set_xlim(1e-2,1e5)
    ax.set_ylim(1e-2,1e5)
    ax.text(4e2, 3e-2, '{:.0f}% have $m_e^*$ < $m_h^*$'.format(len(df.ix[df.memean<df.mhmean])/len(df)*100), size=fs)
    ax.set_title('marker size proportional to bandgap')
    return fig , ax

def plot_eg_m(DF):
    '''plot mean(m_e) and mean(m_h) vs. E_g'''
    df = DF.ix[DF.eg>0] # remove metals
    fig , ax = plt.subplots()
    ax.scatter(df['eg'], df['memean'], c=c1, s=ms, label='', lw=1, alpha=.4)
    ax.scatter(df['eg'], df['mhmean'], c=c2, s=ms, label='', lw=1, alpha=.4)
    ax.scatter(-1, -1, c=c1, s=100, label='$m_e^*$', lw=1, alpha=1)
    ax.scatter(-1, -1, c=c2, s=100, label='$m_h^*$', lw=1, alpha=1)
    ax.set_yscale('log')
    ax.set_xlabel('$E_g$ (eV)')
    ax.set_ylabel('avg($m^*$) $m_e^{-1}$')
    ax.set_xlim(0, 10)
    ax.set_ylim(1e-1,1e2)
    p = ax.add_patch(patches.Rectangle((2, 1.5e-1), 9.5-2, 1-1.5e-1, fc='none', ec='k', ls='--'))
    ax.text(4, 1.8e-1, 'transparent conductor candidates', size=fs)
    p.set_zorder(10)
    ax.legend(loc=1)
    return fig , ax
def old2newmpid(x):
    return('mp-{}'.format(x))
def dict_from_semi(x):
    try:
        if np.isnan(x):
            return(np.nan)
    except:
        pass
    return(eval(x.replace(';',',')))
def read_POSCAR(inlines):
    return(np.array([[float(i) for i in j.strip(' ').split()] for j in inlines.rstrip('\n').split('\n')]))
def vol(x):
    return(np.dot(np.cross(x[0], x[1]), x[2]))