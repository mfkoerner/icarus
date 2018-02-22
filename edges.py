########################################
# written for Python 3                 #
# by Doug Fabini (fabini@mrl.ucsb.edu) #
########################################

'''

Module description...

kd, E, P = el.spaghetti(ax, bd, ct, El=(-2, 2), vbmRef=True)

******** still need to handle path break case, extremely short k-leg case

'''

import electrons as el
import numpy as np
import matplotlib.pyplot as plt
try:
	# filepath to matplotlib stylesheet on work computer
	plt.style.use(['E:/Dropbox/Materials/Code/sty-plot'])
except:
	# filepath to matplotlib stylesheet on home computer
	plt.style.use(['C:/Users/dhfab/Dropbox/Materials/Code/sty-plot'])
from scipy.optimize import curve_fit



# ******************* #
# LOW LEVEL FUNCTIONS #
# ******************* #

def quad(x, A):
	''' One parameter quadratic function about x = 0 '''
	return A*x**2

def curv2effMass(A):
	''' Calculate effective mass (in units of electron rest mass) from parabolic band curvature (in units of eV Ang^2) '''
	hbar = 1.055e-34	# Planck constant (J s)
	m0 = 9.1094e-31  # electron mass (kg)
	Ahat = A*1.602e-19/1e10**2 # parabolic curvature converted to units of J m^2
	return hbar**2/2/Ahat/m0

def fitParabola(x, y, ii, ax=None, n=3, left=False):
	''' Estimate carrier effective mass along band
		x :    Momentum (k) in inverse Angstroms
		y :    Energy (E) in eV
		ii :   Index of band extremum, at momentum x[ii]
		ax :   (optional) Axes handle to plot parabolic fit into
		n :    (default) Number of points to use for fitting, including extremum
		left : (default) Boolean, fit to left along E(k) dispersion diagram?
	'''
	x0, y0 = x[ii], y[ii]
	if left:
		xx, yy = x[ii-n : ii+1], y[ii-n : ii+1]
		xxx = x[ii-10: ii+1]
	else:
		xx, yy = x[ii : ii+n+1], y[ii : ii+n+1]
		xxx = x[ii: ii+11]
	popt, pcov = curve_fit(quad, xx-x0, yy-y0)
	if ax:
		ax.plot(xxx, y0+quad(xxx-x0, popt[0]), '-', c='#A83425')
	return curv2effMass(popt[0])

def fitDirection(Eband, klist, kd, pb, K):
	''' Determine which direction (or both) to fit parabolae to for each extremum k point
		Eband : vector of E(k) for a single band
		klist : list of extremum k point indices
	'''
	left, right = [], []
	Nk = Eband.size
	for i, k in enumerate(klist):
		if k == 0:
			right.append(k)			# case 1, first point along path
		elif k == Nk-1:
			left.append(k)			# case 2, last point along path
		elif kd[k] not in pb:
			if k+1 in klist:
				left.append(k)		# case 3a, special point, direction change (left)
			elif k-1 in klist:
				right.append(k)		# case 3b, special point, direction change (right)
			else:
				left.append(k)		# case 6, general point
				right.append(k)
		elif kd[k] in pb:
			dleft = np.linalg.norm(K[k]-K[k-1])
			dright = np.linalg.norm(K[k]-K[k+1])
			if dright > dleft:
				left.append(k)		# case 4, special point, path break (left)
			elif dleft > dright:
				right.append(k)		# case 5, special point, path break (right)
			else:
				raise Exception('LEFT AND RIGHT DISTANCES EQUAL AT PATH BREAK')
		else:
			raise Exception('COULD NOT DETERMINE K-POINT TYPE FOR FITTING DIRECTION ROUTINE')
	return left, right

def trueExtremum(jj, Eband, maximum=True):
	''' Return only indices of true extrema (eliminate false positives from energy tolerance) 
		jj :      List of momentum indices to check
		Eband :   vector of E(k) for the single band of interest
		maximum : (default) Booelan, check for maximum (valence bands) or minimum (conductions bands)?
	'''
	jjj = []
	other = {True: -np.inf, False: np.inf}[maximum]
	# if maximum:
	# 	other = -np.inf
	# else:
	# 	other = np.inf
	for jji in jj:
		if jji-1 < 0:
			prev = other
			subseq = Eband[jji+1]
		elif jji+1 > Eband.size-1:
			subseq = other
			prev = Eband[jji-1]
		else:
			subseq = Eband[jji+1]
			prev = Eband[jji-1]
		if (Eband[jji] <= prev) & (Eband[jji] <= subseq) & (not maximum):
			jjj.append(jji)
		elif(Eband[jji] >= prev) & (Eband[jji] >= subseq) & maximum:
			jjj.append(jji)
	return jjj



# ******************** #
# HIGH LEVEL FUNCTIONS #
# ******************** #

def masses(kd, E, K, pb, ax=None, Etol=1e-2):
	''' Estimate carrier effective masses using the parabolic approximation, optionally plot fits
		kd :   Momentum coordinate along BZ high symmetry path (output of electrons.momentumCoord)
		E :    Array of energy eigenvalues (optionally relative to VBM)
		K :    Dictionary of k-point indices along path and reciprocal space coordinates
		pb :   Array of k-distances of path break points
		ax :   (optional) Axes to plot parabolic fits into
		Etol : (default) Energy tolerance for identifying band extrema
	'''
	Evbm = np.max(E[E<=0])
	Ecbm = np.min(E[E>0])
	vb = np.where(np.any(np.isclose(E, Evbm, atol=Etol), axis=0))[0]
	cb = np.where(np.any(np.isclose(E, Ecbm, atol=Etol), axis=0))[0]
	bands = np.concatenate((cb, vb))
	bandType = ['cb',]*cb.size
	bandType.extend(['vb',]*vb.size)
	Eextremum = {'vb': Evbm, 'cb': Ecbm}
	maxFlag = {'vb': True, 'cb': False}
	bandLabels = {'vb': 'VALENCE', 'cb': 'CONDUCTION'}
	me, mh = [], []
	print('\n *** PARABOLIC BAND FITTING REPORT ***')
	for b, bt in zip(bands, bandType):
		print('   \n{} BAND INDEX {}'.format(bandLabels[bt], b))
		Eband = E[:, b]
		jj = np.where(np.isclose(Eband, Eextremum[bt], atol=Etol))[0]
		jj = trueExtremum(jj, Eband, maximum=maxFlag[bt])
		left, right = fitDirection(Eband, jj, kd, pb, K)
		for ii in left:
			mii = fitParabola(kd, Eband, ii, ax=ax, left=True)
			if bt is 'cb':
				me.append(mii)
			else:
				mh.append(-mii)
			print( '     k-point #{:>4}  (LEFT) : {:7.2f}'.format(ii, mii) )
		for ii in right:
			mii = fitParabola(kd, Eband, ii, ax=ax, left=False)
			if bt is 'cb':
				me.append(mii)
			else:
				mh.append(-mii)
			print( '     k-point #{:>4} (RIGHT) : {:7.2f}'.format(ii, mii) )
	return me, mh

# fin.
