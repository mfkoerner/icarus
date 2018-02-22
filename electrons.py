########################################
# written for Python 3                 #
# by Doug Fabini (fabini@mrl.ucsb.edu) #
########################################

'''

High-level functions are at bottom of file:

	spaghetti : extracts and plots band structure (optionally, with orbital character) from VASP output files
				requires POSCAR, DOSCAR, OUTCAR, PROCAR files

	spiky     : extracts and plots density of states (optionally, with orbital character) from VASP output files
				requires POSCAR, DOSCAR, OUTCAR files

This module *should* intelligently detect and handle:
 - whether or not spin-orbit coupling is on,
 - whether the orbital projections are n,l,m quantum number-resolved, or only n,l quantum number-resolved, and
 - whether the POSCAR uses the "new" format (e.g. vasp.5.4.4) or the "old" format (e.g. vasp.5.3.3).

This module does *not*:
 - Handle spin-polarized calculations, or
 - Pull spin-texture out of calculations with spin-orbit coupling,
though both of these functionalities would be fairly easily added.

See description of function arguments for 'spaghetti' and 'spiky' in the respective functions at the end of this file.

'''



<<<<<<< HEAD
import config
import os
import numpy as np
import matplotlib.pyplot as plt
plt.style.use(config.plot_style)
=======
import os
import numpy as np
import matplotlib.pyplot as plt
try:
	# filepath to matplotlib stylesheet on work computer
	plt.style.use(['E:/Dropbox/Materials/Code/sty-plot'])
except:
	# filepath to matplotlib stylesheet on home computer
	plt.style.use(['C:/Users/dhfab/Dropbox/Materials/Code/sty-plot'])
>>>>>>> newdoug



# ******* #
# LOOKUPS #
# ******* #

# list of colors for plotting orbital projections, can easily be customized
cList = ('#FF7E0E', '#1F77B3', '#2BA02B')

# dictionary lookup for labels of Brillouin Zone special points, adapted from DOI:10.1016/j.commatsci.2010.05.010
bzLabels = {
	'CUB': ['\Gamma', 'X', 'M', '\Gamma', 'R', 'X', 'M', 'R'],
	'FCC': ['\Gamma', 'X', 'W', 'K', '\Gamma', 'L', 'U', 'W', 'L', 'K', 'U', 'X'],
	'BCC': ['\Gamma', 'H', 'N', '\Gamma', 'P', 'H', 'P', 'N'],
	'TET': ['\Gamma', 'X', 'M', '\Gamma', 'Z', 'R', 'A', 'Z', 'X', 'R', 'M', 'A'],
	'BCT1': ['\Gamma', 'X', 'M', '\Gamma', 'Z', 'P', 'N', 'Z_1', 'M', 'X', 'P'],
	'BCT2': ['\Gamma', 'X', 'Y', '\Sigma', '\Gamma', 'Z', '\Sigma_1', 'N', 'P', 'Y_1', 'Z', 'X', 'P'],
	'ORC': ['\Gamma', 'X', 'S', 'Y', '\Gamma', 'Z', 'U', 'R', 'T', 'Z', 'Y', 'T', 'U', 'X', 'S', 'R'],
	'ORCF1': ['\Gamma', 'Y', 'T', 'Z', '\Gamma', 'X', 'A_1', 'Y', 'T', 'X_1', 'X', 'A', 'Z', 'L', '\Gamma'],
	'ORCF2': ['\Gamma', 'Y', 'C', 'D', 'X', '\Gamma', 'Z', 'D_1', 'H', 'C', 'C_1', 'Z', 'X', 'H_1', 'H', 'Y', 'L', '\Gamma'],
	'ORCF3': ['\Gamma', 'Y', 'T', 'Z', '\Gamma', 'X', 'A_1', 'Y', 'X', 'A', 'Z', 'L', '\Gamma'],
	'ORCI': ['\Gamma', 'X', 'L', 'T', 'W', 'R', 'X_1', 'Z', '\Gamma', 'Y', 'S', 'W', 'L_1', 'Y', 'Y_1', 'Z'],
	'ORCC': ['\Gamma', 'X', 'S', 'R', 'A', 'Z', '\Gamma', 'Y', 'X_1', 'A_1', 'T', 'Y', 'Z', 'T'],
	'HEX': ['\Gamma', 'M', 'K', '\Gamma', 'A', 'L', 'H', 'A', 'L', 'M', 'K', 'H'],
	'RHL1': ['\Gamma', 'L', 'B_1', 'B', 'Z', '\Gamma', 'X', 'Q', 'F', 'P_1', 'Z', 'L', 'P'],
	'RHL2': ['\Gamma', 'P', 'Z', 'Q', '\Gamma', 'F', 'P_1', 'Q_1', 'L', 'Z'],
	'MCL': ['\Gamma', 'Y', 'H', 'C', 'E', 'M_1', 'A', 'X', '\Gamma', 'Z', 'D', 'M', 'Z', 'A', 'D', 'Y', 'X', 'H_1'],
	'MCLC1': ['\Gamma', 'Y', 'F', 'L', 'I', 'I_1', 'Z', 'F_1', 'Y', 'X_1', 'X', '\Gamma', 'N', 'M', '\Gamma'],
	'MCLC2': ['\Gamma', 'Y', 'F', 'L', 'I', 'I_1', 'Z', 'F_1', 'N', '\Gamma', 'M'],
	'MCLC3': ['\Gamma', 'Y', 'F', 'H', 'Z', 'I', 'F_1', 'H_1', 'Y_1', 'X', '\Gamma', 'N', 'M', '\Gamma'],
	'MCLC4': ['\Gamma', 'Y', 'F', 'H', 'Z', 'I', 'H_1', 'Y_1', 'X', '\Gamma', 'N', 'M', '\Gamma'],
	'MCLC5': ['\Gamma', 'Y', 'F', 'L', 'I', 'I_1', 'Z', 'H', 'F_1', 'H_1', 'Y_1', 'X', '\Gamma', 'N', 'M', '\Gamma'],
	'TRI1a': ['X', '\Gamma', 'Y', 'L', '\Gamma', 'Z', 'N', '\Gamma', 'M', 'R', '\Gamma'],
	'TRI2a': ['X', '\Gamma', 'Y', 'L', '\Gamma', 'Z', 'N', '\Gamma', 'M', 'R', '\Gamma'],
	'TRI1b': ['X', '\Gamma', 'Y', 'L', '\Gamma', 'Z', 'N', '\Gamma', 'M', 'R', '\Gamma'],
	'TRI2b': ['X', '\Gamma', 'Y', 'L', '\Gamma', 'Z', 'N', '\Gamma', 'M', 'R', '\Gamma'],

	'CUSTOM1': ['\Gamma', 'X', 'U', 'L', '\Gamma'],
	'CUSTOM2': ['L', '\Gamma', 'Z', 'P$_1$', 'F', 'Q'],
	'CUSTOM3': ['\Gamma', 'Y', 'D', 'Z', '\Gamma'],
	'CUSTOM4': ['\Gamma', 'M', 'L', 'A', '\Gamma']
}

def labelify(raw):
	''' Format BZ special point labels '''
	conv = {
		'\Gamma':	'$\\Gamma$',
		'Z_1':		'Z$_1$',
		'\Sigma':	'$\\Sigma$',
		'\Sigma_1':	'$\\Sigma_1$',
		'Y_1':		'Y$_1$',
		'X_1':		'X$_1$',
		'A_1':		'A$_1$',
		'C_1':		'C$_1$',
		'D_1':		'D$_1$',
		'H_1':		'H$_1$',
		'H_2':		'H$_2$',
		'L_1':		'L$_1$',
		'L_2':		'L$_2$',
		'B_1':		'B$_1$',
		'P_1':		'P$_1$',
		'P_2':		'P$_2$',
		'Q_1':		'Q$_1$',
		'M_1':		'M$_1$',
		'M_2':		'M$_2$',
		'I_1':		'I$_1$',
		'F_1':		'F$_1$'
	}
	if raw in conv.keys():
		return conv[raw]
	else:
		return raw



# ******************************** #
# SHARED DATA EXTRACTION FUNCTIONS #
# ******************************** #

# PARSE POSCAR FOR ATOMS
def atomDict(baseDir):
	''' Parse POSCAR for atoms '''
	with open(os.path.join(baseDir, 'POSCAR'), 'r') as f:
		for i, line in enumerate(f):
			if i is 5:
				atomType = line.split()
			if i is 6:
				atomMult = line.split()
				atomMult = [int(am) for am in atomMult]
				break
	atomList = []
	for at, am in zip(atomType, atomMult):
		atomList.extend([at,]*am)
	atoms = {i+1: a for i , a in enumerate(atomList)}
	nAtoms = len(atomList)
	return atoms

def fermiLevel(baseDir):
	''' Parse DOSCAR for Fermi energy '''
	Ef = None
	with open(os.path.join(baseDir, 'DOSCAR'), 'r') as f:
		for i, line in enumerate(f):
			if i is 5:
				Ef = float(line.split()[3])
				break
	return Ef

def parseOutcar(baseDir):
	''' Parse OUTCAR for spin-orbit coupling flag, reciprocal lattice vectors '''
	outFile = os.path.join(baseDir, 'OUTCAR')
	lines2read = 0
	with open(outFile, 'r') as f:
		for line in f:
			if 'LSORBIT' in line:
				if line.split()[2] is 'T':
					so = True
				else:
					so = False
			if lines2read > 0:
				# L = line.split()
				if lines2read is 3:
					# b1 = np.array([float(L[3]),float(L[4]),float(L[5])])
					b1 = np.array([float(line[45:58]), float(line[58:71]), float(line[72:])])
				elif lines2read is 2:
					# b2 = np.array([float(L[3]),float(L[4]),float(L[5])])
					b2 = np.array([float(line[45:58]), float(line[58:71]), float(line[72:])])
				elif lines2read is 1:
					# b3 = np.array([float(L[3]),float(L[4]),float(L[5])])
					b3 = np.array([float(line[45:58]), float(line[58:71]), float(line[72:])])        
				lines2read -= 1
			if 'reciprocal lattice vectors' in line:
				lines2read = 3
	return so, (b1, b2, b3)

def parseIBZKPT(baseDir):
	''' Parse IBZKPT for total number of k-points, and number in IBZ '''
	fname = os.path.join(baseDir, 'IBZKPT')
	with open(fname, 'r') as f:
		for i, line in enumerate(f):
			if i is 1:
				nkibz = int(line)
				break
	data = np.genfromtxt(fname, skip_header=3, max_rows=nkibz)
	nktot = int(np.sum(data, axis=0)[-1])
	return nkibz, nktot



# **************************************** #
# BAND STRUCTURE DATA EXTRACTION FUNCTIONS #
# **************************************** #

def parseProcar(baseDir, SOC):
	''' Parse PROCAR for (fractional) k-points, energies, occupancies '''
	raw, newFormat = [], False
	extraLine = 0
	with open(os.path.join(baseDir, 'PROCAR'), 'r') as f:
		for i, line in enumerate(f):
			raw.append(line)
			if i is 8 and line[2] is ' ':
				newFormat = True
				print('NEW PROCAR FORMAT DETECTED')
			if i is 8 and len(line.split()) > 11:
				extraLine = 1
				print('f orbitals detected')
			if i is 1:
				NK = int(line.split()[3])
				NB = int(line.split()[7])
				NION = int(line.split()[11])
				print('Found {} K-points, {} bands, {} ions.'.format(NK, NB, NION))
	if SOC:
		NLINES = 4*(NION+1+extraLine)+4
	else:
		NLINES = NION+5+extraLine
	E = np.zeros((NK,NB), dtype=float)
	K = {}
	P = np.zeros((NION,3,NK,NB), dtype=float)
	firstTimeFlag = True
	for i in range(NK):
		i_kheader = i*(NLINES*NB+3)+3
		line = raw[i_kheader]
		if newFormat:
			K[i] = np.array([float(line[19:30]), float(line[30:41]), float(line[41:52])])
		else:
			K[i] = np.array([float(line[18:29]), float(line[29:40]), float(line[40:51])])
		for j in range(NB):
			i_eheader = i_kheader + j*NLINES+2
			line = raw[i_eheader]
			E[i, j] = float(line.split()[4])
			for k in range(NION):
				line = raw[i_eheader+3+k]
				if firstTimeFlag:
					if len(line.split()) > 5:
						print('orbital projections are lm-decomposed')
						endind = 10
					else:
						print('spd projection only')
						endind = 4
					P = np.zeros((NION, endind-1, NK, NB), dtype=float)
					firstTimeFlag = False
				spd = np.array([float(thing) for thing in line.split()[1:endind]])
				P[k, :, i, j] = spd
	return K, E, P



# ******************************** #
# D.O.S. DATA EXTRACTION FUNCTIONS #
# ******************************** #

def getTotalDOS(baseDir):
	''' Parse DOSCAR for total density of states '''
	with open(os.path.join(baseDir, 'DOSCAR'), 'r') as f:
		for i, line in enumerate(f):
			if i is 5:
<<<<<<< HEAD
				Ef = float(line.split()[3])				#Fermi energy from DOSCAR header
				nebin = int(line.split()[2])			#NEDOS from DOSCAR header
				break
	data = np.genfromtxt(os.path.join(baseDir, 'DOSCAR'), skip_header=6, max_rows=nebin)
	E = data[:,0] - Ef 									#Energy from fermi energy
	tdos = data[:,1]									#Dos at each energy (not cumulative, just for that energy)
=======
				Ef = float(line.split()[3])
				nebin = int(line.split()[2])
				break
	data = np.genfromtxt(os.path.join(baseDir, 'DOSCAR'), skip_header=6, max_rows=nebin)
	E = data[:,0] - Ef
	tdos = data[:,1]
>>>>>>> newdoug
	return E, tdos

def getpDOS(baseDir, nAtoms, nebin, so):
	''' Parse DOSCAR for atom- and orbital-projected density of states '''
	pdos = np.zeros((nAtoms, 9, nebin))
	for k in range(nAtoms):
		data = np.genfromtxt(os.path.join(baseDir, 'DOSCAR'), skip_header=6+(nebin+1)*(k+1), skip_footer=(nebin+1)*(nAtoms-k-1))
		if so:
			# 36 column pdos: 9 (s, px, py, pz, ...) times 4 (m_tot, m_x, ...)
			# as written, this function only parses the m_tot values
			pdos[k, 0, :] = data[:, 1]		# s
			pdos[k, 1, :] = data[:, 5]		# p
			pdos[k, 2, :] = data[:, 9]
			pdos[k, 3, :] = data[:, 13]
			pdos[k, 4, :] = data[:, 17]		# d
			pdos[k, 5, :] = data[:, 21]
			pdos[k, 6, :] = data[:, 25]
			pdos[k, 7, :] = data[:, 29]
			pdos[k, 8, :] = data[:, 33]
		else:
			# 9 column pdos (s, px, py, pz, ...)
			pdos[k, :, :] = data[:, 1:].T
	return pdos



# ****************** #
# ANALYSIS FUNCTIONS #
# ****************** #

def scatterPoints(kd, E):
	''' Repeat kd and flatten E to make suitable for scatter plotting '''
	ee = E.flatten(order='C')
	kk = np.repeat(kd, E.shape[1])
	return kk, ee

def momentumCoord(K, ktol=0.02, shuffleVec=None):
	''' Convert k-points (3-D space) to BZ path coordinate, detect special points and discontinuities in path
	     - K:    Absolute K-points dictionary, output of absoluteKPoints()
	     - ktol: Threshold absolute k-point separation for detecting path breaks
	'''
	kdist = [0.]
	dirChange, pathBreak = [], []
	if shuffleVec is not None:
		K = {i: K[shuffleVec[i]] for i in range(len(shuffleVec))}
<<<<<<< HEAD
=======
	# kstepPrev = None
>>>>>>> newdoug
	for i in sorted(K.keys())[1:]:
		kstep = np.linalg.norm(K[i] - K[i-1])
		if kstep > ktol:									# discontinuity
			pathBreak.append(kdist[-1])
			kdist.append(kdist[-1])
		elif kstep == 0.:									# repitition indicates special point
<<<<<<< HEAD
=======
			# if kstepPrev != 0:						# handle edge case where multiple repeats occur for some reason
>>>>>>> newdoug
			dirChange.append(kdist[-1])
			kdist.append(kdist[-1])
		else:
			kdist.append(kdist[-1] + kstep)
<<<<<<< HEAD
=======
		# kstepPrev = kstep
>>>>>>> newdoug
	pathBreak = np.array(pathBreak)
	pathBreak = pathBreak[np.where(pathBreak > 0)]
	return np.array(kdist), np.array(dirChange), pathBreak

def absoluteKPoints(Kfrac, b):
	''' Convert fractional k-point coordinates to absolute coordinates using reciprocal lattice vectors '''
	K = {}
	bmat = np.vstack((b[0],b[1],b[2]))
	for i in sorted(Kfrac.keys()):
		K[i] = np.dot(Kfrac[i], bmat)
	return K

def smoothDOS(x, sigma):
	''' Smooth a DOS (or any 1-D quantity) by convolution with a Gaussian '''
	from scipy.ndimage.filters import gaussian_filter1d
	return gaussian_filter1d(x, sigma)

def mapAtoms(at):
	''' Re-map atom dictionary from parsing POSCAR to allow indexing of orbital projections '''
	atr = {s: [] for s in set(at.values())}
	for k, v in at.items():
		atr[v].append(k-1)
	return atr

def mapOrbs(projType='nlm'):
	''' Return dictionary that maps orbital labels to columns of projection outputs '''
	if projType is 'nlm':
		return {'s': [0], 'p': [1,2,3], 'd': [4,5,6,7,8], 'all': [0,1,2,3,4,5,6,7,8]}
	elif projType is 'nl':
		return {'s': [0], 'p': [1], 'd': [2], 'all': [0,1,2]}
	else:
		raise Exception('Only "n,l,m" or "n,l" quantum number projection allowed')

def bandgap(E, tdos, tol=1e-5):
	''' Return the bandgap given the total DOS as function of energy '''
	x, y = E, tdos
	vbm = x[(y>tol) & (x<0)][-1]
	cbm = x[(y>tol) & (x>0)][0]
	return cbm-vbm



# ****************** #
# PLOTTING FUNCTIONS #
# ****************** #

def labelBZ(ax, kdist, dirChange, pathBreak, labels, fdy=0.05):
	''' Label high symmetry points of the Brillouin zone along momentum axis '''
	i = 0
	yl = ax.get_ylim()
	joined = np.sort(np.concatenate((np.array([np.min(kdist)]), dirChange, pathBreak, np.array([np.max(kdist)]))))
	for j, v in enumerate(joined):
		if v in pathBreak:
			s = '|'.join((labels.pop(0), labels.pop(0)))
		else:
			s = labels.pop(0)
		ax.text(v, yl[0] - (yl[1] - yl[0])*fdy, s, horizontalalignment='center')
	if len(labels)>0:
		print('\n******************************************************')
		print('***** WARNING: LABEL ENTRIES REMAIN, CHANGE KTOL *****')
		print('******************************************************')

def bandPlot(ax, kdist, E, dirChange, pathBreak, bzl=['',]*100, yl=None, fdy=0.03, gray=0.5):
	''' Plot band structure, no orbital projections '''
	if yl is None:
		yl = np.min(E)-1, np.max(E)+1
	else:
		yl = yl[0], yl[1]
	xmin, xmax = np.min(kdist), np.max(kdist)
	ax.plot( kdist , E , '-' , c=(gray,)*3 , lw=3 )
	for k in dirChange:											# plot vertical lines
		ax.plot((k, k), yl, '-', c=(gray,)*3, lw=2, zorder=0)
	for k in pathBreak:
		ax.plot((k, k), yl, 'k-', lw=2, zorder=3)
	ax.plot((xmin, xmax), (0, 0), 'k--', lw=2)				# plot Fermi level
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(yl)
	ax.set_ylabel('$E-E_f$ (eV)')
	ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')	# remove x ticks and ticklabels
	labelBZ(ax, kdist, dirChange, pathBreak, bzl, fdy)
	ax.grid(b='off')
	return



# ********** #
# HIGH-LEVEL #
# ********** #

def spaghetti(bd, ct, ax=None, El=None, orbs=(), ms=200, vbmRef=True, lblOffset=0.03, shuffleVec=None):
	''' Parse band structure calculation and plot
		'ax'  : Destination axes for plot
		'bd'  : Directory with VASP output files
		'ct'  : Label describing Brillouin zone type (e.g. FCC, BCT1) to use lookup table for special point labels
		'El'  : (optional) energy limits for plot
		'orbs': (optional) list or tuple of orbital projections to plot specified as a tuple of (ATOM, ORBITAL)
		        - e.g. [ ('Pt', 'd'), ('O', 'all') ] will plot two orbital projections, the Pt d-states, and all O states
				- currently only orbital options are 's', 'p', 'd', 'all', but woud be easily expanded to do m-quantum number-resolved
				- currently can only sum all atoms of same type, but would be easily expanded to plot only certain atoms of a type
		'ms'  : (optional) marker size for orbital projections
	'''
	so, b = parseOutcar(bd)
	Ef = fermiLevel(bd)
	Kf, E, P = parseProcar(bd, so)
	if vbmRef:
		E = E-np.max(E[E<Ef])
	else:
		E = E-Ef
	K = absoluteKPoints(Kf, b)
	kd, dc, pb = momentumCoord(K, shuffleVec=shuffleVec)
	if shuffleVec is not None:
		E = E[shuffleVec, :]
	if ax:
		bzl = [labelify(bzli) for bzli in bzLabels[ct]]
		bandPlot(ax, kd, E, dc, pb, bzl=bzl, yl=El, fdy=lblOffset)
		at = atomDict(bd)
		atr = mapAtoms(at)
		if P.shape[1] is 9:
			pt = 'nlm'
		else:
			pt = 'nl'
		om = mapOrbs(projType=pt)
		kk, ee = scatterPoints(kd, E)
		for i, (a, o) in enumerate(orbs):
			ax.scatter(kk, ee, s=np.sum(np.sum(P[atr[a], :, :, :], axis=0)[om[o], :, :], axis=0)*ms, zorder=3, facecolor=cList[i], lw=0.5, alpha=0.7)
		if vbmRef:
			ax.set_ylabel('$E-E_{VBM}$ (eV)')
		if El is not None:
			ax.set_ylim(El)
	return kd, E, K, P, pb

def spiky(ax, bd, El=None, orbs=(), smooth=0):
	''' Parse density of states calculation and plot
		'ax'    : Destination axes for plot
		'bd'    : Directory with VASP output files
		'El'    : (optional) energy limits for plot
		'orbs'  : (optional) list or tuple of orbital projections to plot specified as a tuple of (ATOM, ORBITAL)
		        - e.g. [ ('Pt', 'd'), ('O', 'all') ] will plot two orbital projections, the Pt d-states, and all O states
				- currently only orbital options are 's', 'p', 'd', 'all', but woud be easily expanded to do m-quantum number-resolved
				- currently can only sum all atoms of same type, but would be easily expanded to plot only certain atoms of a type
		'smooth': (optional) gaussian convolution smoothing factor for visual clarity of DOS. Units are the energy resolution in the DOSCAR.
		        - 0 means no smoothing
		        - does not affect calculation of parameters from the DOS (e.g. bandgap), only used in plotting
	'''
	at = atomDict(bd)
	so, _ = parseOutcar(bd)
	E, tdos = getTotalDOS(bd)
	pdos = getpDOS(bd, len(at), E.size, so)
	ax.plot(smoothDOS(tdos, sigma=smooth), E, 'k-', label='total', zorder=3)
	atr = mapAtoms(at)
	if pdos.shape[1] is 9:
		pt = 'nlm'
	else:
		pt = 'nl'
	om = mapOrbs(projType=pt)
	for a, o in orbs:
		ax.plot(smoothDOS(np.sum(np.sum(pdos[atr[a], :, :], axis=0)[om[o], :], axis=0), sigma=smooth), E, label='{} ${}$'.format(a, o))
	if El is not None:
		ax.set_ylim(El)
	ax.plot(ax.get_xlim(), (0, 0), 'k--')
	ax.legend()
	ax.set_ylabel('$E-E_{VBM}$ (eV)')
	ax.set_xlabel('DOS')
	ax.grid(b='off')
	ax.set_xticks([])
	return E, tdos, pdos

# fin.
