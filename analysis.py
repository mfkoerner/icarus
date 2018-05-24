

import os
from electrons import np, plt, minDirectTransition
import photons as ph
from matgen import td, ad, getFormulae, getFormula, mpSpaghetti, ed, mpr, unpickle
from pymatgen.analysis import bond_valence as bv
import matplotlib.gridspec as gs
from pymatgen.core.periodic_table import Element, Specie

bd = os.path.join(td, 'General-Direct')

refList = [
	# IV, III-V, II-VI
	'mp-149',		# Si
	'mp-2534',		# GaAs
	'mp-406',		# CdTe
	# 'mp-2133',		# ZnO
	# CIGS-like
	'mp-5238',		# CuGaS2
	'mp-22736',		# CuInS2
	# CZTS-like
	# 'mp-6408',		# Cu2ZnGeS4
	#halide perovskites
	# 'mp-567629',	# o-CsPbBr3
	# 'mp-23037',		# CsPbCl3
	# 'mp-CABB',		# CsAgBiBr6 (not in MP)
	# emerging chalcogenides
	'mp-2231',		# o-SnS
	# 'mp-3102',		# CuTaS3
	'mp-2160',		# Sb2Se3
	'mp-4468',		# CuSbS2
	'mp-5414'		# NaSbS2
	]

anions = set(bv.ELECTRONEG)
radioactive = set([Element('Po'), Element('Fr'), Element('Tc'), Element('Ac')])
toxic = set([Element('Tl'), Element('Hg')])
pnictogens = set([Element('N'), Element('P'), Element('As'), Element('Sb')])
chalcogens = set([Element('S'), Element('Se'), Element('Te')])
halogens = set([Element('F'), Element('Cl'), Element('Br'), Element('I')])
oxygen = set([Element('O')])
hydrogen = set([Element('H')])
mainGroup = anions - set([Element('H')]) | set([Element('Al'), Element('Ga'), Element('Ge'), Element('In'), Element('Sn'), Element('Tl'), Element('Pb'), Element('Bi')])
alkali = set([Element('Li'), Element('Na'), Element('K'), Element('Rb'), Element('Cs')])
alkEarth = set([Element('Be'), Element('Mg'), Element('Ca'), Element('Sr'), Element('Ba')])
II = set([Element('Zn'), Element('Cd'), Element('Hg')])
III = set([Element('B'), Element('Al'), Element('Ga'), Element('In')])
IV = set([Element('Ge'), Element('Sn'), Element('Pb')])
V = pnictogens
VI = chalcogens | set([Element('O')])

def Jmax(mpid):
	''' '''
	hnu, alpha, eps, N, Eg = ph.optical(os.path.join(ad, 'Absorption', mpid, 'absorb'))
	xx, yySun, yyMat, ss = ph.FOM(hnu, alpha, Eg)
	return ss[-1]/1e8, Eg, hnu, alpha

def getClassif():
	''' '''
	raw = []
	with open(os.path.join(ad, 'classification-PRELIM.csv'), 'r') as fh:
		for line in fh:
			raw.append(line.strip('\n'))
	_ = raw.pop(0)
	classif = {r.split(',')[0]: r.split(',')[2] for r in raw}
	out = {}
	for family in set(classif.values()):
		out[family] = set([k for k in classif.keys() if classif[k] == family])
	return out

def plot2Dannotated(K, X, Y, annotate=True):
	''' '''
	fig, ax = plt.subplots()
	ax.plot([X[k] for k in K], [Y[k] for k in K], '.')
	form = getFormulae()
	if annotate:
		for k in K:
			ax.text(X[k], Y[k], form[k], fontdict={'size': 10})
	return fig, ax

def plotFamilies(data, xdata='egopt', ydata='jmax'):
	''' '''
	families = ['III-V', 'II-VI', 'IV-VI', 'pnictide', 'chalcogenide', 'halide', 'oxide', 'heusler', 'zintl']
	lims = {'egopt': (0, 2.5), 'jmax': (0, 14), 'jmax*Eg': (0, 7), 'me_mean': (3e-2, 1e0), 'mh_mean': (3e-2, 1e0)}
	scale = {'egopt': 'linear', 'jmax': 'linear', 'jmax*Eg': 'linear', 'me_mean': 'log', 'mh_mean': 'log'}
	ytextmult = {'jmax': 0.9, 'me_mean': 0.06, 'mh_mean': 0.06, 'jmax*Eg': 0.9}
	ticks = {'egopt': [0, 1, 2], 'jmax*Eg': [0, 2, 4, 6], 'jmax': [0, 4, 8, 12]}
	labels = {'egopt': '$E_g$ (eV)', 'jmax': '$j_{max}$ (a.u.)', 'me_mean': 'avg($m_e^*$) / $m_0$', 'mh_mean': 'avg($m_h^*$) / $m_0$', 'jmax*Eg': '$j_{max} \\times E_g$ (a.u.)'}
	K = list(data[xdata].keys())
	x = {k: data[xdata][k] for k in K}
	if ydata == 'jmax*Eg':
		y = {k: data['egopt'][k]*data['jmax'][k] for k in K}
	else:
		y = {k: data[ydata][k] for k in K}
	raw = []
	with open(os.path.join(ad, 'classification-PRELIM.csv'), 'r') as fh:
		for line in fh:
			raw.append(line.strip('\n'))
	classif = {r.split(',')[0]: r.split(',')[2] for r in raw[1:]}
	fig, ax = plt.subplots()
	grid = gs.GridSpec(3, 3, wspace=0.02, hspace=0.02)
	axList = []
	for i, fam in enumerate(families):
		ax = plt.subplot(grid[i%3, int(i/3)])
		ax.plot([x[k] for k in K], [y[k] for k in K], '.', c=(0.7,)*3)
		ax.plot([x[k] for k in K if classif[k]==fam], [y[k] for k in K if classif[k]==fam], '.')
		ax.grid('off')
		ax.set_xlim(lims[xdata])
		ax.set_xscale(scale[xdata])
		if xdata in ticks.keys():
			ax.set_xticks(ticks[xdata])
		if i%3 != 2:
			ax.set_xticklabels('')
		ax.set_ylim(lims[ydata])
		ax.set_yscale(scale[ydata])
		if ydata in ticks.keys():
			ax.set_yticks(ticks[ydata])
		if int(i/3) != 0:
			ax.set_yticklabels('')
		ax.text(lims[xdata][-1]*0.95, lims[ydata][-1]*ytextmult[ydata], fam, ha='right', va='top')
		axList.append(ax)
		if i is 1:
			ax.set_ylabel(labels[ydata])
		elif i is 5:
			ax.set_xlabel(labels[xdata])
	return fig, axList

def xtalSummary(data, mpid):
	''' '''
	fig = plt.figure(figsize=(20, 10))
	grid = gs.GridSpec(3, 3, width_ratios=(2, 2, 1), height_ratios=(1, 1, 2))
	axL = plt.subplot(grid[:, 0])
	axM1 = plt.subplot(grid[0, 1])
	axM2 = plt.subplot(grid[1, 1])
	axM3 = plt.subplot(grid[2, 1])
	axR1 = plt.subplot(grid[0:2, 2])

	# band structure
	kd, E, K, _, pb = mpSpaghetti(mpid, ax=axL, El=(-4, 6), lblOffset=0.04)
	ed.masses(kd, E, K, pb, ax=axL)
	_ = minDirectTransition(E, kd, ax=axL)



	# absorption
	hnu, alpha, _, _, egopt = ph.optical(os.path.join(ad, 'Absorption', mpid, 'absorb'))
	ph.plotAbsorption(axM3, hnu[hnu>egopt], alpha[hnu>egopt])
	hnu, alpha, _, _, egopt = ph.optical(os.path.join(ad, 'Absorption', 'mp-149', 'absorb'))
	axM3.plot(hnu[hnu>egopt], alpha[hnu>egopt], c=(0.7,)*3, zorder=0, label='Si')
	axM3.text(0.65, 8e2, 'Si', va='top', ha='center')

	# summary
	xl = axR1.get_xlim()
	yl = axR1.get_ylim()
	axR1.text(xl[0], yl[0]+(yl[1]-yl[0])/2, '{} ({})\nS.G. #{}\n{} meV/atom above hull\n$E_g$ = {:.2f} eV\n$E_{{direct}}$ = {:.2f} eV\navg($m_e^*$) = {:.2f} $m_0$\navg($m_h^*$) = {:.2f} $m_0$\n$j_{{max}}$ = {:.2f}'.format(
		getFormula(mpid), mpid,
		mpr.query(mpid, properties=['spacegroup'])[0]['spacegroup']['number'],
		int(data['eh'][mpid]*1000),
		data['eg'][mpid], data['dg'][mpid],
		data['me_mean'][mpid], data['mh_mean'][mpid],
		data['jmax'][mpid]
		), va='center')
	axR1.axis('off')
	return fig




# CANNOT CURRENTLY DETECT SPECIAL ANIONS (e.g. O2(2-), PO4(3-), SO4(2-), etc.)
# CANNOT CURRENTLY DETECT MAIN GROUP - MAIN GROUP BONDING IN ZINTL PHASES, ETC.
def classifyCompound(mpid):
	''' '''
	bva = bv.BVAnalyzer()
	s = unpickle('E:/Dropbox/Materials/Research/Project-Icarus/2018-03-23-structures-n627')[mpid]
	try:
		ss = bva.get_oxi_state_decorated_structure(s)
	except ValueError as err:
		return 'VALUE ERROR: {}'.format(err)
	c = ss.composition
	elem = set([ci.element for ci in c])
	a = elem & anions
	# no anions
	if len(a) is 0:
		if elem <= (alkali | alkEarth | mainGroup):
			return 'zintl'
		return 'intermetallic'
	# anions
	if a <= pnictogens:
		return 'pnictide'
	elif a <= chalcogens:
		return 'chalcogenide'
	elif a <= halogens:
		return 'halide'
	elif a <= oxygen:
		return 'oxide'
	# mixed anion families
	elif a <= (oxygen | hydrogen):
		return 'hydroxide'
	elif a <= (oxygen | pnictogens):
		return 'oxypnictide'
	elif a <= (oxygen | chalcogens):
		return 'oxychalcogenide'
	elif a <= (oxygen | halogens):
		return 'oxyhalide'
	elif a <= (chalcogens | halogens):
		return 'chalcohalide'
	return 'UNCLASSIFIABLE'


def loadCompiled():
	''' '''
	data = unpickle(os.path.join(ad, '2018-03-19-MASTER-n627'))
	refs = unpickle(os.path.join(ad, '2018-03-21-REF-optical-n14'))
	return data, refs

def plotAbsSummary(data, refs=None):
	''' '''
	fig, ax = plt.subplots()
	# grid = gs.GridSpec(2, 1, hspace=0.0)
	# ax1 = plt.subplot(grid[i%3, int(i/3)])
	K = list(data['hnu'].keys())
	for k in K:
		ax.semilogy(data['hnu'][k][data['hnu'][k]>data['egopt'][k]], data['alpha'][k][data['hnu'][k]>data['egopt'][k]], c=(0.7,)*3)
	h = []
	if refs:
		for k in list(refs['hnu'].keys()):
			h.append(ax.semilogy(refs['hnu'][k][refs['hnu'][k]>refs['egopt'][k]], refs['alpha'][k][refs['hnu'][k]>refs['egopt'][k]], lw=3, label=getFormula(k))[0])
	ax.legend(loc='lower right')
	ax.set_xlim(0, 4)
	ax.set_xticks([0,1,2,3,4])
	ax.set_ylim(1e3, 3e6)
	ax.set_xlabel('$h\\nu$ (eV)')
	ax.set_ylabel('$\\alpha$ (cm$^{-1}$)')
	return h

def loadSQ():
	''' '''
	xy = np.loadtxt('data/SQ.dat')
	return xy[:, 0], xy[:, 1]
