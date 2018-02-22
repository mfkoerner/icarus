########################################
# written for Python 3                 #
# by Doug Fabini (fabini@mrl.ucsb.edu) #
########################################

'''

Dependencies:
 - electrons (in-house module)
 - os
 - numpy
 - matplotlib
 - pymatgen
 - json

'''



# IMPORT PACKAGES
import config
import electrons as el, edges as ed
import os, json, pickle
np = el.np
plt = el.plt
from pymatgen import MPRester, Lattice, Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.ext.matproj import MPRestError
from json import JSONDecodeError

mpr = MPRester(config.matprojapi)

td = config.shared_loc

ad = config.local_archive



# ******************************************* #
# QUICK FUNCTIONS FOR COMMON DATABASE QUERIES #
# ******************************************* #

def nEntries():
	''' Count number of entries in MP database '''
	return len(mpr.query(criteria={}, properties=['material_id']))

def nMetals():
	''' Count number of metals in MP database '''
	return len(mpr.query(criteria={'band_gap': {'$eq': 0}}, properties=['material_id']))

def nInsulators():
	''' Count number of insulators in MP database '''
	return len(mpr.query(criteria={'band_gap': {'$gt': 0}}, properties=['material_id']))

def nEntriesWithBS():
	''' Count number of entries with computed bandstructures in MP database '''
	return len(mpr.query(criteria={'has_bandstructure': True}, properties=['material_id']))

def nMetalsWithBS():
	''' Count number of metals with computed bandstructures in MP database '''
	return len(mpr.query(criteria={'has_bandstructure': True, 'band_gap': {'$eq': 0}}, properties=['material_id']))

def nInsulatorsWithBS():
	''' Count number of insulators with computed bandstructures in MP database '''
	return len(getInsulatorsWithBS())

def getInsulatorsWithBS():
	''' Return material_id's of all insulators with computed bandstructures in MP database '''
	return mpr.query(criteria={'has_bandstructure': True, 'band_gap': {'$gt': 0}}, properties=['material_id', 'pretty_formula'])

def nSemiconductorsWithBS():
	''' Count number of insulators with computed bandstructures in MP database '''
	return len(getSemiconductorsWithBS())

def getSemiconductorsWithBS():
	''' Return material_id's of all semiconductors (Eg < 2.5 eV) with computed bandstructures in MP database '''
	return mpr.query(criteria={'has_bandstructure': True, 'band_gap': {'$gt': 0, '$lte': 2.5}}, properties=['material_id', 'pretty_formula', 'band_gap', 'e_above_hull'])

def getSpaceGroup(mpid):
	''' Get space group symbol and number for entry with material_id "mpid" '''
	return mpr.query(criteria={'material_id': mpid}, properties=['structure'])[0]['structure'].get_space_group_info()

def getFormula(mpid):
	''' Get pretty formula for entry with material_id "mpid" '''
	return mpr.query(criteria={'material_id': mpid}, properties=['pretty_formula'])[0]['pretty_formula']

def getRealLatt(mpid):
	''' Get the real lattice for a single material from the 'structure' object'''
	return mpr.get_structure_by_material_id(mpid).lattice

# *****************************************************************************************************************************************************************************************
# POTENTIALLY NEEDS FIXING FOR BZ SIZE ISSUE **********************************************************************************************************************************************
def getGapMomentum(mpid):
	''' '''
	bs = mpr.get_bandstructure_by_material_id(mpid)
	return np.linalg.norm(bs.get_cbm()['kpoint'].cart_coords - bs.get_vbm()['kpoint'].cart_coords)



# ******************* #
# LOW LEVEL UTILITIES #
# ******************* #

def removeAdjacent(a):
	''' Remove adjacent duplicates from list '''
	i = 1
	while i < len(a):    
		if a[i] == a[i-1]:
			a.pop(i)
			i -= 1  
		i += 1
	return a

def jsonLoad(f):
	''' Load json-serialized data from filepath f '''
	with open(f, 'r') as fh:
		return json.load(fh)

def unpickle(f):
	''' Load pickle-serialized data from filepath f '''
	with open(f, 'rb') as fh:
		return pickle.load(fh)


# ************************************ #
# FLEXIBLE (LOCAL / DATABASE) QUERYING #
# ************************************ #

def getBSRecipLatt(mpid, local=True):
	'''
	Return reciprocal lattice associated with the band structure calculation for a given material id
	 - if local, get from local archive
	 - else, get from database
	'''
	if local:
		return Lattice.from_dict(jsonLoad(os.path.join(td, 'autoparse/BS-Recip-Latt/{}'.format(mpid))))
	else:
		return latticeFromTaskID(getBSTaskID(mpid)[mpid]).reciprocal_lattice

def getBS(mpid, local=True):
	'''
	Return band structure calculation for a given material id
	 - if local, get from local archive
	 - else, get from database
	'''
	if local:
		return unpickle(os.path.join(ad, 'BS/{}'.format(mpid)))
	else:
		return mpr.get_bandstructure_by_material_id(mpid)



# *************************** #
# FUNCTIONS FOR BZ WORKAROUND #
# *************************** #

def getBSTaskID(mpids):
	''' '''
	if type(mpids) is str:
		mpids = [mpids]
	data = mpr.query({'task_id': {'$in': mpids}}, ['material_id', 'band_structure'], mp_decode=False)
	out = {}
	for d in data:
		try:
			if 'GGA+U' in d['band_structure'].keys():
				rt = 'GGA+U'
				out[d['material_id']] = d['band_structure'][rt]['task_id']
			elif 'GGA' in d['band_structure'].keys():
				rt = 'GGA'
				out[d['material_id']] = d['band_structure'][rt]['task_id']
			else:
				out[d['material_id']] = None
		except:
			out[d['material_id']] = None
	return out

def latticeFromTaskID(tid):
	''' Construct and return the real lattice associated with task ID tid '''
	td = mpr.get_task_data(tid)
	realLatt = Structure.from_str(td[0]['cif'], fmt='cif').lattice
	return realLatt

def getRunTypes(mpids):
	''' Returns run type(s) (e.g. 'GGA') for all materials in mpids in dictionary format {mpid: set(run types)}'''
	if type(mpids) is str:
		mpids = [mpids]
	data = mpr.query({'task_id': {'$in': mpids}}, ['material_id', 'band_structure'], mp_decode=False)
	out = {}
	for d in data:
		try:
			out[d['material_id']] = set(d['band_structure'].keys())
		except:
			out[d['material_id']] = []
	return out








# ************************************************ #
# HIGH-LEVEL DATA GRABBING / VISUALIZING FUNCTIONS #
# ************************************************ #

def getMPbands(mpid, scissor=None, workaround=True, local=True, warnings=False):
	''' '''
	# (Optionally) print code limitations warning
	if warnings:
		print('\n***** CURRENTLY ONLY GETS SPIN UP BANDS *****')
	# Load band structure data and (optionally) apply scissor
	bs = getBS(mpid, local=local)
	if scissor:
		bs = bs.apply_scissor(scissor)
		if bs is None:
			raise Exception('Scissor operation failed (metallic?)')
	# Reference energy depending on metal / insulator character
	if bs.is_metal():
		E = bs.bands[Spin.up].T-bs.efermi	# metallic, reference energies to Fermi energy
	else:
		E = bs.bands[Spin.up].T - bs.get_vbm()['energy']			# insulating, reference energies to VBM
	# Transform k-points appropriately, optionally applying BZ size issue workaround
	if workaround:
		# WORKAROUND: LOAD RECIPROCAL LATTICE EXPLICITLY FROM BS TASK ID (STILL SOMETIMES WRONG)
		RLmat = getBSRecipLatt(mpid, local=local).matrix
		K = {i: np.dot(bs.kpoints[i].frac_coords, RLmat) for i in range(len(bs.kpoints))}
	else:
		K = {i: bs.kpoints[i].cart_coords for i in range(len(bs.kpoints))}
	# Transform high-symmetry point labels
	br = [bs.branches[i]['name'] for i in range(len(bs.branches))]
	bzl = []
	for b in br:
		bzl.extend(b.split('-'))
	bzl = [el.labelify(bzli) for bzli in removeAdjacent(bzl)]
	return E, K, bzl

def mpSpaghetti(mpid, ax=None, El=None, lblOffset=0.03, ktol=1e-1, scissor=None, workaround=True, local=True):
	''' '''
	E, K, bzl = getMPbands(mpid, scissor=scissor, workaround=workaround, local=local)
	kd, dc, pb = el.momentumCoord(K, ktol=ktol)
	if ax:
		el.bandPlot(ax, kd, E, dc, pb, bzl=bzl, yl=El, fdy=lblOffset)
		ax.set_ylabel('$E-E_{VBM}$ (eV)')
		if El is not None:
			ax.set_ylim(El)
	return kd, E, K, dc, pb

def showRandom(A, ax=None, ed=False, ktol=1e-1):
	''' (Optionally) plots random band structure from input list and fits effective masses '''
	ii = int(np.random.rand()*len(A))
	mpid = A[ii]['material_id']
	name = A[ii]['pretty_formula']
	if ax:
		kd, E, K, _, pb = mpSpaghetti(mpid, El=(-4, 8), ktol=ktol, ax=ax)
	print(name)
	if ed:
		import edges as ed
		ed.masses(kd, E, K, pb, ax=ax)










# *************************** #
# HIGH-LEVEL DATABASE PARSING #
# *************************** #

def parseRealRecipLatts(A):
	''' '''
	dirName = 'Real-Recip-Latt'
	import time
	t0 = time.time()
	iSucc = 0
	iFail = 0
	print('Attempting to parse for a total of {} entries.'.format(len(A)))
	raw = []
	if 'COMPLETE.txt' in os.listdir(os.path.join(td, 'autoparse', dirName)):
		with open(os.path.join(td, 'autoparse', dirName, 'COMPLETE.txt'), 'r') as f:
			for line in f:
				raw.append(line.strip('\n'))
		print('Found list of {} completed phases.'.format(len(raw)))
	else:
		print('No list of completed phases detected.')
	for Ai in A:
		mpid = Ai['material_id']
		if mpid not in raw:
			try:
				recipBS = mpr.get_bandstructure_by_material_id(mpid).lattice_rec
				realStruc = mpr.get_structure_by_material_id(mpid).lattice
				recipStruc = realStruc.reciprocal_lattice

				with open(os.path.join(td, 'autoparse', dirName, 'Recip-BS/{}'.format(mpid)), 'w') as f:
					json.dump(recipBS.as_dict(1), f)
				with open(os.path.join(td, 'autoparse', dirName, 'Recip-Struc/{}'.format(mpid)), 'w') as f:
					json.dump(recipStruc.as_dict(1), f)
				with open(os.path.join(td, 'autoparse', dirName, 'Real-Struc/{}'.format(mpid)), 'w') as f:
					json.dump(realStruc.as_dict(1), f)

				print(mpid)
				with open(os.path.join(td, 'autoparse', dirName, 'COMPLETE.txt'), 'a') as f:
					f.write(mpid + '\n')
				iSucc += 1
			except:
				print('FAILED FOR {}'.format(mpid))
				with open(os.path.join(td, 'autoparse', dirName, 'FAILED.txt'), 'a') as f:
					f.write(mpid + '\n')
					iFail += 1
	print('Parsed {} phases ({} successful, {} failed) in {} minutes'.format(iSucc+iFail, iSucc, iFail, (time.time()-t0)/60))






def parseProp(A, prop):
	''' General database parsing function '''
	td = {
		'BS': os.path.join(ad, 'BS')
	}[prop]
	import time
	t0 = time.time()
	iSucc, iFail = 0, 0
	print('Attempting to parse for a total of {} entries.'.format(len(A)))
	raw, failedRaw = [], []
	if 'COMPLETE.txt' in os.listdir(td):
		with open(os.path.join(td, 'COMPLETE.txt'), 'r') as f:
			for line in f:
				raw.append(line.strip('\n'))
		print('Found list of {} completed phases.'.format(len(raw)))
	else:
		print('No list of completed phases detected.')
	if 'FAILED.txt' in os.listdir(td):
		with open(os.path.join(td, 'FAILED.txt'), 'r') as f:
			for line in f:
				failedRaw.append(line.strip('\n'))
		print('Found list of {} failed phases.'.format(len(failedRaw)))
	else:
		print('No list of completed phases detected.')
	for Ai in A:
		mpid = Ai['material_id']
		if mpid not in raw:
			try:
				if prop is 'BS':
					data = mpr.get_bandstructure_by_material_id(mpid)
					dd = data.as_dict()
					for s in ('vbm', 'cbm'):
						if len(dd[s]['kpoint_index']) is 1:
							print(s)
				else:
					raise Exception('Property {} not implemented yet'.format(prop))
				
				with open(os.path.join(td, mpid), 'w') as f:
					json.dump(dd, f)

				with open(os.path.join(td, 'COMPLETE.txt'), 'a') as f:
					f.write(mpid + '\n')
				iSucc += 1
			except (TypeError, MPRestError) as e:
				recordFailure(mpid, td, failedRaw=failedRaw)
				print(e)
				iFail += 1
	print('Parsed {} phases ({} successful, {} failed) in {} minutes'.format(iSucc+iFail, iSucc, iFail, (time.time()-t0)/60))

def recordFailure(mpid, bd, failedRaw=[]):
	''' Record failed parsing in logfile if entry not already present '''
	print('FAILED FOR {}'.format(mpid))
	if mpid not in failedRaw:
		with open(os.path.join(bd, 'FAILED.txt'), 'a') as f:
			f.write(mpid + '\n')








# **************** #
# SCRIPTY ANALYSES #
# **************** #

def masterDirectList():
	''' Returns effective mass, stability, and bandgap information for ~4000 direct gap semiconductors (Eg < 2.5 eV)'''
	import pandas as pd
	d1 = pd.read_csv(os.path.join(td, 'autoparse/Eff-Mass/MASS.txt'), index_col=0)
	d2 = pd.read_csv(os.path.join(td, 'autoparse/Gap-and-Type/GAP.txt'), index_col=0)
	d2 = d2.drop('formula', axis=1)
	A = getSemiconductorsWithBS()
	Ehull = {Ai['material_id']: Ai['e_above_hull'] for Ai in A}
	d3 = pd.DataFrame.from_dict(Ehull, orient='index')
	d3 = d3.rename(columns={0: 'E_hull'})
	d = d1.merge(d2, how='inner', left_index=True, right_index=True)
	d = d.merge(d3, how='inner', left_index=True, right_index=True)
	return d

def masterIndirectList():
	''' Returns effective mass, stability, gap momentum, and bandgap information for ~??? indirect gap semiconductors (Eg < 2.5 eV)'''
	import pandas as pd
	d1 = pd.read_csv(os.path.join(td, 'autoparse/Eff-Mass/MASS.txt'), index_col=0)
	d2 = pd.read_csv(os.path.join(td, 'autoparse/Gap-and-Type/GAP.txt'), index_col=0)
	d2 = d2.drop('formula', axis=1)
	d3 = pd.read_csv(os.path.join(td, 'autoparse/Gap-Momentum/MOMENTUM.txt'), index_col=0)
	d3 = d3.drop('formula', axis=1)
	A = getSemiconductorsWithBS()
	Ehull = {Ai['material_id']: Ai['e_above_hull'] for Ai in A}
	d4 = pd.DataFrame.from_dict(Ehull, orient='index')
	d4 = d4.rename(columns={0: 'E_hull'})
	d = d1.merge(d2, how='inner', left_index=True, right_index=True)
	d = d.merge(d3, how='inner', left_index=True, right_index=True)
	d = d.merge(d4, how='inner', left_index=True, right_index=True)
	return d

def calcMasses(mpids, ktol=1e-1):
	''' '''
	for j, mpid in enumerate(mpids):
		try:
			E, K, _ = getMPbands(mpid)
			print('\n############################# {} / {} #############################'.format(j+1, len(mpids)))
			kd, _, pb = el.momentumCoord(K, ktol=ktol)
			me, mh = ed.masses(kd, E, K, pb, ax=None)
			data = {
				'me_min': np.min(me), 'me_max': np.max(me), 'me_mean': np.mean(me),
				'mh_min': np.min(mh), 'mh_max': np.max(mh), 'mh_mean': np.mean(mh)
				}
			with open(os.path.join(ad, 'MASS/{}'.format(mpid)), 'w') as fh:
				json.dump(data, fh)
		except:
			print('FAILED TO CALCULATE FOR {}'.format(mpid))

# fin.

















def getRecipLatt(mpid, bs_stored=True, struc_stored=True, struc_calc=True):
	''' Get the reciprocal lattice for a single material by any of three methods:
	 - directly from the 'bandstructure' object ("bs_stored")
	 - directly from the 'structure' object ("struc_stored")
	 - calculated from the real lattice in the 'structure' object ("struc_calc")
	'''
	out = []
	if bs_stored:
		out.append(mpr.get_bandstructure_by_material_id(mpid).lattice_rec)
	else:
		out.append(None)
	if struc_stored or struc_calc:
		s = mpr.get_structure_by_material_id(mpid).lattice
	if struc_stored:
		out.append(s.reciprocal_lattice)
	else:
		out.append(None)
	if struc_calc:
		aMat = s.matrix
		out.append(Lattice(recipLatt([aMat[0,:], aMat[1,:], aMat[2,:]])))
	else:
		out.append(None)
	return out




# *************************************** #
# OBSOLETE PARSING (USE LOCAL COPY OF BS) #
# *************************************** #

def parseEffMasses(mpids):
	''' '''
	dirName = 'Eff-Mass'
	import edges as ed
	import time
	t0 = time.time()
	iSucc = 0
	iFail = 0
	print('Attempting to parse for a total of {} entries.'.format(len(mpids)))
	raw = []
	if 'COMPLETE.txt' in os.listdir(os.path.join(td, 'autoparse', dirName)):
		with open(os.path.join(td, 'autoparse', dirName, 'COMPLETE.txt'), 'r') as f:
			for line in f:
				raw.append(line.strip('\n'))
		print('Found list of {} completed phases.'.format(len(raw)))
	else:
		print('No list of completed phases detected.')
	for mpid in mpids:
		form = getFormula(mpid)
		if mpid not in raw:
			try:
				kd, E, K, _, pb = mpSpaghetti(mpid, ax=None, ktol=1e-2)
				me, mh = ed.masses(kd, E, K, pb, ax=None)

				with open(os.path.join(td, 'autoparse', dirName, 'MASS.txt'), 'a') as f:
					f.write('{},{},{},{},{},{},{},{}\n'.format(mpid, form, np.mean(me), min(me), max(me), np.mean(mh), min(mh), max(mh)))

				print('{:>20} {:>20} {:>6.2f} {:>6.2f}'.format(mpid, form, np.mean(me), np.mean(mh)))
				with open(os.path.join(td, 'autoparse', dirName, 'COMPLETE.txt'), 'a') as f:
					f.write(mpid + '\n')
				iSucc += 1
			except:
				print('FAILED FOR {}'.format(mpid))
				with open(os.path.join(td, 'autoparse', dirName, 'FAILED.txt'), 'a') as f:
					f.write(mpid + '\n')
					iFail += 1
	print('Parsed {} phases ({} successful, {} failed) in {} minutes'.format(iSucc+iFail, iSucc, iFail, (time.time()-t0)/60))

def parseGaps(A):
	''' '''
	dirName = 'Gap-and-Type'
	import time
	t0 = time.time()
	iSucc = 0
	iFail = 0
	print('Attempting to parse for a total of {} entries.'.format(len(A)))
	raw = []
	if 'COMPLETE.txt' in os.listdir(os.path.join(td, 'autoparse', dirName)):
		with open(os.path.join(td, 'autoparse', dirName, 'COMPLETE.txt'), 'r') as f:
			for line in f:
				raw.append(line.strip('\n'))
		print('Found list of {} completed phases.'.format(len(raw)))
	else:
		print('No list of completed phases detected.')
	for Ai in A:
		mpid = Ai['material_id']
		form = Ai['pretty_formula']
		if mpid not in raw:
			try:
				bs = mpr.get_bandstructure_by_material_id(mpid)
				Eg = bs.get_band_gap()['energy']
				Egdir = bs.get_direct_band_gap()
				isDir = bs.get_band_gap()['direct']

				with open(os.path.join(td, 'autoparse', dirName, 'GAP.txt'), 'a') as f:
					f.write('{},{},{},{},{}\n'.format(mpid, form, isDir, Eg, Egdir))

				print('{:>20} {:>20} {:>5} {:>20} {:>20}'.format(mpid, form, isDir, Eg, Egdir))
				with open(os.path.join(td, 'autoparse', dirName, 'COMPLETE.txt'), 'a') as f:
					f.write(mpid + '\n')
				iSucc += 1
			except:
				print('FAILED FOR {}'.format(mpid))
				with open(os.path.join(td, 'autoparse', dirName, 'FAILED.txt'), 'a') as f:
					f.write(mpid + '\n')
					iFail += 1
	print('Parsed {} phases ({} successful, {} failed) in {} minutes'.format(iSucc+iFail, iSucc, iFail, (time.time()-t0)/60))

def parseGapMomentum(A):
	''' '''
	import pandas as pd
	d = pd.read_csv(os.path.join(td, 'autoparse/Gap-and-Type/GAP.txt'), index_col=0)
	dirName = 'Gap-Momentum'
	import time
	t0 = time.time()
	iSucc = 0
	iFail = 0
	print('Attempting to parse for a total of {} entries.'.format(len(A)))
	raw = []
	if 'COMPLETE.txt' in os.listdir(os.path.join(td, 'autoparse', dirName)):
		with open(os.path.join(td, 'autoparse', dirName, 'COMPLETE.txt'), 'r') as f:
			for line in f:
				raw.append(line.strip('\n'))
		print('Found list of {} completed phases.'.format(len(raw)))
	else:
		print('No list of completed phases detected.')
	for Ai in A:
		mpid = Ai['material_id']
		try:
			direct = d.ix[mpid].direct_gap
			if ~direct:
				form = Ai['pretty_formula']
				if mpid not in raw:
					try:
						bs = mpr.get_bandstructure_by_material_id(mpid)
						kcbm = bs.get_cbm()['kpoint'].cart_coords
						kvbm = bs.get_vbm()['kpoint'].cart_coords
						dk = np.linalg.norm(kcbm-kvbm)

						with open(os.path.join(td, 'autoparse', dirName, 'MOMENTUM.txt'), 'a') as f:
							f.write('{},{},{}\n'.format(mpid, form, dk))

						print('{:>20} {:>20} {:>20}'.format(mpid, form, dk))
						with open(os.path.join(td, 'autoparse', dirName, 'COMPLETE.txt'), 'a') as f:
							f.write(mpid + '\n')
						iSucc += 1
					except:
						print('FAILED FOR {}'.format(mpid))
						with open(os.path.join(td, 'autoparse', dirName, 'FAILED.txt'), 'a') as f:
							f.write(mpid + '\n')
							iFail += 1
		except:
			print('FAILED AT GETTING NATURE OF GAP FOR {}'.format(mpid))
	print('Parsed {} phases ({} successful, {} failed) in {} minutes'.format(iSucc+iFail, iSucc, iFail, (time.time()-t0)/60))
