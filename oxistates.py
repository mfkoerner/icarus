# Oxistates.py
# written by Mitchell Koerner
# 5/16/18

import config
import matgen as mg
from pymatgen import MPRester
import pymatgen as pmg
from pymatgen.analysis.bond_valence import BVAnalyzer

mpr = MPRester(config.matprojapi)
bva = BVAnalyzer()

def oxidation(mpid):
	"""
	returns oxidation states of atoms in a structure specified by materials project id
	inputs: mpid string, example: "mp-546794" (SiO)
	outputs:
		oxistates:	ordered list of oxidation states as sites are labled in mystruct
		mystruct:	structure object which contains sites that correspond to oxistates
	use mystruct.sites to create a matching site to oxidation state
	also, mystruct.sites[i].specie (no s) can be used to get species of each site
	"""
	mystruct = mpr.get_structure_by_material_id(mpid)
	try:
		oxistates = bva.get_valences(mystruct)
	except ValueError:
		oxistates = None
	return(oxistates, mystruct)

def printOxi(mpids, display=False):
	form = mg.getFormulae()
	out = {}
	for mpid in mpids:
		o, s = oxidation(mpid)
		if display:
			print('{:<20}{:<20}'.format(mpid, mg.getFormula(mpid)))
		if o:
			pairs = set([(oi, si.specie) for oi, si in zip(o, s.sites)])
			for p in pairs:
				if display:
					print(' ... {:>2} {}'.format(p[0], p[1].name))
		else:
			pairs = set()
			if display:
				print(' ... UNABLE TO COMPUTE OXIDATION STATES')
		out[mpid] = pairs
	return out
