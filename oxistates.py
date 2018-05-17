# Oxistates.py
# written by Mitchell Koerner
# 5/16/18

from pymatgen import MPRester
import pymatgen as pmg
from pymatgen.analysis.bond_valence import BVAnalyzer

mpr = MPRester("YOUR_API_KEY")
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