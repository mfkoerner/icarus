# deltasol.py
# written by Mitchell Koerner
# 5/24/18

import numpy as np
import pymatgen as pmg


# total number of valence electrons for each atom (counting octet for main group, outermost s and d for transition metals)
atomic_valence={'H':  1,'He': 0,'Li': 1,'Be': 2,'B':  3,'C':  4,'N':  5,'O':  6,'F':  7,'Ne': 0,'Na': 1,'Mg': 2,'Al': 3,'Si': 4,'P':  5,'S':  6,
'Cl': 7,'Ar': 0,'K':  1,'Ca': 2,'Sc': 3,'Ti': 4,'V':  5,'Cr': 6,'Mn': 7,'Fe': 8,'Co': 9,'Ni': 10,'Cu': 11,'Zn': 12,'Ga': 3,'Ge': 4,'As': 5,
'Se': 6,'Br': 7,'Kr': 0,'Rb': 1,'Sr': 2,'Y':  3,'Zr': 4,'Nb': 5,'Mo': 6,'Tc': 7,'Ru': 8,'Rh': 9,'Pd': 10,'Ag': 11,'Cd': 12,'In': 3,'Sn': 4,
'Sb': 5,'Te': 6,'I':  7,'Xe': 0,'Cs': 1,'Ba': 2,'La': 3,'Ce': None,'Pr': None,'Nd': None,'Pm': None,'Sm': None,'Eu': None,'Gd': None,'Tb': None,
'Dy': None,'Ho': None,'Er': None,'Tm': None,'Yb': None,'Lu': None,'Hf': 4,'Ta': 5,'W':  6,'Re': 7,'Os': 8,'Ir': 9,'Pt': 10,'Au': 11,'Hg': 12,
'Tl': 3,'Pb': 4,'Bi': 5,'Po': 6,'At': 7,'Rn': 0,'Fr': None,'Ra': None,'Ac': None,'Th': None,'Pa': None,'U':  None,'Np': None,'Pu': None,'Am': None,
'Cm': None}
bad_valences = {i for i in atomic_valence.keys() if atomic_valence[i] is None}


def count_valence(structure, return_all = False):
    """
    Gives the total nubmer of valence electrons in a structure according to rules for atomic_valence dictionary

    Inputs:
        structure:  pymatgen structure object
        return_all: when true, returns (total number of valence electrons, species list, atomic valences list)

    Outputs:
        total:      integer of total number of valence electrons

        if return_all
        total:      same
        species:    list of pymatgen specie objects (it's just structure.species)
        valences:   list of valences from atomic_valence for each atom in species
    """
    species = structure.species
    if len(set(species).intersection(bad_valences)) > 0:
        print (set(species).intersection(bad_valences) != {})
        raise ValueError("compound contains F-Block elements")
    valences = [atomic_valence[i.name] for i in species]
    total = sum(valences)
    if return_all:
        return(total, species, valences)
    else:
        return(total)

def find_NELECT(fp):
    """
    Finds NELECT from an OUTCAR file

    Inputs:
        fp:             filepath to OUTCAR file which should be read

    Outputs:
        n_electrons:    Number of electrons in the run OUTCAR is representing
                        comes from NELECT = ???.???? line
    """
    with open(fp, "r") as f:
        inlines = f.readlines()
    NELECT_lines = [i for i in inlines if "NELECT" in i]
    assert len(NELECT_lines) == 1, "There should be exactly one line in OUTCAR with 'NELECT'"
    line = NELECT_lines[0]
    n_electrons = float(line.split("=")[1].split()[0].strip())
    return(n_electrons)