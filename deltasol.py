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

    Bugs:
        error handling is broken :( I'm so bad at python error handling... (possibly fixed)
    """
    species = [str(i) for i in structure.species]
    if len(set(species).intersection(bad_valences)) > 0:     #somehow broken (maybe not anymore)
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

def product(array):
    """
    Returns the product of all values in the array
    """
    result = 1
    for i in array:
        result *= i
    return(result)

def multiply_cell(structure, valence, accept = (60, 80), verbose = False):
    """
    Multiplies a cell until it has the proper number of valence electrons
     for 1/60 - 1/80 valence removed

    Inputs:
        structure:          pymatgen.structure to be multiplied
        valence:            valence number of structure (from count_valence)
        accept:             (min, max) for 1/min to 1/max
        verbose:            verbose output for longer runs with >1 electron removed

    Outputs:
        final_structure:    new multiplied structure to be run
        ncells:             atoms in final structure / atoms in initial
        nelectrons:         correct number of electrons to add or remove for calculation
        mult_cell:          multiplication of the structure that was used [2,3,1] means
                             x*2 y*3 z*1

    Bugs:
        Prefers higher values of accept (towards 1/80 valence electrons) if the cell could
         somehow choose either 80 or 70, it would pick 80 even though 70 is closer to middle
    """
    assert type(valence) == int, "valence must be an integer"
    ncells = 1
    mult_cell = [1,1,1]
    final_structure = structure
    while True:
        nelectrons = 1
        test = valence * ncells / nelectrons
        while test >= min(accept):
            if (test >= min(accept)) and (test <= max(accept)):
                return(final_structure, ncells, nelectrons, test, mult_cell)
            if verbose:
                print("{} {} {} {}".format(ncells, nelectrons, test, mult_cell))
            nelectrons += 1
            test = valence * ncells / nelectrons
        abc = final_structure.lattice.abc
        index = abc.index(min(abc))
        mult_cell[index] += 1
        ncells = product(mult_cell)
        final_structure = structure * mult_cell