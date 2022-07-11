"""Convienence functions for building dilute structures"""
from pymatgen.core import Structure
from dfttk import PRLStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from itertools import permutations, product, chain
import dfttk.structure_builders.substitutions as substitutions
import dfttk.utils
from collections import Counter
def dilute_substitution(structures, element_dict):
	"""
	Return dilute structures base on endmembers and elements

	Parameters
	----------
	structures: list of endmember structures from endmember.py

	elements: list of elements for substitutions

	Returns
	-------
	dilute structures: list
	"""
    dilute_strs=list()
    for strs in structures:
        sga = SpacegroupAnalyzer(strs)
        wyckoff_letters = sga.get_symmetry_dataset()['wyckoffs']
        equal_atoms=sga.get_symmetry_dataset()['equivalent_atoms']
        num=Counter(equal_atoms).keys()
        for i in num:
            dictstr=strs.as_dict()
            old_ele=dictstr['sites'][i]['species'][0]['element']
            for j in element_dict:
                if old_ele != j:
                    dictstr['sites'][i]['species'][0]['element']=j
                    dilute_strs.append(Structure.from_dict(dictstr))
    return dilute_strs
