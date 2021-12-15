"""Convienence functions for building endmembers"""
from pymatgen.core import Structure
from dfttk import PRLStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from itertools import permutations, product, chain
import dfttk.structure_builders.substitutions as substitutions
import dfttk.utils

def get_sublattice_information(structure, equivalent_wyckoff_sites=None):
	"""
	Return defined sublattice_name based on wyckoff positions and other sublattice information

	Parameters
	----------
	structure : pymatgen.Structure

	equivalent_wyckoff_sites: dict
		Dict of Wyckoff sites that are treated as the same sublattice, e.g. {'b': 'f'} will
		give combine Wyckoff site 'b' and Wyckoff site 'f' into one sublattice named as 'b'. For
		multiple equivalent sites, one can use {'b': 'f', 'c': 'i'}

	Returns
	-------

	site_list: list

	subl_model_name: list
		Name for each sublattice based on the name of wyckoff position
	"""
	structure.replace_species({sp.name: "H" for sp in structure.species})
	sga = SpacegroupAnalyzer(structure)
	wyckoff_sites = sga.get_symmetry_dataset()['wyckoffs']
	if equivalent_wyckoff_sites is not None:
		true_sublattice_sites = [equivalent_wyckoff_sites[i] if i in equivalent_wyckoff_sites else i for i in wyckoff_sites]
	else:
		true_sublattice_sites = wyckoff_sites
	site_list=[]
	for i in true_sublattice_sites:
		site_dict={'sublattice_sites': i}
		site_list.append(site_dict)
		subl_model_name = sorted(set(true_sublattice_sites))
	return site_list, subl_model_name

def get_templates(structure, equivalent_wyckoff_sites=None):
	"""
	Return templates of structure and configuration for substitution of endmembers

	Parameters
	----------
	structure : pymatgen.Structure

	equivalent_wyckoff_sites: dict
		Dict of Wyckoff sites that are treated as the same sublattice, e.g. {'b': 'f'} will
		give combine Wyckoff site 'b' and Wyckoff site 'f' into one sublattice named as 'f'. For
		multiple equivalent sites, one can use {'b': 'f', 'c': 'i'}
	
	Returns
	-------
	template_structure: pymatgen.Structure

	template_configuration: list
		Template configuration from the template structure
	"""
	if equivalent_wyckoff_sites is not None:
		site_list, true_sublattices=get_sublattice_information(structure, equivalent_wyckoff_sites=equivalent_wyckoff_sites)
	else:
		site_list, true_sublattices=get_sublattice_information(structure)
	dict_struct=structure.as_dict()
	for i in range(0,len(dict_struct['sites'])):
		dict_struct['sites'][i]['properties'].update(site_list[i])
		dict_struct['sites'][i]['species'][0]['element']=site_list[i]['sublattice_sites']
	template_configuration=[]
	for i in range(0,len(true_sublattices)):
		for element in dict_struct['sites']:
			if element['species'][0]['element']==true_sublattices[i]:
				element['species'][0]['element'] = str(Element.from_Z(i+1))
		template_configuration.append(str(Element.from_Z(i+1)))
	template_structure=Structure.from_dict(dict_struct)
	return template_structure, template_configuration


def get_endmembers_with_templates(template_structure, template_configuration, sublattice_configuration):
	"""
	Return endmembers of the structure

	Parameters
	----------
	template_structure: pymatgen.Structure

	template_configuration: list
		Configuration in the template structure, should be consitent with sublattice model, e.g., one unique
		element for each sublattice. If input manually, make sure be consistent with the sublattice_model_name 

	sublattice_configuration: list
		List of configurations for each sublattice

	Returns
	-------
	endmembers: list 
		List of structures in pymatgen.Structure
	"""
	element_dict=dict(map(lambda x, y: [x, y], template_configuration, sublattice_configuration))
	subl_combination=[]
	for comb in product(*map(element_dict.get, sorted(element_dict))):
		subl_combination.append(list(comb))
	endmembers=[]
	for i in subl_combination:
		endmembers.append(substitutions.substitute_configuration(template_structure, [template_configuration], [i]))
	return endmembers
