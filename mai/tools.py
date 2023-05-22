import re

def get_refcode(atoms_filename):
	"""
	Get the name of the MOF

	Args:
		atoms_filename (string): filename of the ASE Atoms object (accepts
		CIFS, POSCARs, and CONTCARs)
		
	Return:
		refcode (string): name of MOF (defaults to 'mof' if the original
		filename is just named CONTCAR or POSCAR)
	"""
	if '.cif' in atoms_filename:
		refcode = atoms_filename.split('.cif')[0]
	elif 'CAR_' in atoms_filename:
		refcode = atoms_filename.split('CAR_')[-1]
	elif '_CAR' in atoms_filename:
		refcode = atoms_filename.split('_CAR')[0]
	elif atoms_filename in ['CONTCAR', 'POSCAR']:
		refcode = 'mof'
	else:
		raise ValueError('Unknown file naming scheme')

	return refcode

def string_to_formula(species_string):
	"""
	Convert a species string to a chemical formula

	Args:
		species_string (string): string of atomic/molecular species
		
	Return:
		formula (string): stoichiometric chemical formula
	"""
	syms = re.findall('[A-Z][a-z]?',species_string)
	formula = []
	for sym in syms:
		end_idx = re.search(sym,species_string).end()
		try:
			diff = len(species_string)-len(species_string[end_idx:])
			next_idx = (diff-1)+re.search('[A-Z][a-z]?',species_string[end_idx:]).end()
			end_char = species_string[end_idx:next_idx]
		except:
			end_char = species_string[end_idx:]	
		try:
			stoich = int(end_char)
		except:
			stoich = 1
		formula.extend([sym]*stoich)
	return formula