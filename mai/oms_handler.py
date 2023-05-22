import numpy as np
import os
import json
from ase.io import read

def get_ase_NN_idx(atoms,coords):
	"""
	Get the ASE indices for the coordinating atoms

	Args:
		atoms (Atoms object): ASE Atoms object for the MOF
		
		coords (numpy array): coordinates of the coordinating atoms
	
	Returns:
		ase_NN_idx (list of ints): ASE indices of the coordinating atoms
	"""
	ase_NN_idx = []
	tol = 0.01
	for i in range(np.shape(coords)[0]):
		nn_fail = True
		for j, element in enumerate(atoms):
			if (sum(coords[i,:] >= element.position-tol) == 3 and 
			sum(coords[i,:] <= element.position+tol) == 3):
				ase_NN_idx.append(j)
				nn_fail = False
				break
		if nn_fail:
			raise ValueError('Could not match to ASE index')

	return ase_NN_idx

def get_ase_oms_idx(atoms,coords):
	"""
	Get the ASE index of the OMS

	Args:
		atoms (Atoms object): ASE Atoms object for the MOF
		
		coords (numpy array): coordinates of the OMS
	
	Returns:
		ase_oms_idx (int): ASE index of OMS
	"""
	tol = 0.01
	oms_idx_fail = True
	for i, element in enumerate(atoms):
		if (sum(coords >= element.position-tol) == 3 and 
			sum(coords <= element.position+tol) == 3):
			ase_oms_idx = i
			oms_idx_fail = False
			break
	if oms_idx_fail:
		raise ValueError('Could not match to ASE index')
	return ase_oms_idx

def get_omd_data(oms_data_path,name,atoms):
	"""
	Get info about the open metal site from OpenMetalDetector
	results files

	Args:
		oms_data_path (string): path to the OpenMetalDetector results

		name (string): name of the MOF
		
		atoms (ASE Atoms object): Atoms object for the MOF
	
	Returns:
		omsex_dict (dict): dictionary of data from the OpenMetalDetector results
	"""
	cnums_all = []
	oms_coords_all = []
	oms_idx_all = []
	oms_sym_all = []
	NN_coords_all = []
	NN_idx_all = []
	oms_sym_all = []

	json_path = os.path.join(oms_data_path, name, f'{name}.json')
	if os.stat(json_path).st_size == 0:
		return None
	with open(json_path) as f:
		oms_results = json.load(f)
	if not oms_results['cif_okay'] or oms_results['problematic']:
		print(f'{name}: Open Metal Detector failed')
		return None
	if not oms_results['has_oms']:
		return None

	for i, site in enumerate(oms_results['metal_sites']):
		if site['problematic'] or not site['is_open'] or not site['unique']:
			continue
		cnum = site['number_of_linkers']
		sphere = read(
			os.path.join(
				oms_data_path, name, f'first_coordination_sphere{str(i)}.cif'
			)
		)
		oms_coords = sphere[0].position
		oms_sym = sphere[0].symbol
		oms_idx = get_ase_oms_idx(atoms,oms_coords)
		NN_coords = sphere[1:].positions
		NN_idx = get_ase_NN_idx(atoms,NN_coords)

		cnums_all.append(cnum)
		oms_coords_all.append(oms_coords)
		oms_idx_all.append(oms_idx)
		oms_sym_all.append(oms_sym)
		NN_coords_all.append(NN_coords)
		NN_idx_all.append(NN_idx)
	return {
		'cnums': cnums_all,
		'oms_coords': oms_coords_all,
		'oms_idx': oms_idx_all,
		'oms_sym': oms_sym_all,
		'NN_coords': NN_coords_all,
		'NN_idx': NN_idx_all,
	}