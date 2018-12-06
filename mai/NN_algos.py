from pymatgen.analysis.local_env import (MinimumVIRENN, VoronoiNN,
	JmolNN, MinimumDistanceNN, MinimumOKeeffeNN, BrunnerNN_real,
	BrunnerNN_relative, BrunnerNN_reciprocal, EconNN)
from pymatgen.io import ase as pm_ase

def get_NNs_pm(atoms,site_idx,NN_method):
	"""
	Get coordinating atoms to the adsorption site

	Args:
		atoms (Atoms object): atoms object of MOF

		site_idx (int): ASE index of adsorption site
		
		NN_method (string): string representing the desired Pymatgen
		nearest neighbor algorithm (accepts 'vire','voronoi','jmol',
		'min_dist','okeeffe','brunner_real','brunner_recpirocal',
		'brunner_relative', and 'econ')
	
	Returns:
		neighbors_idx (list of ints): ASE indices of coordinating atoms
	"""
	#Convert ASE Atoms object to Pymatgen Structure object
	bridge = pm_ase.AseAtomsAdaptor()
	struct = bridge.get_structure(atoms)

	#Select Pymatgen NN algorithm
	NN_method = NN_method.lower()
	if NN_method == 'vire':
		nn_object = MinimumVIRENN()
	elif NN_method == 'voronoi':
		nn_object = VoronoiNN()
	elif NN_method == 'jmol':
		nn_object = JmolNN()
	elif NN_method == 'min_dist':
		nn_object = MinimumDistanceNN()
	elif NN_method == 'okeeffe':
		nn_object = MinimumOKeeffeNN()
	elif NN_method == 'brunner_real':
		nn_object = BrunnerNN_real()
	elif NN_method == 'brunner_recpirocal':
		nn_object = BrunnerNN_reciprocal()
	elif NN_method == 'brunner_relative':
		nn_object = BrunnerNN_relative()
	elif NN_method == 'econ':
		nn_object = EconNN()
	else:
		raise ValueError('Invalid NN algorithm specified')

	#Find coordinating atoms
	neighbors = nn_object.get_nn_info(struct,site_idx)
	neighbors_idx = []
	for neighbor in neighbors:
		neighbors_idx.append(neighbor['site_index'])

	return neighbors_idx