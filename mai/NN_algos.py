from pymatgen.analysis.local_env import (MinimumVIRENN, VoronoiNN,
	JmolNN, MinimumDistanceNN, MinimumOKeeffeNN, BrunnerNN_real,
	BrunnerNN_relative, BrunnerNN_reciprocal, EconNN, CutOffDictNN,
	Critic2NN, OpenBabelNN, CovalentBondNN, CrystalNN)
from pymatgen.io import ase as pm_ase
import warnings

def get_NNs_pm(atoms,site_idx,NN_method):
	"""
	Get coordinating atoms to the adsorption site

	Args:
		atoms (Atoms object): atoms object of MOF

		site_idx (int): ASE index of adsorption site
		
		NN_method (string): string representing the desired Pymatgen
		nearest neighbor algorithm:
		refer to http://pymatgen.org/_modules/pymatgen/analysis/local_env.html
	
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
	elif NN_method == 'dict':
		#requires a cutoff dictionary located in the pwd
		nn_object = CutOffDictNN(cut_off_dict='cut_off_dict.txt')
	elif NN_method == 'critic2':
		nn_object = Critic2NN()
	elif NN_method == 'openbabel':
		nn_object = OpenBabelNN()
	elif NN_method == 'covalent':
		nn_object = CovalentBondNN()
	elif NN_method == 'crystal':
		nn_object = CrystalNN(porous_adjustment=True)
	elif NN_method == 'crystal_nonporous':
		nn_object = CrystalNN(porous_adjustment=False)
	else:
		raise ValueError('Invalid NN algorithm specified')

	#Find coordinating atoms
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		neighbors = nn_object.get_nn_info(struct,site_idx)
	return [neighbor['site_index'] for neighbor in neighbors]