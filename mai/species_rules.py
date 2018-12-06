from ase.build import molecule
from ase.geometry import get_distances
from ase import Atoms, Atom
import numpy as np

def add_monoatomic(mof,ads_species,ads_pos):
	"""
	Add adsorbate to the ASE atoms object

	Args:
		ads_pos (numpy array): 1D numpy array for the proposed
		adsorption position
	
	Returns:
		new_mof (ASE Atoms object): Atoms object for new
		structure with adsorbate
	"""

	try:
		adsorbate = Atoms([Atom(ads_species,ads_pos)])
	except:
		raise ValueError('Unsupported adsorbate: '+ads_species)
	mof.extend(adsorbate)
	mof.wrap()
	
	return mof

def add_diatomic(mof,ads_species,ads_pos,site_idx,d_bond=1.25,angle=180.0,eta=1):
	"""
	Add diatomic to the structure

	Args:
		site_idx (int): ASE index of site
		
		ads_pos (array): 1D numpy array for the best adsorbate position
		
		atoms (ASE Atoms object): Atoms object of structure
	
		X1 (string): chemical symbol for atom 1 of diatomic

		X2 (string): chemical symbol for atom 2 of diatomic (defualts to X1)
	
		d_bond (float): X1-X2 bond length (defaults to 2 A)

		angle (float): site-X1-X2 angle (defaults to 180 degrees)

		center (bool): if center of diatomic should be put at ads_site,
		as opposed to X1 (defaults to False)

	Returns:
		atoms (ASE Atoms object): new ASE Atoms object with adsorbate
	
		n_new_atoms (int): number of atoms in adsorbate
	"""
	
	mol = molecule('H2')
	X1 = ads_species[0]
	if ads_species[1] == '2':
		X2 = X1
	else:
		X2 = ads_species[1]
	try:
		mol[0].symbol = X1
		mol[1].symbol = X2
	except:
		raise ValueError('Invalid chemical species: '+ads_species)

	D,D_len = get_distances([ads_pos],mof[site_idx].position,cell=mof.cell,pbc=mof.pbc)
	r_vec = D[0,0]
	r_bond = (r_vec/np.linalg.norm(r_vec))*d_bond

	if eta == 1:
		mol[0].position = ads_pos
		mol[1].position = ads_pos+r_bond
	elif eta == 2:#fixme
		mol[0].position = ads_pos
		mol[1].position = ads_pos+r_bond
	else:
		raise ValueError('Invalid value for eta: '+str(eta))

	#Add diatomic to the structure
	if site_idx is None:
		raise ValueError('Site index must not be None')
	mof.extend(mol)
	mof.set_angle(site_idx,-2,-1,angle=angle)
	mof.wrap()

	return mof

def add_CH4_SS(mof,site_idx,ads_pos):
	"""
	Add CH4 to the structure

	Args:
		site_idx (int): ASE index of site based on single-site model

		ads_pos (array): 1D numpy array for the best adsorbate position
		
		atoms (ASE Atoms object): Atoms object of structure
	
	Returns:
		atoms (ASE Atoms object): new ASE Atoms object with adsorbate
	"""
	
	#Get CH4 parameters
	CH4 = molecule('CH4')
	CH_length = CH4.get_distance(0,1)
	CH_angle = CH4.get_angle(1,0,2)
	CH2_dihedral = CH4.get_dihedral(2,1,0,4)
	CH_length = CH4.get_distance(0,1)
	CH_angle = CH4.get_angle(1,0,2)
	CH_dihedral = CH4.get_dihedral(2,1,0,4)

	#Add CH4 to ideal adsorption position
	CH4[0].position = ads_pos

	#Make one of the H atoms colinear with adsorption site and C
	D,D_len = get_distances([ads_pos],mof[site_idx].position,cell=mof.cell,pbc=mof.pbc)
	r_vec = D[0,0]
	r = (r_vec/np.linalg.norm(r_vec))*CH_length

	#Construct rest of CH4 using Z-matrix format
	CH4[1].position = ads_pos+r
	CH4.set_distance(0,2,CH_length,fix=0)
	CH4.set_angle(1,0,2,CH_angle)
	CH4.set_distance(0,3,CH_length,fix=0)
	CH4.set_angle(1,0,3,CH_angle)
	CH4.set_dihedral(2,1,0,3,-CH_dihedral)
	CH4.set_distance(0,4,CH_length,fix=0)
	CH4.set_angle(1,0,4,CH_angle)
	CH4.set_dihedral(2,1,0,4,CH2_dihedral)

	#Add CH4 molecule to the structure
	mof.extend(CH4)
	mof.wrap()
	
	return mof