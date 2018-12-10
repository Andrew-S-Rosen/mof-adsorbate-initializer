from ase.build import molecule
from ase.geometry import get_distances
from ase import Atoms, Atom
from mai.tools import string_to_formula
import numpy as np
from ase.io import write

def construct_mof(ads_pos_optimizer,mof,ads_pos,site_idx):
	full_ads_species = ads_pos_optimizer.full_ads_species
	ads_species = ads_pos_optimizer.ads_species
	r_cut = ads_pos_optimizer.r_cut
	overlap_tol = ads_pos_optimizer.overlap_tol
	d_bond = ads_pos_optimizer.d_bond
	d_bond2 = ads_pos_optimizer.d_bond2
	angle = ads_pos_optimizer.angle
	angle2 = ads_pos_optimizer.angle2
	eta = ads_pos_optimizer.eta
	connect = ads_pos_optimizer.connect

	n_new_atoms = len(string_to_formula(ads_species))
	if '_grid' in full_ads_species:
		if ads_species == 'CH4':
			new_mof = add_CH4_SS(mof,site_idx,ads_pos)
		else:
			raise ValueError('Unsupported species for grid method')
	else:
		if n_new_atoms == 1:
			new_mof = add_monoatomic(mof,ads_species,ads_pos)
		elif n_new_atoms == 2:
			new_mof = add_diatomic(mof,ads_species,ads_pos,site_idx,
				d_bond=d_bond,angle=angle,eta=eta,connect=connect,
				r_cut=r_cut,overlap_tol=overlap_tol)
		elif n_new_atoms == 3:
			new_mof = add_triatomic(mof,ads_species,ads_pos,site_idx,
				d_bond1=d_bond,d_bond2=d_bond2,angle1=angle,angle2=angle2,
				connect=connect,r_cut=r_cut,overlap_tol=overlap_tol)	
		else:
			raise ValueError('Too many atoms in adsorbate')

	return new_mof

def add_monoatomic(mof,ads_species,ads_pos):
	"""
	Add adsorbate to the ASE atoms object

	Args:
		mof (ASE Atoms object): starting ASE Atoms object of structure

		ads_species (string): adsorbate species

		ads_pos (numpy array): 1D numpy array for the proposed
		adsorption position
	
	Returns:
		mof (ASE Atoms object): ASE Atoms object with adsorbate
	"""

	try:
		adsorbate = Atoms([Atom(ads_species,ads_pos)])
	except:
		raise ValueError('Unsupported adsorbate: '+ads_species)
	mof.extend(adsorbate)
	mof.wrap()
	
	return mof

def add_diatomic(mof,ads_species,ads_pos,site_idx,d_bond=1.25,angle=None,eta=1,
	connect=1,r_cut=2.5,overlap_tol=0.75):
	"""
	Add diatomic to the structure

	Args:
		mof (ASE Atoms object): starting ASE Atoms object of structure

		ads_species (string): adsorbate species

		ads_pos (array): 1D numpy array for the best adsorbate position
		
		site_idx (int): ASE index of site

		d_bond (float): X1-X2 bond length (defaults to 1.25)

		angle (float): site-X1-X2 angle (defaults to 180 degrees except for
		side-on in which it defaults to 90 or end-on O2 in which it defaults
		to 120)

		eta (int): denticity of end-on (1) or side-on (2) (defaults to 1)

		r_cut (float): cutoff distance for calculating nearby atoms when
		ranking adsorption sites (defualts to 2.5)

		overlap_tol (float): distance below which atoms are assumed to be
		overlapping (defualts to 0.75)

	Returns:
		mof (ASE Atoms object): ASE Atoms object with adsorbate
	"""
	
	#Construct X2 diatomic from string
	mol = Atoms([Atom('X',ads_pos),Atom('X',ads_pos)])
	ads_formula = string_to_formula(ads_species)
	if connect == 2:
		ads_formula[0], ads_formula[1] = ads_formula[1], ads_formula[0]
	X1 = ads_formula[0]
	X2 = ads_formula[1]
	try:
		mol[0].symbol = X1
		mol[1].symbol = X2
	except:
		raise ValueError('Invalid chemical species: '+ads_species)

	#Set default bond angle
	if angle is None:
		if eta == 1 and ads_species in ['O2','OO']:
			angle = 120.0
		elif eta == 2:
			angle = 90.0
		else:
			angle = 180.0
	while angle > 180:
		angle -= 180

	#Set X1 position and placeholder X2 position
	r_vec = mof.get_distance(2,site_idx,vector=True,mic=True)
	r_bond = d_bond*(r_vec/np.linalg.norm(r_vec))
	mol[0].position = ads_pos
	mol[1].position = ads_pos+r_bond

	#Add diatomic to the structure
	if site_idx is None:
		raise ValueError('Site index must not be None')
	mof.extend(mol)
	mof.wrap()
	mof.set_angle(site_idx,-2,-1,angle)
	mof.wrap()

	#Make adsorption mode side-on if requested
	if eta == 2:
		line = mof.get_distance(-1,-2,vector=True,mic=True)
		shift = (line/np.linalg.norm(line))*d_bond/2
		mof[-2].position += shift
		mof[-1].position += shift
		mof.wrap()
	elif eta > 2:
		raise ValueError('Wrong value for eta: '+str(eta))

	#If not linear, sweep possible angles for X2 to minimize sterics
	if angle != 180:
		dtheta = 10
		n_angles = int(360/dtheta)
		for i in range(n_angles):
			mof_temp = mof.copy()
			
			if eta == 1:
				temp_atoms = 2
				mof_temp.extend(Atoms([Atom('X',ads_pos+r_bond+1e-6)]))
				mof_temp.wrap()
				mof_temp.set_angle(site_idx,-3,-1,360-angle)
				write('mof.cif',mof_temp)
			elif eta == 2:
				temp_atoms = 1

			r_temp = mof_temp.get_distance(-1,-2,vector=True,mic=True)
			pos_temp = mof_temp[-1].position+r_temp/2
			mof_temp.extend(Atoms([Atom('X',pos_temp)]))

			ads_temp = mof_temp.copy()[-(1+temp_atoms):]
			ads_temp.rotate(i*dtheta,mof_temp.get_distance(site_idx,-1,
				vector=True,mic=True),pos_temp)
			del mof_temp[-temp_atoms:]
			mof_temp[-1].position = ads_temp[-(1+temp_atoms)].position

			dist_mat = mof_temp.get_distances(-1,np.arange(0,
				len(mof_temp)-1).tolist(),vector=False,mic=True)
			NNs = sum(dist_mat <= r_cut)

			if i == 0:
				ads_pos2 = mof_temp[-1].position
				old_min_NNs = NNs
			elif sum(dist_mat <= overlap_tol) == 0 and NNs < old_min_NNs:
				ads_pos2 = mof_temp[-1].position
				old_min_NNs = NNs
		mof[-1].position = ads_pos2
	mof.wrap()

	return mof

def add_triatomic(mof,ads_species,ads_pos,site_idx,d_bond1=1.25,d_bond2=None,
	angle1=None,angle2=None,connect=1,r_cut=2.5,overlap_tol=0.75):
	"""
	Add triatomic to the structure
	"""
	ads_formula = string_to_formula(ads_species)
	X3 = ads_formula[2]

	if d_bond2 is None:
		d_bond2 = d_bond1 
	if angle1 is None:
		if ads_species in ['H2O','HOO','OHH']:
			angle1 = 104.5
		else:
			angle1 = 180.0
	if angle2 is None:
		if connect == 1 or connect == 3:
			angle2 = 180.0
		else:
			angle2 = angle1
	if connect == 2:
		ads_formula[0], ads_formula[1] = ads_formula[1], ads_formula[0]
	elif connect == 3:
		ads_formula[0], ads_formula[2] = ads_formula[2], ads_formula[0]

	di_ads_species = ''.join(ads_formula[0:2])
	mof = add_diatomic(mof,di_ads_species,ads_pos,site_idx,d_bond=d_bond1,
		angle=angle1,r_cut=r_cut,overlap_tol=overlap_tol)

	if connect == 1 or connect == 3:
		r_vec = mof.get_distance(-2,-1,vector=True,mic=True)
		pos_temp = mof[-1].position+d_bond2*(r_vec/np.linalg.norm(r_vec))
		mof.extend(Atoms([Atom(X3,pos_temp)]))
		mof.wrap()
		mof.set_angle(-3,-2,-1,angle2)
	elif connect == 2:
		if angle1 == 180 and angle2 == 180:
			raise ValueError('It is not possible to have a linear triatomic '+
				'with connect=2')
		r_vec = mof.get_distance(site_idx,-2,vector=True,mic=True)
		r_bond = d_bond1*(r_vec/np.linalg.norm(r_vec))
		mof.extend(Atoms([Atom(X3,ads_pos+r_bond+1e-6)]))
		mof.wrap()
		mof.set_angle(-2,-3,-1,angle2)
	else:
		raise ValueError('Connecting atom must have value of <= 3')
	mof.wrap()

	return mof

def add_CH4_SS(mof,site_idx,ads_pos):
	"""
	Add CH4 to the structure
 ASE Atoms object of structure

		site_idx (int): ASE index of site based on single-site model

	Args:
		mof (ASE Atoms object): starting
		ads_pos (array): 1D numpy array for the best adsorbate position
			
	Returns:
		mof (ASE Atoms object): ASE Atoms object with adsorbate
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
	CH4.set_distance(0,2,CH_length,fix=0,mic=True)
	CH4.set_angle(1,0,2,CH_angle)
	CH4.set_distance(0,3,CH_length,fix=0,mic=True)
	CH4.set_angle(1,0,3,CH_angle)
	CH4.set_dihedral(2,1,0,3,-CH_dihedral)
	CH4.set_distance(0,4,CH_length,fix=0,mic=True)
	CH4.set_angle(1,0,4,CH_angle)
	CH4.set_dihedral(2,1,0,4,CH2_dihedral)

	#Add CH4 molecule to the structure
	mof.extend(CH4)
	mof.wrap()
	
	return mof
