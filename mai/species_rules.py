from ase.build import molecule
from ase.geometry import get_distances
from ase import Atoms, Atom
from mai.tools import string_to_formula
import numpy as np

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

	return mof

def add_diatomic(mof,ads_species,ads_pos,site_idx,d_X1X2=1.25,ang_MX1X2=None,eta=1,
	connect=1,r_cut=2.5,overlap_tol=0.75):
	"""
	Add diatomic to the structure

	Args:
		mof (ASE Atoms object): starting ASE Atoms object of structure

		ads_species (string): adsorbate species

		ads_pos (array): 1D numpy array for the best adsorbate position
		
		site_idx (int): ASE index of site

		d_X1X2 (float): X1-X2 bond length (defaults to 1.25)

		ang_MX1X2 (float): site-X1-X2 angle (defaults to 180 degrees except for
		side-on in which it defaults to 90)

		eta (int): denticity of end-on (1) or side-on (2) (defaults to 1)

		connect (int): the connecting atom in the species string (defaults to 1)

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

	#Set default bond ang_MX1X2
	if ang_MX1X2 is None:
		if eta == 1:
			ang_MX1X2 = 180.0
		elif eta == 2:
			ang_MX1X2 = 90.0
	while ang_MX1X2 > 180:
		ang_MX1X2 -= 180

	#Set X1 position and placeholder X2 position
	n_start = len(mof)
	mof.extend(mol[0])
	r_vec = mof.get_distance(site_idx,-1,vector=True,mic=True)
	r_bond = d_X1X2*(r_vec/np.linalg.norm(r_vec))
	mol[1].position = ads_pos+r_bond
	mof.extend(mol[1])
	mof[-1].position += 1e-6
	mof.set_angle(site_idx,-2,-1,ang_MX1X2)

	#Make adsorption mode side-on if requested
	if eta == 2:
		line = mof.get_distance(-1,-2,vector=True,mic=True)
		shift = (line/np.linalg.norm(line))*d_X1X2/2
		mof[-2].position += shift
		mof[-1].position += shift
	elif eta > 2:
		raise ValueError('Wrong value for eta: '+str(eta))

	#If not linear, sweep possible angles for X2 to minimize sterics
	if ang_MX1X2 != 180:
		dtheta = 10
		n_angles = int(360/dtheta)

		if eta == 1:
			temp_atoms = 2
			mof.extend(Atoms([Atom('X',ads_pos+r_bond)]))
			mof[-1].position += 1e-6
			mof.set_angle(site_idx,-3,-1,360-ang_MX1X2)
		elif eta == 2:
			temp_atoms = 1
		r_temp = mof.get_distance(-1,-2,vector=True,mic=True)
		pos_temp = mof[-1].position+r_temp/2
		mof.extend(Atoms([Atom('X',pos_temp)]))

		for i in range(n_angles):
			mof_temp = mof.copy()
			ads_temp = mof.copy()[-3:]
			ads_temp.rotate(i*dtheta,mof_temp.get_distance(site_idx,-1,
				vector=True,mic=True),pos_temp)
			del mof_temp[-temp_atoms:]

			mof_temp[-1].position = ads_temp[0].position
			if eta == 2:
				mof_temp[-2].position = ads_temp[1].position
			
			dist_mat = mof_temp.get_distances(-1,np.arange(0,
				n_start-1).tolist(),vector=False,mic=True)
			n_overlap = sum(dist_mat <= overlap_tol)
			NNs = sum(dist_mat <= r_cut)
			min_dist = np.min(dist_mat)

			if eta == 2:
				dist_mat2 = mof_temp.get_distances(-2,np.arange(0,
					n_start-1).tolist(),vector=False,mic=True)
				n_overlap2 = sum(dist_mat2 <= overlap_tol)
				NNs2 = sum(dist_mat2 <= r_cut)
				min_dist2 = np.min(dist_mat2)				
				if min_dist2 < min_dist:
					min_dist = min_dist2
				if NNs2 > NNs:
					NNs = NNs2
				if n_overlap2 > n_overlap:
					n_overlap = n_overlap2

			if i == 0:
				best_mof = mof_temp.copy()
				old_min_NNs = NNs
				old_min_dist = min_dist
				old_overlap = n_overlap
			elif n_overlap <= old_overlap and (NNs < old_min_NNs or 
				(NNs == old_min_NNs and min_dist > old_min_dist)):
				best_mof = mof_temp.copy()
				old_min_NNs = NNs
				old_min_dist = min_dist
				old_overlap = n_overlap

		mof = best_mof

	return mof

def add_triatomic(mof,ads_species,ads_pos,site_idx,d_X1X2=1.25,d_X2X3=None,
	ang_MX1X2=None,ang_triads=None,connect=1,r_cut=2.5,overlap_tol=0.75):
	"""
	Add triatomic to the structure

	Args:
		mof (ASE Atoms object): starting ASE Atoms object of structure

		ads_species (string): adsorbate species

		ads_pos (array): 1D numpy array for the best adsorbate position
		
		site_idx (int): ASE index of site

		d_X1X2 (float): X1-X2 bond length (defaults to 1.25)

		d_X2X3 (float): X2-X3 bond length for connect == 1 or
		X1-X3 bond length for connect == 2 (defaults to d_X1X2)

		ang_MX1X2 (float): site-X1-X2 angle (defaults to 180 degrees)

		ang_triads (float): triatomic angle (defaults to 180 degrees for connect == 1
		and ang_MX1X2 for connect == 2)

		connect (int): the connecting atom in the species string (defaults to 1)

		r_cut (float): cutoff distance for calculating nearby atoms when
		ranking adsorption sites (defualts to 2.5)

		overlap_tol (float): distance below which atoms are assumed to be
		overlapping (defualts to 0.75)

	Returns:
		mof (ASE Atoms object): ASE Atoms object with adsorbate
	"""
	ads_formula = string_to_formula(ads_species)
	X3 = ads_formula[2]

	if d_X2X3 is None:
		d_X2X3 = d_X1X2 
	if ang_MX1X2 is None:
		ang_MX1X2 = 180.0
	if ang_triads is None:
		if connect == 1 or connect == 3:
			ang_triads = 180.0
		else:
			ang_triads = ang_MX1X2
	if connect == 2:
		ads_formula[0], ads_formula[1] = ads_formula[1], ads_formula[0]
	elif connect == 3:
		ads_formula[0], ads_formula[2] = ads_formula[2], ads_formula[0]

	di_ads_species = ''.join(ads_formula[0:2])
	mof = add_diatomic(mof,di_ads_species,ads_pos,site_idx,d_X1X2=d_X1X2,
		ang_MX1X2=ang_MX1X2,r_cut=r_cut,overlap_tol=overlap_tol)

	if connect == 1 or connect == 3:
		r_vec = mof.get_distance(-2,-1,vector=True,mic=True)
		pos_temp = mof[-1].position+d_X2X3*(r_vec/np.linalg.norm(r_vec))
		mof.extend(Atoms([Atom(X3,pos_temp)]))
		mof[-1].position += 1e-6
		mof.set_angle(-3,-2,-1,ang_triads)
	elif connect == 2:
		if ang_MX1X2 == 180 and ang_triads == 180:
			raise ValueError('It is not possible to have a linear triatomic '+
				'with connect=2')
		r_vec = mof.get_distance(site_idx,-2,vector=True,mic=True)
		r_bond = d_X2X3*(r_vec/np.linalg.norm(r_vec))
		mof.extend(Atoms([Atom(X3,ads_pos+r_bond)]))
		mof[-1].position += 1e-6
		mof.set_angle(-2,-3,-1,ang_triads)
	else:
		raise ValueError('Connecting atom must have value of <= 3')

	return mof

def add_CH4_SS(mof,site_idx,ads_pos):
	"""
	Add CH4 to the structure from single-site model

	Args:
		mof (ASE Atoms object): starting ASE Atoms object of structure

		site_idx (int): ASE index of site based on single-site model

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
	
	return mof
