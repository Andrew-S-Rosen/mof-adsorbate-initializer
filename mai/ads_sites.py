import numpy as np
from ase import Atoms, Atom
from ase.io import write
from mai.regression import OLS_fit, TLS_fit
from mai.tools import string_to_formula
from mai.species_rules import add_monoatomic, add_diatomic, add_triatomic, add_CH4_SS
import os

"""
This module provides a class to identify ideal adsorption sites
"""
class ads_pos_optimizer():
	"""
	This identifies ideal adsorption sites
	"""
	def __init__(self,adsorbate_constructor,write_file=True,
		new_mofs_path=None,error_path=None):
		"""
		Initialized variables

		Args:
			adsorbate_constructor (class): adsorbate_constructor class
			containing many relevant defaults
			
			write_file (bool): if True, the new ASE atoms object should be
			written to a CIF file (defaults to True)
			
			new_mofs_path (string): path to store the new CIF files if
			write_file is True (defaults to /new_mofs)
			
			error_path (string): path to store any adsorbates flagged as
			problematic (defaults to /errors)
		"""
		self.full_ads_species = adsorbate_constructor.ads
		self.ads = adsorbate_constructor.ads.split('_')[0]
		self.d_MX1 = adsorbate_constructor.d_MX1
		self.overlap_tol = adsorbate_constructor.overlap_tol
		self.r_cut = adsorbate_constructor.r_cut
		self.sum_tol = adsorbate_constructor.sum_tol
		self.rmse_tol = adsorbate_constructor.rmse_tol
		self.site_idx = adsorbate_constructor.site_idx
		self.d_X1X2 = adsorbate_constructor.d_X1X2
		self.d_X2X3 = adsorbate_constructor.d_X2X3
		self.ang_MX1X2 = adsorbate_constructor.ang_MX1X2
		self.ang_triads = adsorbate_constructor.ang_triads
		self.eta = adsorbate_constructor.eta
		self.connect = adsorbate_constructor.connect
		self.start_atoms = adsorbate_constructor.start_atoms
		self.name = adsorbate_constructor.name
		self.write_file = write_file
		
		if new_mofs_path is None:
			new_mofs_path = os.path.join(os.getcwd(),'new_mofs')
		self.new_mofs_path = new_mofs_path
		if error_path is None:
			error_path = os.path.join(new_mofs_path,'errors')
		self.error_path = error_path

	def get_dist_planar(self,normal_vec):
		"""
		Get distance vector for planar adsorption site

		Args:
			normal_vec (numpy array): 1D numpy array for normal vector
		
		Returns:
			dist (float): distance vector scaled to d_MX1
		"""
		#Scale normal vector to desired bond length
		d_MX1 = self.d_MX1
		unit_normal = normal_vec/np.linalg.norm(normal_vec)
		dist = unit_normal*d_MX1

		return dist

	def get_NNs(self,ads_pos,site_idx):
		"""
		Get the number of atoms nearby the proposed adsorption site within r_cut

		Args:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
			
			site_idx (int): ASE index for adsorption site
		
		Returns:
			NN (int): number of neighbors within r_cut
			
			min_dist (float): distance from adsorbate to nearest atom
		"""
		r_cut = self.r_cut

		#Add proposed adsorbate to a copy of the MOF
		mof_temp = self.start_atoms.copy()
		adsorbate = Atoms([Atom('X',ads_pos)])
		mof_temp.extend(adsorbate)

		#Determine the number of nearby atoms and minimum distance
		compare_with = np.arange(0,len(mof_temp)-1).tolist()
		del compare_with[site_idx]
		neighbor_dist = mof_temp.get_distances(len(mof_temp)-1,compare_with,
			mic=True)
		NN = sum(neighbor_dist <= r_cut)
		min_dist = np.min(neighbor_dist)
		
		return NN, min_dist

	def get_planar_ads_pos(self,center_coord,dist,site_idx):
		"""
		Get adsorption site for planar structure

		Args:
			center_coord (numpy array): 1D numpy array for adsorption site
			(i.e. the central atom)
			
			site_idx (int): ASE index for adsorption site
		
		Returns:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		"""
		NN = []
		min_dist = []

		#Get +/- normal vector
		for i in range(2):
			if i == 0:
				ads_pos_temp = center_coord + dist
			elif i == 1:
				ads_pos_temp = center_coord - dist
			NN_temp, min_dist_temp = self.get_NNs(ads_pos_temp,site_idx)
			NN.append(NN_temp)
			min_dist.append(min_dist_temp)

		#Select best direction for normal vector
		if NN[0] == NN[1]:
			if min_dist[0] >= min_dist[1]:
				ads_pos = center_coord + dist
			else:
				ads_pos = center_coord - dist
		elif NN[0] <= NN[1]:
			ads_pos = center_coord + dist
		else:
			ads_pos = center_coord - dist

		return ads_pos

	def get_nonplanar_ads_pos(self,scaled_sum_dist,center_coord):
		"""
		Get adsorption site for non-planar structure

		Args:
			scaled_sum_dist (numpy array): 2D numpy array for the
			scaled Euclidean distance vectors between each coordinating
			atom and the central atom (i.e. the adsorption site)
			
			center_coord (numpy array): 1D numpy array for adsorption site
		
		Returns:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		"""

		#Sum up Euclidean vectors and scale to bond distance
		d_MX1 = self.d_MX1
		dist = d_MX1*scaled_sum_dist/np.linalg.norm(scaled_sum_dist)
		ads_pos = center_coord - dist

		return ads_pos

	def get_bi_ads_pos(self,normal_vec,center_coord,site_idx):
		"""
		Get adsorption site for a 2-coordinate site

		Args:
			normal_vec (numpy array): 1D numpy array for the
			normal vector to the line
			
			center_coord (numpy array): 1D numpy array for adsorption site
			
			site_idx (int): ASE index of adsorption site
		
		Returns:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		"""
		r_cut = self.r_cut
		overlap_tol = self.overlap_tol
		dist = self.get_dist_planar(normal_vec)

		#Prepare 2 temporary adsorbates
		try_angles = np.arange(0,360,10)
		ads_pos_temp_unrotated1 = center_coord + dist
		ads_pos_temp_unrotated2 = center_coord - dist
		ads_temp1 = Atoms([Atom('X',ads_pos_temp_unrotated1)])
		ads_temp2 = Atoms([Atom('X',ads_pos_temp_unrotated2)])
		mof_temp_orig = self.start_atoms.copy()

		#Rotate one of the adsorbates about the axis
		#defined by the two coordinating atoms
		for i, angle in enumerate(try_angles):

			#Add 2 overlapping atoms
			mof_temp = mof_temp_orig.copy()
			mof_temp.extend(ads_temp1)
			mof_temp.extend(ads_temp2)

			#Offset one of the atoms by small value
			#(just to be higher than machine epsilon)
			mof_temp[-1].position += 1e-6

			#Set the new angle
			mof_temp.set_angle(-1,site_idx,-2,angle)

			#Determine number of nearby atoms
			dist_mat = mof_temp.get_distances(-2,np.arange(0,
				len(mof_temp)-2).tolist(),mic=True)
			NNs = sum(dist_mat <= r_cut)
			n_overlap = sum(dist_mat <= overlap_tol)
			min_dist = np.min(dist_mat)

			#Select best option
			if i == 0:
				ads_pos = mof_temp[-2].position
				old_min_NNs = NNs
				old_min_dist = min_dist
				old_overlap = n_overlap
			elif n_overlap <= old_overlap and (NNs < old_min_NNs or
				(NNs == old_min_NNs and min_dist > old_min_dist)):
				ads_pos = mof_temp[-2].position
				old_min_NNs = NNs
				old_min_dist = min_dist
				old_overlap = n_overlap

		return ads_pos

	def get_tri_ads_pos(self,normal_vec,scaled_sum_dist,center_coord,site_idx):
		"""
		Get adsorption site for a 3-coordinate site

		Args:
			normal_vec (numpy array): 1D numpy array for the
			normal vector to the line
			
			scaled_sum_dist (numpy array): 2D numpy array for the
			scaled Euclidean distance vectors between each coordinating
			atom and the central atom (i.e. the adsorption site)
			
			center_coord (numpy array): 1D numpy array for adsorption site
			
			site_idx (int): ASE index of adsorption site
		
		Returns:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		"""
		dist = self.get_dist_planar(normal_vec)
		ads_pos_planar = self.get_planar_ads_pos(center_coord,dist,site_idx)
		NN_planar, min_dist_planar = self.get_NNs(ads_pos_planar,site_idx)
		ads_pos_nonplanar = self.get_nonplanar_ads_pos(scaled_sum_dist,center_coord)
		NN_nonplanar, min_dist_nonplanar = self.get_NNs(ads_pos_nonplanar,
			site_idx)

		#3-coordinate can be something like trigonal planar or T-shaped.
		#Use planar algorithm (normal vector) for trigonal planar structures
		#andd sum of Euclidean vectors for the other shapes
		if NN_planar == NN_nonplanar:
			if min_dist_planar >= min_dist_nonplanar:
				ads_pos = ads_pos_planar
			else:
				ads_pos = ads_pos_nonplanar
		elif NN_planar <= NN_nonplanar:
			ads_pos = ads_pos_planar
		else:
			ads_pos = ads_pos_nonplanar
			
		return ads_pos

	def get_opt_ads_pos(self,mic_coords,site_idx):
		"""
		Get the optimal adsorption site

		Args:
			mic_coords (numpy array): 2D numpy array for the
			coordinates of each coordinating atom using the central
			atom (i.e. adsorption site) as the origin
			
			site_idx (int): ASE index of adsorption site
		
		Returns:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		"""
		sum_tol = self.sum_tol
		rmse_tol = self.rmse_tol
		scale_factor = 2.0

		#Get coordinates of center atom and coordination number
		if mic_coords.ndim == 1:
			cnum = 1
		else:
			cnum = np.shape(mic_coords)[0]
		start_atoms = self.start_atoms.copy()
		center_coord = start_atoms[site_idx].position

		#Calculate relevant quantities based on coordination number
		if cnum == 1:
			#Define direction as colinear with center atom and
			#coordinating atom
			normal_vec = mic_coords
		elif cnum == 2:
			#Define direction based on line perpendicular to
			#the two coordinating atoms
			normal_vec = OLS_fit(mic_coords)
		elif cnum >= 3:
			#Calculate normal vector to best-fit plane and the
			#sum of Euclidean vectors
			scaled_mic_coords = mic_coords*scale_factor/np.linalg.norm(
				mic_coords,axis=1)[np.newaxis].T
			scaled_sum_dist = np.sum(scaled_mic_coords,axis=0)
			norm_scaled = np.linalg.norm(scaled_sum_dist)
			rmse, normal_vec = TLS_fit(mic_coords)

		#Get adsorption site based on coordination number
		if cnum == 1:
			dist = self.get_dist_planar(normal_vec)
			ads_pos = center_coord-dist
		elif cnum == 2:
			ads_pos = self.get_bi_ads_pos(normal_vec,center_coord,site_idx)
		elif cnum == 3 and norm_scaled > sum_tol:
			ads_pos = self.get_tri_ads_pos(normal_vec,scaled_sum_dist,center_coord,
				site_idx)
		elif norm_scaled <= sum_tol or rmse <= rmse_tol:
			dist = self.get_dist_planar(normal_vec)
			ads_pos = self.get_planar_ads_pos(center_coord,dist,site_idx)
		else:
			ads_pos = self.get_nonplanar_ads_pos(scaled_sum_dist,center_coord)

		return ads_pos

	def get_new_atoms(self,ads_pos,site_idx):
		"""
		Get new ASE atoms object with adsorbate from pymatgen analysis

		Args:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		
		Returns:
			new_mof (ASE Atoms object): new ASE Atoms object with adsorbate		
		"""
		full_ads_species = self.full_ads_species
		name = self.name
		full_name = name+'_site'+str(site_idx)
		new_name = full_name+'_'+full_ads_species
		mof = self.start_atoms.copy()
		new_mof = self.construct_mof(mof,ads_pos,site_idx)
		overlap = self.check_and_write(new_mof,new_name)
		if overlap:
			return None

		return new_mof

	def get_new_atoms_grid(self,site_pos,ads_pos):
		"""
		Get new ASE atoms object with adsorbate from energy grid

		Args:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		
		Returns:
			new_mof (ASE Atoms object): new ASE Atoms object with adsorbate		
		"""
		full_ads_species = self.full_ads_species
		ads = self.ads
		site_idx = self.site_idx
		name = self.name
		full_name = name+'_site'+str(site_idx)
		new_name = full_name+'_'+full_ads_species
		
		#Add molecule to structure
		mof = self.start_atoms.copy()
		if full_ads_species == 'CH4_grid':
			new_mof = self.construct_mof(mof,ads_pos,site_idx)
		else:
			raise ValueError('Unsupported adsorbate: '+ads)

		#Confirm no overlapping atoms and write file
		overlap = self.check_and_write(new_mof,new_name)
		if overlap:
			return None

		return new_mof

	def check_and_write(self,new_mof,new_name):
		"""
		Check for overlapping atoms and write CIF file

		Args:
			new_mof (ASE Atoms object): the new MOF-adsorbate complex

			new_name (string): the name of the new CIF file to write
		
		Returns:
			overlap (boolean): True or False for overlapping atoms		
		"""
		overlap_tol = self.overlap_tol
		write_file = self.write_file
		error_path = self.error_path
		ads = self.ads
		new_mofs_path = self.new_mofs_path

		n_new_atoms = len(string_to_formula(ads))
		overlap = False

		#Loop over each new atom and confirm no overlap
		for i in range(n_new_atoms):
			dist = new_mof.get_distances(-(i+1),
				np.arange(0,len(new_mof))[0:len(new_mof)-n_new_atoms],mic=True)
			if np.sum(dist <= overlap_tol) > 0:
				overlap = True
				if write_file:
					if not os.path.exists(error_path):
						os.makedirs(error_path)
					write(os.path.join(error_path,new_name+'.cif'),new_mof)
				break
		if overlap:
			return True
		else:
			if write_file:
				if not os.path.exists(new_mofs_path):
					os.makedirs(new_mofs_path)
				write(os.path.join(new_mofs_path,new_name+'.cif'),new_mof)
		return False


	def construct_mof(self,mof,ads_pos,site_idx):
		"""
		Construct the MOF-adsorbate complex

		Args:
			ads_pos_optimizer (class): see ads_sites.py for details

			ads (string): adsorbate species

			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		
		Returns:
			mof (ASE Atoms object): ASE Atoms object with adsorbate
		"""
		full_ads_species = self.full_ads_species
		ads = self.ads
		r_cut = self.r_cut
		overlap_tol = self.overlap_tol
		d_X1X2 = self.d_X1X2
		d_X2X3 = self.d_X2X3
		ang_MX1X2 = self.ang_MX1X2
		ang_triads = self.ang_triads
		eta = self.eta
		connect = self.connect

		n_new_atoms = len(string_to_formula(ads))

		#Construct MOF based on type of adsorbate
		if '_grid' in full_ads_species:
			if ads == 'CH4':
				new_mof = add_CH4_SS(mof,site_idx,ads_pos)
			else:
				raise ValueError('Unsupported species for grid method')
		else:
			if n_new_atoms == 1:
				new_mof = add_monoatomic(mof,ads,ads_pos)
			elif n_new_atoms == 2:
				new_mof = add_diatomic(mof,ads,ads_pos,site_idx,
					d_X1X2=d_X1X2,ang_MX1X2=ang_MX1X2,eta=eta,connect=connect,
					r_cut=r_cut,overlap_tol=overlap_tol)
			elif n_new_atoms == 3:
				new_mof = add_triatomic(mof,ads,ads_pos,site_idx,
					d_X1X2=d_X1X2,d_X2X3=d_X2X3,ang_MX1X2=ang_MX1X2,ang_triads=ang_triads,
					connect=connect,r_cut=r_cut,overlap_tol=overlap_tol)	
			else:
				raise ValueError('Too many atoms in adsorbate')

		new_mof.wrap()
		return new_mof