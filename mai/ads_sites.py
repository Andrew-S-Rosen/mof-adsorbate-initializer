import numpy as np
from ase import Atoms, Atom
from ase.io import read, write
from ase.build import molecule
from mai.regression import OLS_fit, TLS_fit
from mai.janitor import get_refcode
from mai.gases import add_CH4_SS, add_N2
import os

"""
This module provides a class to identify ideal adsorption sites
"""
class ads_pos_optimizer():
	"""
	This identifies ideal adsorption sites
	"""
	def __init__(self,adsorbate_constructor,atoms_filepath,write_file=True,
		new_mofs_path=None,error_path=None):
		"""
		Initialized variables

		Args:
			adsorbate_constructor (class): adsorbate_constructor class
			containing many relevant defaults

			atoms_filepath (string): path to the structure file
			
			write_file (bool): if True, the new ASE atoms object should be
			written to a CIF file (defaults to True)
			
			new_mofs_path (string): path to store the new CIF files if
			write_file is True (defaults to atoms_filepath/new_mofs)
			
			error_path (string): path to store any adsorbates flagged as
			problematic (defaults to atoms_filepath/errors)
		"""
		self.ads_species = adsorbate_constructor.ads_species
		self.bond_dist = adsorbate_constructor.bond_dist
		self.overlap_tol = adsorbate_constructor.overlap_tol
		self.r_cut = adsorbate_constructor.r_cut
		self.sum_tol = adsorbate_constructor.sum_tol
		self.rmse_tol = adsorbate_constructor.rmse_tol
		self.site_idx = adsorbate_constructor.site_idx
		self.atoms_filepath = atoms_filepath
		self.write_file = write_file

		if new_mofs_path is None:
			new_mofs_path = os.path.join(os.path.dirname(atoms_filepath),
				'new_mofs')
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
			dist (float): distance vector scaled to bond_dist
		"""
		#Scale normal vector to desired bond length
		bond_dist = self.bond_dist
		unit_normal = normal_vec/np.linalg.norm(normal_vec)
		dist = unit_normal*bond_dist

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
		ads_species = self.ads_species
		r_cut = self.r_cut
		atoms_filepath = self.atoms_filepath

		#Add proposed adsorbate
		mof_temp = read(atoms_filepath)
		adsorbate = Atoms([Atom(ads_species,ads_pos)])
		mof_temp.extend(adsorbate)

		#Determine the number of nearby atoms and minimum distance
		compare_with = np.arange(0,len(mof_temp)-1).tolist()
		del compare_with[site_idx]
		neighbor_dist = mof_temp.get_distances(len(mof_temp)-1,compare_with,
			mic=True)
		NN = sum(neighbor_dist <= r_cut)
		min_dist = np.min(neighbor_dist)

		return NN, min_dist

	def get_best_to_worst_idx(self,ads_poss,site_idx_list):
		"""
		Sort the potential adsorption sites by best to worst

		Args:
			ads_poss (numpy array): 2D numpy array for the proposed
			adsorption positions
			
			site_idx_list (list of ints): ASE indices for adsorption sites
		
		Returns:
			best_to_worst_idx (list of ints): sorted adsorption sites from best
			to worst
		"""
		NN = []
		min_dist = []
		i_vec = []
		best_to_worst_idx = []
		if len(site_idx_list) != np.shape(ads_poss)[0]:
			raise ValueError('Incompatible lengths of lists')

		#Cycle through proposed adsorbates sort by best
		for i, ase_cus_idx in enumerate(site_idx_list):
			NN_temp, min_dist_temp = self.get_NNs(ads_poss[i,:],ase_cus_idx)
			NN.append(NN_temp)
			min_dist.append(min_dist_temp)
			i_vec.append(i)
		merged_list = list(zip(i_vec,NN,min_dist))
		merged_list.sort(key=lambda x: x[2],reverse=True)
		merged_list.sort(key=lambda x: x[1])
		for item in merged_list:
			best_to_worst_idx.append(item[0])

		return best_to_worst_idx

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
		bond_dist = self.bond_dist
		dist = bond_dist*scaled_sum_dist/np.linalg.norm(scaled_sum_dist)
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
		bond_dist = self.bond_dist
		r_cut = self.r_cut
		overlap_tol = self.overlap_tol
		ads_species = self.ads_species
		atoms_filepath = self.atoms_filepath
		dist = self.get_dist_planar(normal_vec)

		#Prepare 2 overlapping adsorbates
		try_angles = np.arange(0,360,10)
		ads_pos_temp_unrotated1 = center_coord + dist
		ads_pos_temp_unrotated2 = center_coord - dist
		ads_temp1 = Atoms([Atom(ads_species,ads_pos_temp_unrotated1)])
		ads_temp2 = Atoms([Atom(ads_species,ads_pos_temp_unrotated2)])
		mof_temp_orig = read(atoms_filepath)

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
			mof_temp.set_distance(site_idx,len(mof_temp)-1,bond_dist,
				fix=0,mic=True)
			mof_temp.set_distance(site_idx,len(mof_temp)-2,bond_dist,
				fix=0,mic=True)	
			mof_temp.set_angle(len(mof_temp)-1,site_idx,len(mof_temp)-2,angle)

			#Determine number of nearby atoms
			dist_mat = mof_temp.get_distances(len(mof_temp)-2,
				np.arange(0,len(mof_temp)-2).tolist(),mic=True)
			NNs = sum(dist_mat <= r_cut)

			#Select best option
			if i == 0:
				ads_pos = mof_temp[-2].position
				old_min_NNs = NNs
			elif sum(dist_mat <= overlap_tol) == 0 and NNs < old_min_NNs:
				ads_pos = mof_temp[-2].position
				old_min_NNs = NNs

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
		atoms_filepath = self.atoms_filepath
		scale_factor = 2.0

		#Get coordinates of center atom and coordination number
		if mic_coords.ndim == 1:
			cnum = 1
		else:
			cnum = np.shape(mic_coords)[0]
		center_coord = read(atoms_filepath)[site_idx].position

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

	def add_ads_species(self,ads_pos):
		"""
		Add adsorbate to the ASE atoms object

		Args:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		
		Returns:
			new_mof (ASE Atoms object): Atoms object for new
			structure with adsorbate
		"""
		ads_species = self.ads_species
		atoms_filepath = self.atoms_filepath
		new_mof = read(atoms_filepath)
		adsorbate = Atoms([Atom(ads_species,ads_pos)])
		new_mof.extend(adsorbate)
		new_mof.wrap()
		
		return new_mof

	def get_new_atoms_zeo_oms(self,ads_poss,best_to_worst_idx,cluster):
		"""
		Construct new Atoms object with adsorbate based on Zeo++ OMS
		detection

		Args:
			ads_poss (numpy array): 2D numpy array for the proposed
			adsorption positions
			
			best_to_worst_idx (list of ints): sorted adsorption sites from best
			to worst
			
			cluster (list of ints): atomic numbers for each atom in coordination
			environment (useful for debugging) 
		
		Returns:
			new_mof (ASE Atoms object): new ASE Atoms object with adsorbate
			
			name (string): name of new structure with adsorbate
		"""
		overlap_tol = self.overlap_tol
		ads_species = self.ads_species
		write_file = self.write_file
		new_mofs_path = self.new_mofs_path
		error_path = self.error_path
		atoms_filepath = self.atoms_filepath

		atoms_filename = os.path.basename(atoms_filepath)
		name = get_refcode(atoms_filename)
		basename = name+'_'+ads_species
		success = False

		#Cycle through all proposed adsorbates and write the best
		#one that does not overlap with other atoms
		for idx in best_to_worst_idx:
			new_mof = self.add_ads_species(ads_poss[idx,:])
			dist_mat = new_mof.get_distances(len(new_mof)-1,
				np.arange(0,len(new_mof)-1).tolist(),mic=True)
			
			if sum(dist_mat <= overlap_tol) == 0:
				new_name = basename+'_OMS'+str(idx)
				if write_file == True:
					write(os.path.join(new_mofs_path,new_name+'.cif'),new_mof)
				success = True
				break
			else:
				del new_mof[-1]

		if success == False:
			if write_file == True:
				write(os.path.join(error_path,basename+'_'+str(cluster)+'.cif'),
					new_mof)
			new_mof = None
			new_name = None

		return new_mof, new_name

	def get_new_atoms(self,ads_pos,site_idx):
		"""
		Get new ASE atoms object with adsorbate

		Args:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		
		Returns:
			new_mof (ASE Atoms object): new ASE Atoms object with adsorbate
		
			name (string): name of new structure with adsorbate
		"""
		atoms_filepath = self.atoms_filepath
		ads_species = self.ads_species
		overlap_tol = self.overlap_tol
		write_file = self.write_file
		new_mofs_path = self.new_mofs_path
		error_path = self.error_path

		atoms_filename = os.path.basename(atoms_filepath)
		name = get_refcode(atoms_filename)
		new_name = name+'_'+ads_species

		new_mof = self.add_ads_species(ads_pos)
		compare_to = np.arange(0,len(new_mof)-1).tolist()
		compare_to.remove(site_idx)
		dist_mat = new_mof.get_distances(len(new_mof)-1,compare_to,mic=True)

		if sum(dist_mat <= overlap_tol) == 0:
			if write_file == True:
				write(os.path.join(new_mofs_path,new_name+'.cif'),new_mof)
		else:
			if write_file == True:
				write(os.path.join(error_path,new_name+'.cif'),new_mof)
			return None, None

		return new_mof, new_name

	def get_new_atoms_grid(self,site_pos,ads_pos):
		"""
		Get new ASE atoms object with adsorbate from energy grid

		Args:
			ads_pos (numpy array): 1D numpy array for the proposed
			adsorption position
		
		Returns:
			new_mof (ASE Atoms object): new ASE Atoms object with adsorbate
		
			name (string): name of new structure with adsorbate
		"""
		atoms_filepath = self.atoms_filepath
		ads_species = self.ads_species
		overlap_tol = self.overlap_tol
		write_file = self.write_file
		new_mofs_path = self.new_mofs_path
		error_path = self.error_path
		site_idx = self.site_idx

		#Add molecule to structure
		atoms_filename = os.path.basename(atoms_filepath)
		mof = read(atoms_filepath)
		name = get_refcode(atoms_filename)
		new_name = name+'_'+ads_species
		if ads_species == 'CH4':
			new_mof, n_new_atoms = add_CH4_SS(site_idx,ads_pos,mof)
		elif ads_species == 'N2':
			N2 = molecule('N2')
			NN_length = N2.get_distance(0,1)
			new_mof, n_new_atoms = add_N2(site_idx,ads_pos,mof)
		else:
			raise ValueError('Unsupported molecular adsorbate')

		#Confirm no overlapping atoms and write file
		overlap = False
		for i in range(n_new_atoms):
			dist = new_mof.get_distances(-(i+1),np.arange(0,n_new_atoms).tolist(),
				mic=True)
			if np.sum(dist <= overlap_tol) > 0:
				overlap = True
				if write_file == True:
					write(os.path.join(error_path,name+'_'+ads_species+'.cif'),new_mof)
				break
		if overlap == True:
			return None, None
		else:
			if write_file == True:
				write(os.path.join(new_mofs_path,new_name+'.cif'),new_mof)

		return new_mof, new_name
