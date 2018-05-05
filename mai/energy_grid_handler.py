import os
import pandas as pd
import numpy as np
from ase.build import molecule
from ase.geometry import get_distances

def read_grid(grid_filepath):
	"""
	Convert ASCII energy grid to pandas dataframe

	Args:
		grid_filepath (string): path to ASCII energy grid
	Returns:
		df (pandas df object): df containing energy grid details (x,y,z,E)
	"""
	if not os.path.isfile(grid_filepath):
		return None
	df = pd.read_csv(grid_filepath,delim_whitespace=True,na_values='?',
		usecols=[0,1,2,3])
	df.columns = ['x','y','z','E']
	
	return df

def grid_within_cutoff(df,atoms,max_dist,site_pos,partition=1e6):
	"""
	Reduces grid dataframe into data within max_dist of active site

	Args:
		df (pandas df object): df containing energy grid details (x,y,z,E)
		atoms (ASE Atoms object): Atoms object of structure
		max_dist (float): maximum distance from active site to consider
		site_pos (array): numpy array of the adsorption site
		partition (float): how many data points to partition the df for. This
		is used to prevent memory overflow errors. Decrease if memory errors
		arise.
	Returns:
		new_df (pandas df object): modified df only around max_dist from active
		site and also with a new distance (d) column
	"""
	df['d'] = ''
	n_loops = int(np.ceil(len(df)/partition))
	new_df = pd.DataFrame()

	#Only onsider data within max_dist of active site. Cut up the original
	#dataframe into chunks defined by partition to prevent memory issues
	for i in range(n_loops):
		if i == n_loops-1:
			idx = np.arange(i*int(partition),len(df))
		else:
			idx = np.arange(i*int(partition),(i+1)*int(partition))
		D,D_len = get_distances([site_pos],df.loc[idx,['x','y','z']].as_matrix(),
			cell=atoms.cell,pbc=atoms.pbc)
		D_len.shape = (-1,)
		df.loc[idx,'d'] = D_len
	new_df = df[df['d'] <= max_dist]

	return new_df

def get_best_grid_pos(atoms,max_dist,site_idx,grid_filepath):
	"""
	Finds minimum energy position in grid dataframe

	Args:
		atoms (ASE Atoms object): Atoms object of structure
		max_dist (float): maximum distance from active site to consider
		site_idx (int): ASE index of adsorption site
		grid_filepath (string): path to ASCII energy grid
	Returns:
		ads_pos (array): 1D numpy array for the ideal adsorption position
	"""
	df = read_grid(grid_filepath)
	if df is None:
		return 'nogrid'
	site_pos = atoms[site_idx].position
	cut_df = grid_within_cutoff(df,atoms,max_dist,site_pos)
	if np.sum(cut_df['E']) == 0:
		return 'invalid'
	best = cut_df.loc[cut_df['E'].idxmin()]
	ads_pos = [best['x'],best['y'],best['z']]

	return ads_pos

def add_CH4(site_idx,ads_pos,atoms):
	"""
	Add CH4 to the structure

	Args:
		site_idx (int): ASE index of site based on single-site model
		ads_pos (array): 1D numpy array for the best adsorbate position
		atoms (ASE Atoms object): Atoms object of structure
	Returns:
		atoms (ASE Atoms object): new ASE Atoms object with adsorbate
		n_new_atoms (int): number of atoms in adsorbate
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
	D,D_len = get_distances([ads_pos],atoms[site_idx].position,cell=atoms.cell,pbc=atoms.pbc)
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
	atoms.extend(CH4)
	atoms.wrap()
	
	return atoms, len(CH4)

def add_N2(site_idx,ads_pos,atoms):
	"""
	Add N2 to the structure based on single-site model

	Args:
		site_idx (int): ASE index of site
		ads_pos (array): 1D numpy array for the best adsorbate position
		atoms (ASE Atoms object): Atoms object of structure
	Returns:
		atoms (ASE Atoms object): new ASE Atoms object with adsorbate
		n_new_atoms (int): number of atoms in adsorbate
	"""
	#Get N2 parameters
	N2 = molecule('N2')
	NN_length = N2.get_distance(0,1)

	D,D_len = get_distances([ads_pos],atoms[site_idx].position,cell=atoms.cell,pbc=atoms.pbc)
	r_vec = D[0,0]
	r_N = (r_vec/np.linalg.norm(r_vec))*NN_length

	#Construct N2
	N2[0].position = ads_pos-r_N/2
	N2[1].position = ads_pos+r_N/2

	#Add N2 molecule to the structure
	if site_idx is None:
		raise ValueError('Site index must not be None')
	atoms.extend(N2)
	atoms.wrap()

	return atoms, len(N2)