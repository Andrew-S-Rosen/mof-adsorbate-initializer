import os
import pandas as pd
import numpy as np
from ase.build import molecule
from ase.geometry import get_distances

def read_grid(grid_filepath):
	"""
	Convert energy grid to pandas dataframe

	Args:
		grid_filepath (string): path to energy grid (must be .cube or .grid)

	Returns:
		df (pandas df object): df containing energy grid details (x,y,z,E)
	"""
	if not os.path.isfile(grid_filepath):
		return None
	grid_ext = grid_filepath.split('.')[-1]
	if grid_ext == 'cube':
		df = cube_to_xyzE(grid_filepath)
	elif grid_ext == 'grid':
		df = pd.read_csv(grid_filepath,delim_whitespace=True,na_values='?',
			header=None,usecols=[0,1,2,3])
	else:
		raise ValueError('Unsupported grid format for file '+grid_filepath)
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

	#Only consider data within max_dist of active site. Cut up the original
	#dataframe into chunks defined by partition to prevent memory issues
	for i in range(n_loops):
		if i == n_loops-1:
			idx = np.arange(i*int(partition),len(df))
		else:
			idx = np.arange(i*int(partition),(i+1)*int(partition))
		D,D_len = get_distances([site_pos],df.loc[idx,['x','y','z']].values,
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
		
		grid_filepath (string): path to energy grid

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

def cube_to_xyzE(cube_file):
	"""
	Converts cube to ASCII file
	Adopted from code by Julen Larrucea
	Original source: https://github.com/julenl/molecular_modeling_scripts

	Args:
		cube_file (string): path to cube file
	
	Returns:
		pd_data (Pandas dataframe): dataframe of (x,y,z,E) grid
	"""
	at_coord=[]
	spacing_vec=[]
	nline = 0
	values=[]
	data = []
	with open(cube_file,'r') as rw:
		for line in rw:
			nline += 1
			if nline == 3:
				nat=int(line.split()[0])
			elif nline >3 and nline <= 6:
				spacing_vec.append(line.split())
			elif nline > 6 and nline <= 6+nat:
				at_coord.append(line.split())
			elif nline > 5 and nline > 6+nat:
				for i in line.split():
					values.append(float(i))
	idx = -1
	for i in range(0,int(spacing_vec[0][0])):
		for j in range(0,int(spacing_vec[1][0])):
			for k in range(0,int(spacing_vec[2][0])):
				idx += 1
				x,y,z = i*float(spacing_vec[0][1]),j*float(spacing_vec[1][2]),k*float(spacing_vec[2][3])
				data.append([x,y,z,values[idx]])
	pd_data = pd.DataFrame(data)
	return pd_data