import os
import pandas as pd
import numpy as np
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