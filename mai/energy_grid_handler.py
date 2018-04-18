import os
import pandas as pd
import numpy as np
from ase.build import molecule
from ase.geometry import get_distances

def read_grid(grid_filepath):
	if not os.path.isfile(grid_filepath):
		return None
	df = pd.read_csv(grid_filepath,delim_whitespace=True,na_values='?',
		usecols=[0,1,2,3])
	df.columns = ['x','y','z','E']
	
	return df

def grid_within_cutoff(df,atoms,max_dist,site_pos,partition=1e6):
	df['d'] = ''
	n_loops = int(np.ceil(len(df)/partition))
	new_df = pd.DataFrame()
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
	df = read_grid(grid_filepath)
	if df is None:
		return None
	site_pos = atoms[site_idx].position
	cut_df = grid_within_cutoff(df,atoms,max_dist,site_pos)
	best = cut_df.loc[cut_df['E'].idxmin()]
	ads_pos = [best['x'],best['y'],best['z']]

	return ads_pos

def add_CH4(site_pos,ads_pos,atoms):
	CH4 = molecule('CH4')
	CH_length = CH4.get_distance(0,1)
	CH_angle = CH4.get_angle(1,0,2)
	CH2_dihedral = CH4.get_dihedral(2,1,0,4)
	CH_length = CH4.get_distance(0,1)
	CH_angle = CH4.get_angle(1,0,2)
	CH_dihedral = CH4.get_dihedral(2,1,0,4)
	CH4[0].position = ads_pos
	D,D_len = get_distances([ads_pos],[site_pos],cell=atoms.cell,pbc=atoms.pbc)
	r_vec = D[0,0]
	r = (r_vec/np.linalg.norm(r_vec))*CH_length
	CH4[1].position = ads_pos+r
	CH4.set_distance(0,2,CH_length,fix=0)
	CH4.set_angle(1,0,2,CH_angle)
	CH4.set_distance(0,3,CH_length,fix=0)
	CH4.set_angle(1,0,3,CH_angle)
	CH4.set_dihedral(2,1,0,3,-CH_dihedral)
	CH4.set_distance(0,4,CH_length,fix=0)
	CH4.set_angle(1,0,4,CH_angle)
	CH4.set_dihedral(2,1,0,4,CH2_dihedral)
	atoms.extend(CH4)

	return atoms