import os
from mai.adsorbate_constructor import adsorbate_constructor
from ase.build import molecule
from ase.io import write
from mai.janitor import get_refcode
import numpy as np

#NOTE: Still needs to be commented

N2O = molecule('N2O')
NN_length = N2O.get_distance(0,1)
NO_length = N2O.get_distance(1,2)
mof_path = 'oxygenated_MOFs'
N2_mofs_path = 'add_N2/'
N2O_mofs_path = 'add_N2O'
max_dist = 3.0
overlap_tol = 1.3
mol_species = 'N2'
site_species = 'O'

if not os.path.exists(N2O_mofs_path):
	os.makedirs(N2O_mofs_path)
filenames = os.listdir(mof_path)
filenames.reverse()
nonmetal_list = [1,2,6,7,8,9,10,15,16,17,18,34,35,36,53,54,86]
MO_min_dist = 2.75
for filename in filenames:
	filepath = os.path.join(mof_path,filename)
	ads_const = adsorbate_constructor(mol_species,max_dist,
		site_species=site_species,overlap_tol=overlap_tol)
	mof_adsorbate, mof_name = ads_const.get_adsorbate_raspa(filepath,
		new_mofs_path=N2_mofs_path)
	if mof_adsorbate is None:
		continue
	refcode = get_refcode(filepath)
	O_idx = [atom.index for atom in mof_adsorbate if atom.symbol == site_species][-1]
	N_dist1 = mof_adsorbate.get_distance(O_idx,-1,mic=True,vector=False)
	N_dist2 = mof_adsorbate.get_distance(O_idx,-2,mic=True,vector=False)
	if N_dist1 < N_dist2:
		N1_idx = -1
		N2_idx = -2
	else:
		N1_idx = -2
		N2_idx = -1
	vec = mof_adsorbate.get_distance(O_idx,N1_idx,mic=True,vector=True)
	mag_vec = np.linalg.norm(vec)
	norm_vec = vec/mag_vec
	mof_adsorbate[O_idx].position += norm_vec*0.75
	mof_adsorbate.set_distance(O_idx,N1_idx,NO_length,fix=0,mic=True)
	mof_adsorbate.set_distance(N1_idx,N2_idx,NN_length,fix=0,mic=True)
	compare_with = np.arange(0,len(mof_adsorbate)-2).tolist()
	neighbor_dist = mof_adsorbate.get_distances(O_idx,compare_with,
		mic=True)
	NN = sum(neighbor_dist < MO_min_dist)-1
	sorted_neighbors = np.argsort(neighbor_dist).tolist()
	del sorted_neighbors[0]
	min_indices = sorted_neighbors[0:NN]
	indices = [O_idx,N1_idx,N2_idx]
	flag=False
	if NN > 0:
		for min_idx in min_indices:
			if mof_adsorbate[min_idx].number not in nonmetal_list:
				vec = mof_adsorbate.get_distance(min_idx,O_idx,mic=True,vector=True)
				mag_vec = np.linalg.norm(vec)
				norm_vec = vec/mag_vec
				for idx in indices:
					mof_adsorbate[idx].position += norm_vec*(MO_min_dist-mag_vec)
	mof_adsorbate.wrap()
	compare_with.remove(O_idx)
	for idx in indices:
		neighbor_dist = mof_adsorbate.get_distances(idx,compare_with,
			mic=True)
		overlap = sum(neighbor_dist <= overlap_tol)
		if overlap > 0:
			print('ERROR: '+refcode)
			break
	write(os.path.join(N2O_mofs_path,refcode+'_N2O.cif'),mof_adsorbate)