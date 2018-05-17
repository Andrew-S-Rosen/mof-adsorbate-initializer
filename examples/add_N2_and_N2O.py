import os
from mai.adsorbate_constructor import adsorbate_constructor
from ase.build import molecule
from ase.io import write
from mai.janitor import get_refcode

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
for filename in os.listdir(mof_path):
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
	mof_adsorbate.set_distance(O_idx,N1_idx,max_dist-0.5,fix=0,mic=True)
	mof_adsorbate.set_distance(N1_idx,O_idx,NO_length,fix=0,mic=True)
	mof_adsorbate.set_distance(N1_idx,N2_idx,NN_length,fix=0,mic=True)
	write(os.path.join(N2O_mofs_path,refcode+'_N2O.cif'),mof_adsorbate)