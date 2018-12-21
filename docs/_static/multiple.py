import os
from mai.adsorbate_constructor import adsorbate_constructor
from ase.io import read, write

starting_mof_path = os.path.join('example_MOFs','Ni2Cl2-BBTA.cif') #path to CIF of MOF

start_mof = read(starting_mof_path)
Ni_idx = [atom.index for atom in start_mof if atom.symbol == 'Ni']
ads_species = 'HOH' #adsorbate species (eta1-O)

#add adsorbate
atoms = start_mof
for i, site_idx in enumerate(Ni_idx):
	ads = adsorbate_constructor(ads='HOH',d_MX1=2.0,d_X1X2=0.96,d_MX2X3=0.96,
		ang_MX1X2=120,ang_triads=104.5,connect=2)
	atoms = ads.get_adsorbate(atoms=atoms,site_idx=site_idx,write_file=False)

write('Ni2Cl2-BBTA_allH2O.cif',atoms)