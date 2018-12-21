import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','Sc-MIL-88B.cif') #path to CIF of MOF

#add adsorbate
ads_species = 'N2O' #adsorbate species (eta1-N)
ads = adsorbate_constructor(ads='N2O',d_MX1=2.0,d_X1X2=1.13,
	d_X2X3=1.19)
new_mof_atoms = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)

ads_species = 'ON2' #adsorbate species (eta1-O)
ads = adsorbate_constructor(ads='ON2',d_MX1=2.0,d_X1X2=1.13,
	d_X2X3=1.19)
new_mof_atoms = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)
