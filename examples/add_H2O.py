import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','Ni2Cl2-BBTA.cif') #path to CIF of MOF

#add adsorbate
ads_species = 'HOH' #adsorbate species (eta1-O)
ads = adsorbate_constructor(ads='HOH',d_MX1=2.0,
	d_X1X2=0.96,d_X2X3=0.96,ang_MX1X2=120,ang_triads=104.5,connect=2)
new_mof_atoms = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)
