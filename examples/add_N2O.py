import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','Ni2Cl2-BTDD.cif') #path to CIF of MOF

#add N2O adsorbate in η1-N mode
ads = adsorbate_constructor(ads='N2O',d_MX1=2.0,d_X1X2=1.13,d_X2X3=1.19)
new_mof_atoms1 = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)

#add N2O adsorbate in η1-O mode
ads = adsorbate_constructor(ads='ON2',d_MX1=2.0,d_X1X2=1.19,d_X2X3=1.13)
new_mof_atoms2 = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)
