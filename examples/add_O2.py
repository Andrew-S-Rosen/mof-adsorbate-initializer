import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','Ni-BTP.cif') #path to CIF of MOF

#add adsorbate
ads = adsorbate_constructor(ads='O2_end',d_MX1=1.75,eta=1)
new_mof_atoms = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)

ads = adsorbate_constructor(ads='O2_side',d_MX1=1.75,eta=2)
new_mof_atoms = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)
