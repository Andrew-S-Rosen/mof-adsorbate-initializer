import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','Ni-BTP.cif') #path to CIF of MOF

#add O2 adsorbate in η1-O mode
ads = adsorbate_constructor(ads='O2_end',d_MX1=1.5,eta=1,
	d_X1X2=1.2,ang_MX1X2=120)
new_mof_atoms1 = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)

#add O2 adsorbate in η2-O mode
ads = adsorbate_constructor(ads='O2_side',d_MX1=1.5,eta=2,
	d_X1X2=1.2,ang_MX1X2=90)
new_mof_atoms2 = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)
