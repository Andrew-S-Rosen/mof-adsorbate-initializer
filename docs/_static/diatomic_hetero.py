from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Ni-BTP.cif' #path to CIF of MOF

#add adsorbate
ads = adsorbate_constructor(ads='CO',d_MX1=1.5,site_idx=0,
	d_X1X2=1.13)
new_mof_atoms1 = ads.get_adsorbate(atoms_path=mof_path)

ads = adsorbate_constructor(ads='OC',d_MX1=1.5,site_idx=0,
	d_X1X2=1.13)
new_mof_atoms2 = ads.get_adsorbate(atoms_path=mof_path)
