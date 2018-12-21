from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Ni2Cl2-BTDD.cif' #path to CIF of MOF

#add adsorbate
ads = adsorbate_constructor(ads='N2O',d_MX1=2.0,site_idx=0,
	d_X1X2=1.13,d_X2X3=1.19)
new_mof_atoms1 = ads.get_adsorbate(atoms_path=mof_path)

ads = adsorbate_constructor(ads='ON2',d_MX1=2.0,site_idx=0,
	d_X1X2=1.19,d_X2X3=1.13)
new_mof_atoms2 = ads.get_adsorbate(atoms_path=mof_path)
