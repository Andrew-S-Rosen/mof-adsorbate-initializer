from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Ni2Cl2-BTDD.cif' #path to CIF of MOF

#add H2O adsorbate
ads = adsorbate_constructor(ads='HOH',d_MX1=2.0,d_X1X2=0.96,d_X2X3=0.96,
	ang_MX1X2=120.0,ang_triads=104.5,connect=2)
new_mof_atoms = ads.get_adsorbate(atoms_path=mof_path,site_idx=0)
