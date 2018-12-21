from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Ni-BTP.cif' #path to CIF of MOF

#add adsorbate
ads_species = 'O2_end' #adsorbate species
ads = adsorbate_constructor(ads='O2_end',d_MX1=1.5,site_idx=0,
	eta=1,d_X1X2=1.2,ang_MX1X2=120)
new_mof_atoms1 = ads.get_adsorbate(atoms_filepath=mof_path)

ads_species = 'O2_side' #adsorbate species
ads = adsorbate_constructor(ads='O2_side',d_MX1=1.5,site_idx=0,
	eta=2,d_X1X2=1.2,ang_MX1X2=90)
new_mof_atoms2 = ads.get_adsorbate(atoms_filepath=mof_path)
