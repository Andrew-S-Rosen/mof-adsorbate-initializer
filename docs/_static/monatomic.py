from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'Cu-BTC.cif' #path to CIF of MOF

#add adsorbate
ads = adsorbate_constructor(ads='O',d_MX1=1.75,site_idx=0)
new_mof_atoms = ads.get_adsorbate(atoms_path=mof_path)