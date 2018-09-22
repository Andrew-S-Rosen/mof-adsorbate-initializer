import os
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = os.path.join('example_MOFs','AHOKIR01-O.cif') #path to CIF of MOF
adsorbate = 'CH4' #adsorbate species
site = 'O' #adsorption site
max_dist = 3.0 #maximum distance between O and C
ads = adsorbate_constructor(adsorbate,max_dist,site_species=site)
mof_adsorbate, mof_name = ads.get_adsorbate_grid(mof_path)