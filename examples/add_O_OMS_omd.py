from omsdetector import MofCollection
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'example_MOFs' #path to CIFs
analysis_path = mof_path #path to store the OMS results
ads_species = 'O' #adsorbate species
bond_length = 2.0 #desired distance between OMS and ads_species

#Run the Open Metal Detector
mof_coll = MofCollection.from_folder(collection_folder=mof_path,analysis_folder=analysis_path)
mof_coll.analyse_mofs()

#add adsorbate
ads = adsorbate_constructor(ads_species,bond_length)
new_mof_atoms, new_mof_name = ads.get_adsorbate_oms(mof_path,oms_format='OMD')