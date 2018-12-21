import os
from omsdetector import MofCollection
from mai.adsorbate_constructor import adsorbate_constructor

mofs_path = 'example_MOFs' #path to folder of CIFs
oms_data_path = os.getcwd() #path to store the OMS results

#Run the Open Metal Detector
mof_coll = MofCollection.from_folder(collection_folder=mofs_path,
	analysis_folder=oms_data_path)
mof_coll.analyse_mofs()

#add adsorbate for every CIF in mofs_path
for f in os.listdir(mofs_path):
	if '.cif' not in f:
		continue
	mof_path = os.path.join(mofs_path,f)
	ads = adsorbate_constructor(ads='O',d_MX1=1.75)
	new_mof_atoms = ads.get_adsorbate(atoms_path=mof_path)
