import os
from omsdetector import MofCollection
from mai.adsorbate_constructor import adsorbate_constructor

mofs_path = 'example_MOFs' #path to folder of CIFs
oms_analysis_folder = os.getcwd() #path to store the OMS results
oms_data_path = os.path.join(oms_analysis_folder,'oms_results') #path to oms_results folder

#Run the Open Metal Detector
mof_coll = MofCollection.from_folder(collection_folder=mofs_path,
	analysis_folder=oms_data_path)
mof_coll.analyse_mofs()

#add adsorbate for every CIF in mofs_path
for file in os.listdir(mofs_path):
	if '.cif' not in file:
		continue
	mof_path = os.path.join(mofs_path,file)
	ads = adsorbate_constructor(ads='O',d_MX1=1.75)
	new_mof_atoms = ads.get_adsorbate(atoms_path=mof_path,
		omd_path=oms_data_path)
