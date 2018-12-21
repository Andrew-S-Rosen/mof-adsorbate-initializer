import os
from pymatgen.io.cif import CifParser
from pymatgen.analysis.structure_matcher import StructureMatcher

success_basepath = os.path.join('test','success')
test_basepath = os.path.join('examples','new_mofs')
jobs = ['add_O','add_O2','add_CO','add_N2O','add_H2O','add_all_H2O','add_CH4_PEG','add_O_auto',]
tol = 1E-4
for job in jobs:
	success_path = os.path.join(success_basepath,job)
	for file in os.listdir(success_path):
		if '.cif' in file:
			mof_test = CifParser(os.path.join(test_basepath,file)).get_structures(primitive=False)[0]
			mof_real = CifParser(os.path.join(success_path,file)).get_structures(primitive=False)[0]
			sm = StructureMatcher(primitive_cell=False)
			rms = sm.get_rms_dist(mof_test,mof_real)
			if not rms or rms[0] > tol:
				print(rms[0])
				raise ValueError('Error with: '+file)
			else:
				print('SUCCESS for '+job)
			break