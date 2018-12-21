import os
from ase.io import read
import numpy as np
success_basepath = os.path.join('test','success')
test_basepath = os.path.join('examples','new_mofs')
jobs = ['add_O','add_O2','add_CO','add_N2O','add_H2O','add_all_H2O','add_CH4_PEG','add_O_auto',]
tol = 1E-3
for job in jobs:
	success_path = os.path.join(success_basepath,job)
	for file in os.listdir(success_path):
		if '.cif' in file:
			mof_test = read(os.path.join(test_basepath,file))
			mof_real = read(os.path.join(success_path,file))
			diff = np.abs(mof_real.get_positions()-mof_test.get_positions())
			if np.sum(diff >= tol) != 0:
				print('MOF (real):')
				print(mof_real.get_positions())
				print('MOF (test):')
				print(mof_test.get_positions())
				print('Difference:')
				print(diff)
				raise ValueError('Error with: '+file)
			else:
				print('SUCCESS for '+job)
			break