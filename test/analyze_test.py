import os
from ase.io import read
import numpy as np

success_basepath = 'test/success/add_'
test_basepath = 'examples/add_'
species = ['O','H','CH4']
tol = 1E-3
for specie in species:
	path = test_basepath+specie+'/'
	success_path = success_basepath+specie+'/'
	for file in os.listdir(success_path):
		if not os.path.isfile(os.path.join(path,file)):
			continue
		if '.png' in file:
			continue
		mof_test = read(os.path.join(path,file))
		mof_real = read(os.path.join(success_path,file))
		diff = np.abs(mof_real.get_positions()-mof_test.get_positions())
		if np.sum(diff >= tol) != 0:
			print(diff)
			raise ValueError('Error with: '+file)