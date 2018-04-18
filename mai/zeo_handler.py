import numpy as np
import os

def get_n_oms(oms_path):
	"""
	Get the number of open metal sites from Zeo++ OMS file

	Args:
		oms_path (string): path to the Zeo++ .oms file
	Returns:
		n_OMS (float): number of open metal sites
	"""
	f = open(oms_path,'r')
	oms_file = f.read()
	n_OMS = int(oms_file.split('OMS=')[1].split('\n')[0])
	f.close()

	return n_OMS

def get_omsex_line(line):
	"""
	Read a line in the OMSEX file

	Args:
		line (string): line in OMSEX file
	Returns:
		oms_symbol (string): atomic symbol for the OMS
		cnum (int): coordination number for the OMS
		oms_coord (numpy array): array of coordinates for the OMS
		coords (numpy array): array of coordinates for the nearest neighbors
	"""
	oms_symbol = line.split('|')[0].strip()
	cnum = int(line.split('CNUM: ')[1].split('|')[0])
	oms_coords = np.asarray(np.matrix(
		line.split('COORD: ')[1].split('|')[0][0:-1]))
	coords = np.asarray(np.matrix(line.split('NN: ')[1][0:-4]))

	return oms_symbol, cnum, oms_coords, coords

def get_ase_NN_idx(atoms,coords):
	"""
	Get the ASE indices for the coordinating atoms

	Args:
		atoms (Atoms object): ASE Atoms object for the MOF
		coords (numpy array): coordinates of the coordinating atoms
	Returns:
		ase_NN_idx (list of ints): ASE indices of the coordinating atoms
	"""
	ase_NN_idx = []
	zeo_tol = 0.1
	for i in range(np.shape(coords)[0]):
		nn_fail = False
		for j, element in enumerate(atoms):
			if (sum(coords[i,:] >= element.position-zeo_tol) == 3 and 
			sum(coords[i,:] <= element.position+zeo_tol) == 3):
				ase_NN_idx.append(j)
				nn_fail = True
				break
		if nn_fail == False:
			raise ValueError('Zeo++ index does not match ASE index')

	return ase_NN_idx

def get_omsex_data(oms_data_path,name,atoms):
	"""
	Get info about the open metal site from Zeo++ OMSEX file using modified
	'network.cc' file

	Args:
		oms_data_path (string): path to the Zeo++ OMS and OMSEX data
		name (string): name of the MOF
		atoms (ASE Atoms object): Atoms object for the MOF
	Returns:
		omsex_dict (dict): dictionary of data from the Zeo++ OMSEX file
	"""
	if os.stat(os.path.join(oms_data_path,name+'.omsex')).st_size == 0:
		return None
	n_OMS = get_n_oms(os.path.join(oms_data_path,name+'.oms'))
	oms_coords_all = np.zeros((n_OMS,3))
	oms_idx_all = []
	oms_sym_all = []
	cnums_all = []
	NN_idx_all = []
	zeo_tol = 0.01
	with open(os.path.join(oms_data_path,name+'.omsex'),'r') as rf:
		for i, line in enumerate(rf):
			(oms_sym_temp, cnum_temp, oms_coords_all[i,:],
				NN_coords_temp) = get_omsex_line(line)
			oms_sym_all.append(oms_sym_temp)
			cnums_all.append(cnum_temp)
			if i == 0:
				NN_coords_all = NN_coords_temp
			else:
				NN_coords_all = np.vstack((NN_coords_all,NN_coords_temp))
			NN_idx_temp = get_ase_NN_idx(atoms,NN_coords_temp)
			NN_idx_all.append(NN_idx_temp)
			for j, element in enumerate(atoms):
				if (sum(oms_coords_all[i,:] >= element.position-zeo_tol) == 3
					and sum(oms_coords_all[i,:] <= element.position+zeo_tol) == 3):
					oms_idx_all.append(j)
					break
			if len(oms_idx_all) < i+1:
				raise ValueError('ERROR with '+name+': a zeo++ OMS (#'+
					str(i)+') is not in same spot as in ASE')
			if len(NN_idx_temp) != cnum_temp:
				raise ValueError('ASE/Zeo++ atom mismatch in NNs')
	omsex_dict = {'cnums':cnums_all,'oms_coords':oms_coords_all,
	'oms_idx':oms_idx_all,'oms_sym':oms_sym_all,'NN_coords':NN_coords_all,
	'NN_idx':NN_idx_all}

	return omsex_dict