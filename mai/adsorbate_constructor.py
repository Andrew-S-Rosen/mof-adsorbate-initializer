import os
import numpy as np
from ase.io import read
from ase.atoms import Atoms
from mai.ads_sites import ads_pos_optimizer
from mai.tools import get_refcode
from mai.oms_handler import get_omd_data
from mai.NN_algos import get_NNs_pm
from mai.grid_handler import get_best_grid_pos
"""
This module provides classes to add adsorbates to a MOF
"""
class adsorbate_constructor():
	"""
	This class constructs an ASE atoms object with an adsorbate
	Initialized variables

	Args:
		ads (string): string of element or molecule for adsorbate
		(defaults to 'X')

		d_MX1 (float): distance between adsorbate and surface atom. If
		used with get_adsorbate_grid, it represents the maximum distance
		(defaults to 2.0)

		eta (int): denticity of end-on (1) or side-on (2) (defaults to 1)

		connect (int): the connecting atom in the species string (defaults to 1)

		d_X1X2 (float): X1-X2 bond length (defaults to 1.25)

		d_X2X3 (float): X2-X3 bond length for connect == 1 or
		X1-X3 bond length for connect == 2 (defaults to d_bond1)

		ang_MX1X2 (float): site-X1-X2 angle (for diatomics, defaults to 180 degrees
		except for side-on in which it defaults to 90 or end-on O2 in which
		it defaults to 120; for triatomics, defaults to 180 except for H2O
		in which it defaults to 104.5)

		ang_triads (float): X3-X1-X2 angle (defaults to 180 degrees for connect == 1
		and 90 degrees for connect == 2)

		r_cut (float): cutoff distance for calculating nearby atoms when
		ranking adsorption sites

		sum_tol (float): threshold to determine planarity. when the sum
		of the Euclidean distance vectors of coordinating atoms is less
		than sum_tol, planarity is assumed
		
		rmse_tol (float): second threshold to determine planarity. when the 
		root mean square error of the best-fit plane is less than rmse_tol,
		planarity is assumed
		
		overlap_tol (float): distance below which atoms are assumed to be
		overlapping
		"""

	def __init__(self,ads='X',d_MX1=2.0,eta=1,connect=1,
		d_X1X2=1.25,d_X2X3=None,ang_MX1X2=None,ang_triads=None,
		r_cut=2.5,sum_tol=0.5,rmse_tol=0.25,overlap_tol=0.75):

		if not isinstance(ads,str):
			raise ValueError('ads must be a string')
		if eta not in [1,2]:
			raise ValueError('eta must be 1 or 2')
		if connect not in [1,2,3]:
			raise ValueError('connect must be 1, 2, or 3')
		if overlap_tol >= d_MX1:
			raise ValueError('The tolerance for overlapping atoms is greater than d_MX1')
		self.ads = ads
		self.d_MX1 = d_MX1
		self.r_cut = r_cut
		self.sum_tol = sum_tol
		self.rmse_tol = rmse_tol
		self.overlap_tol = overlap_tol
		self.d_X1X2 = d_X1X2
		self.d_X2X3 = d_X2X3
		self.ang_MX1X2 = ang_MX1X2
		self.ang_triads = ang_triads
		self.eta = eta
		self.connect = connect

	def get_adsorbate(self,atoms_path=None,site_idx=None,omd_path=None,
		NN_method='crystal',allowed_sites=None,write_file=True,
		new_mofs_path=None,error_path=None,NN_indices=None,atoms=None,
		new_atoms_name=None):
		"""
		Add an adsorbate using PymatgenNN or OMD

		Args:

			atoms_path (string): filepath to the CIF file
			
			site_idx (int): ASE index for the adsorption site

			omd_path (string): filepath to OMD results folder (defaults to
			'/oms_results')

			NN_method (string): string representing the desired Pymatgen
			nearest neighbor algorithm. options include 'crystal',vire','okeefe',
			and others. See NN_algos.py (defaults to 'crystal')

			allowed_sites (list of strings): list of allowed site species
			for use with automatic OMS detection

			write_file (bool): if True, the new ASE atoms object should be
			written to a CIF file (defaults to True)
			
			new_mofs_path (string): path to store the new CIF files if
			write_file is True (defaults to /new_mofs within the directory
			containing the starting CIF file)
			
			error_path (string): path to store any adsorbates flagged as
			problematic (defaults to /errors within the directory
			containing the starting CIF file)

			NN_indices (list of ints): list of indices for first coordination
			sphere (these are usually automatically detected via the default
			of None)

			atoms (ASE Atoms object): the ASE Atoms object of the MOF to add
			the adsorbate to (only include if atoms_path is not specified)

			new_atoms_name (string): the name of the MOF used for file I/O
			purposes (defaults to the basename of atoms_path if 
			provided)

		Returns:
			new_atoms (Atoms object): ASE Atoms object of MOF with adsorbate
		"""

		#Error handling of variables and setting of defaults
		if site_idx is not None and omd_path is not None:
			raise ValueError('Cannot provide site index and OMD results path')
		if atoms_path is not None and not os.path.isfile(atoms_path):
			print(f'WARNING: No MOF found for {atoms_path}')
			return None
		if atoms is not None and not isinstance(atoms,Atoms):
			raise ValueError('atoms argument must be an ASE Atoms object')
		if atoms_path is not None:
			new_atoms_name = os.path.basename(atoms_path)
			name = get_refcode(new_atoms_name)
			atoms = read(atoms_path)
		elif atoms is None:
			raise ValueError('Must specify either atoms or atoms_path args')
		if new_mofs_path is None:
			new_mofs_path = os.path.join(os.getcwd(),'new_mofs')
		if error_path is None:
			error_path = os.path.join(os.getcwd(),'errors')
		if allowed_sites is not None and not isinstance(allowed_sites,list):
			allowed_sites = [allowed_sites]
		if atoms is not None and new_atoms_name is None:
			name = atoms.get_chemical_formula()
		if site_idx is not None and site_idx < 0:
			site_idx = len(atoms)+site_idx

		self.start_atoms = atoms.copy()
		self.name = name
		self.site_idx = site_idx

		#Get ASE indices of coordinating atoms and vectors from adsorption site
		if site_idx is None:
			if omd_path is None:
				omd_path = os.path.join(os.getcwd(),'oms_results')
			omsex_dict = get_omd_data(omd_path,name,atoms)
			if omsex_dict is None:
				return None
		else:
			if NN_indices is None:
				NN_idx = get_NNs_pm(atoms,site_idx,NN_method)
			else:
				NN_idx = NN_indices
			omsex_dict = {'cnums':[len(NN_idx)],'oms_coords':[atoms[site_idx].position],
			'oms_idx':[site_idx],'oms_sym':[atoms[site_idx].symbol],'NN_coords':[atoms[NN_idx].positions],
			'NN_idx':[NN_idx]}

		new_atoms = []
		for i, oms_idx in enumerate(omsex_dict['oms_idx']):
			if allowed_sites is not None and omsex_dict['oms_sym'][i] not in allowed_sites:
				continue
			NN_idx = omsex_dict['NN_idx'][i]
			mic_coords = np.squeeze(atoms.get_distances(oms_idx,NN_idx,
				mic=True,vector=True))

			#Get the optimal adsorption site
			ads_optimizer = ads_pos_optimizer(self,new_mofs_path=new_mofs_path,
				error_path=error_path,write_file=write_file)
			ads_pos = ads_optimizer.get_opt_ads_pos(mic_coords,oms_idx)
			new_atoms_i = ads_optimizer.get_new_atoms(ads_pos,oms_idx)
			new_atoms.append(new_atoms_i)

		if len(new_atoms) == 1:
			new_atoms = new_atoms[0]

		return new_atoms

	def get_adsorbate_grid(self,atoms_path=None,site_idx=None,grid_path=None,
		grid_format='ASCII',write_file=True,new_mofs_path=None,error_path=None,
		atoms=None,new_atoms_name=None):
		"""
		This function adds a molecular adsorbate based on a potential energy grid

		Args:
			atoms_path (string): filepath to the CIF file

			site_idx (int): ASE index for the adsorption site

			grid_path (string): path to the directory containing the PEG
			(defaults to /energy_grids)

			grid_format (string): accepts either 'ASCII' or 'cube' and
			is the file format for the PEG (defaults to 'ASCII')

			write_file (bool): if True, the new ASE atoms object should be
			written to a CIF file (defaults to True)
			
			new_mofs_path (string): path to store the new CIF files if
			write_file is True (defaults to /new_mofs)
			
			error_path (string): path to store any adsorbates flagged as
			problematic (defaults to /errors)

			atoms (ASE Atoms object): the ASE Atoms object of the MOF to add
			the adsorbate to (only include if atoms_path is not specified)

			new_atoms_name (string): the name of the MOF used for file I/O
			purposes (defaults to the basename of atoms_path if 
			provided)

		Returns:
			new_atoms (Atoms object): ASE Atoms object of MOF with adsorbate
		"""

		#Error handling of variables and setting of defaults
		if atoms_path is not None and not os.path.isfile(atoms_path):
			print(f'WARNING: No MOF found for {atoms_path}')
			return None
		if atoms is not None and not isinstance(atoms,Atoms):
			raise ValueError('atoms argument must be an ASE Atoms object')
		if atoms_path is not None:
			new_atoms_name = os.path.basename(atoms_path)
			name = get_refcode(new_atoms_name)
			atoms = read(atoms_path)
		elif atoms is None:
			raise ValueError('Must specify either atoms or atoms_path args')
		if site_idx is None:
			raise ValueError('Must specify site index')
		if atoms is not None and new_atoms_name is None:
			name = atoms.get_chemical_formula()
		if site_idx < 0:
			site_idx = len(atoms)+site_idx

		grid_format = grid_format.lower()
		if grid_format == 'ascii':
			grid_ext = '.grid'
		elif grid_format == 'cube':
			grid_ext = '.cube'
		else:
			raise ValueError(f'Unsupported grid_format {grid_format}')

		if grid_path is None:
			grid_path = os.path.join(os.getcwd(),'energy_grids')
		if new_mofs_path is None:
			new_mofs_path = os.path.join(os.getcwd(),'new_mofs')
		if error_path is None:
			error_path = os.path.join(os.getcwd(),'errors')

		self.start_atoms = atoms.copy()
		self.site_idx = site_idx
		self.ads += '_grid'
		self.name = name
		max_dist = self.d_MX1

		grid_filepath = os.path.join(grid_path,name+grid_ext)

		site_pos = atoms[site_idx].position
		ads_pos = get_best_grid_pos(atoms,max_dist,site_idx,grid_filepath)
		if ads_pos == 'nogrid':
			print(f'WARNING: no grid for {name}')
			return None
		elif ads_pos == 'invalid':
			print(f'WARNING: all NaNs within cutoff for {name}')
			return None
		ads_optimizer = ads_pos_optimizer(self,new_mofs_path=new_mofs_path,
			error_path=error_path,write_file=write_file)
		return ads_optimizer.get_new_atoms_grid(site_pos,ads_pos)