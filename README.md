# MOF Adsorbate Initializer (MAI)
Python code to initialize the position of adsorbates in MOFs for high-throughput DFT screening workflows. Relevant details for the code can be found in the following manuscript (once published):

A.S. Rosen, J.M. Notestein, R.Q. Snurr. "Identifying Promising Metal-Organic Frameworks for Heterogeneous Catalysis via High-Throughput Periodic Density Functional Theory." Under review. 

## What is MAI?
High-throughput computational catalysis typically involves the calculation of many different adsorption energies. MAI is a tool meant to automate the process of constructing reasonable initial guesses for adsorbates in MOFs for the purposes of DFT screening workflows. Currently, MAI supports two main types of calculations, which are described below.

### Adding Adsorbates via Molecular Symmetry
Open metal sites (OMSs), or coordinatively unsaturated metal sites, are commonly studied for catalytic applications. MAI makes it possible to systematically initialize the positions of adsorbates to OMSs in MOFs for further DFT analysis by orienting the adsorbate in such a way that the symmetry of the new coordination environment is maximized.

The simplest method is to specify the atom index where the adsorbate should be added (more complex methods are discussed below). The following is minimal example using MAI:

```python
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'ANUGIA.cif' #path to the CIF file
site_idx = 0 #atom index of adsorption site
adsorbate = 'O' #adsorbate species
r = 2.0 #desired bond length

#add adsorbate
ads = adsorbate_constructor(adsorbate,r,site_idx=site_idx)
new_mof, new_mof_name = ads.get_adsorbate_pm(mof_path)
```

Here, we have used MAI's `adsorbate_constructor` class to tell the program what kind of adsorbate we'd like and where we want it (i.e. an 'O' atom 2.0 Å away from `site_idx=0` in `ANUGIA.cif`). The `adsorbate_constructor` class has a variety of functions to add adsorbates. The one in this example is `get_adsorbate_pm`, which uses the [local environment features](https://www.frontiersin.org/articles/10.3389/fmats.2017.00034/full) of Pymatgen ("pm") to add an adsorbate based on the detected bonding topology of the metal center.

### Adding Adsorbates via a Potential Energy Grid
Occasionally, symmetry arguments may not be sufficient for initializing the position of an adsorbate. For instance, if we want to know where to consider putting a molecule (like CH4) around the vacinity of an active site in a MOF, this may not be immediately obvious, especially for a wide range of MOF structures. The use molecular mechanics-based potential energy grids (PEGs) can be used to provide an inexpensive, robust, and systematic initial guess based on local minima in the PEG.

The minimal example for this scenario is fundamentally similar to that in the prior section, except we replace `get_adsorbate_pm` with `get_adsorbate_grid`. The `adsorbate_constructor` class is therefore very flexible and allows for a variety of methods to initialize adsorbate positions. Details for the various functions within the `adsorbate_constructor` class are described in the following sections. 

## Ready-to-Run Examples
To get started, sample scripts are provided in the `examples` base directory of MAI. They are summarized below:
1. `add_O.py`. Adds an O adsorbate to `ANUGIA.cif` given a pre-determined Cu site.
2. `add_O_OMS_zeo.py`. Adds an O adsorbate to `ANUGIA.cif` using the OMS detection algorithm in Zeo++.
3. `add_O_OMS_omd.py`. Adds an O adsorbate to `ANUGIA.cif` using the OMS detectiong algorithm in [Open_Metal_Detector.py](https://github.com/emmhald/open_metal_detector).
4. `add_CH4_PEG_ASCII.py`. Adds a CH4 adsorbate to a terminal Cu-oxo site in `AHOKIR01-O.cif` using an ASCII-based PEG.
5. `add_CH4_PEG_cube.py`. Adds a CH4 adsorbate to a terminal Cu-oxo site in `ANUGIA_clean_min_O.cif` using a cube-based PEG (`todo`).

## `adsorbate_constructor`
As previously discussed, the main tool to initialize adsorbate positions is the `adsorbate_constructor`, which is described below:

```python
class adsorbate_constructor():
	"""
	This class constructs an ASE atoms object with an adsorbate
	"""
	def __init__(self,ads_species,bond_dist,site_species=None,site_idx=None,
		r_cut=2.5,sum_tol=0.5,rmse_tol=0.25,overlap_tol=0.75):
		"""
		Initialized variables

		Args:
			ads_species (string): string of atomic element for adsorbate (e.g. 'O')
			
			bond_dist (float): distance between adsorbate and surface atom. If
			used with get_adsorbate_raspa, it represents the maximum distance
			for the adsorbate from the surface atom
			
			site_species (string): string of atomic element for the adsorption
			site species. not needed when calling get_adsorbate_zeo_oms
			
			site_idx (int): ASE index for the adsorption site (defaults to
			the last element of element type site_species). not needed when
			calling get_adsorbate_zeo_oms
			
			r_cut (float): cutoff distance for calculating nearby atoms when
			ranking different adsorption sites
			
			sum_tol (float): threshold to determine planarity. when the sum
			of the Euclidean distance vectors of coordinating atoms is less
			than sum_tol, planarity is assumed
			
			rmse_tol (float): second threshold to determine planarity. when the 
			root mean square error of the best-fit plane is less than rmse_tol,
			planarity is assumed
			
			overlap_tol (float): distance below which atoms are assumed to be
			overlapping so that the structure can be flagged as erroneous
		"""
```
Once the `adsorbate_constructor` is instanced, one of three routines can be called: `get_adsorbate_pm`, `get_adsorbate_oms`,`get_adsorbate_grid`. These are described below.

## `get_adsorbate_pm`
This function uses Pymatgen's [`local_env`](http://pymatgen.org/_modules/pymatgen/analysis/local_env.html) features to determind the coordination environment of a given active site and then adds an adsorbate based on molecular symmetry and the degree of planarity. By combining MAI with other toolkits such as [ASE](https://gitlab.com/ase/ase) or [Pymatgen](http://pymatgen.org), it is possible to carry out significantly more complex tasks by automatically determining the desired `site_idx` based on your particular needs.

```python
def get_adsorbate_pm(self,atoms_filepath,NN_method='vire',write_file=True,
	new_mofs_path=None,error_path=None):
	"""
	Use Pymatgen's nearest neighbors algorithms to add an adsorbate

	Args:
		atoms_filepath (string): filepath to the structure file (accepts
		CIFs, POSCARs, and CONTCARs)
		NN_method (string): string representing the desired Pymatgen
		nearest neighbor algorithm (options include 'vire','voronoi',
		'jmol','min_dist','okeeffe','brunner_relative','brunner_reciprocal',
		'brunner_real', and 'econ')
		write_file (bool): if True, the new ASE atoms object should be
		written to a CIF file (defaults to True)
		new_mofs_path (string): path to store the new CIF files if
		write_file is True (defaults to /new_mofs within the
		directory of the starting structure files)
		error_path (string): path to store any adsorbates flagged as
		problematic (defaults to /errors within the directory of the
		starting structure files)
	Returns:
		new_atoms (Atoms object): ASE Atoms object of MOF with adsorbate
		new_name (string): name of MOF with adsorbate
	"""
```

## `get_adsorbate_oms`
Since OMSs are widely investigated for catalysis in MOFs, `get_adsorbate_oms` will automatically detect chemically unique OMSs in a given MOF and will add an adsorbate accordingly. This means that neither `site_idx` nor `site_species` need to be specified in `adsorbate_constructor`. The details for `get_adsorbate_oms` are described below:

```python
def get_adsorbate_oms(self,atoms_filepath,oms_data_path=None,
		write_file=True,new_mofs_path=None,error_path=None):
		"""
		This function adds an adsorbate to each unique OMS in a given
		structure. In cases of multiple identical OMS, the adsorbate with
		fewest nearest neighbors is selected. In cases of the same number
		of nearest neighbors, the adsorbate with the largest minimum distance
		to extraframework atoms (excluding the adsorption site) is selected.

		Args:
			atoms_filepath (string): filepath to the structure file (accepts
			CIFs, POSCARs, and CONTCARs)
			oms_data_path (string): path to the Zeo++ open metal site data
			containing .oms and .omsex files (defaults to /oms_data within the
			directory of the starting structure files)
			write_file (bool): if True, the new ASE atoms object should be
			written to a CIF file (defaults to True)
			new_mofs_path (string): path to store the new CIF files if
			write_file is True (defaults to /new_mofs within the
			directory of the starting structure files)
			error_path (string): path to store any adsorbates flagged as
			problematic (defaults to /errors within the directory of the
			starting structure files)
		Returns:
			new_atoms_list (list): list of ASE Atoms objects with an adsorbate
			added to each unique OMS
			new_name_list (list): list of names associated with each atoms
			object in new_atoms_list
		"""
```
To effectively use `get_adsorbate_oms`, one of two programs must be used to generate the OMS data needed for MAI to properly work: [OpenMetalDetector (OMD)](https://github.com/emmhald/open_metal_detector) or [Zeo++](http://www.zeoplusplus.org/about.html). For ease-of-use, we recommend OMD. A minimal example using MAI with OMS is shown below:

```python
from omsdetector import MofCollection
from mai.adsorbate_constructor import adsorbate_constructor

mof_path = 'MyMOFFolder' #path to CIFs
analysis_path = mof_path #path to store the OMS results
mof_coll = MofCollection.from_folder(collection_folder=mof_path,analysis_folder=analysis_path) #import MOFs
ads_species = 'O' #adsorbate species
bond_length = 2.0 #desired distance between OMS and ads_species

#add adsorbate
ads = adsorbate_constructor(ads_species,bond_length)
new_mof_atoms, new_mof_name = ads.get_adsorbate_oms(mof_path,oms_format='OMD')
```

## `get_adsorbate_grid`
The `get_adsorbate_grid` function is used to add an adsorbate to the adsorption site of a MOF based on a molecular mechanics PEG, such as that comptued by [RASPA](https://www.tandfonline.com/doi/abs/10.1080/08927022.2015.1010082) or [PorousMaterials.jl](https://github.com/SimonEnsemble/PorousMaterials.jl). The molecular adsorbate is initialized in the lowest energy position within some cutoff distance from the proposed adsorption site. It currently only supports CH4 adsorption but can be readily extended to support other molecular adsorbates. The `get_adsorbate_grid` function is described below:

```python
def get_adsorbate_grid(self,atoms_filepath,grid_path=None,
	write_file=True,new_mofs_path=None,error_path=None):
	"""
	This function adds a molecular adsorbate based on an energy grid
	generated using RASPA

	Args:
		atoms_filepath (string): filepath to the structure file (accepts
		CIFs, POSCARs, and CONTCARs)
		
		grid_path (string): path to the directory containing RASPA energy
		grids (defaults to /energy_grids within the directory of the
		starting structure files)
		
		write_file (bool): if True, the new ASE atoms object should be
		written to a CIF file (defaults to True)
		
		new_mofs_path (string): path to store the new CIF files if
		write_file is True (defaults to /new_mofs within the
		directory of the starting structure files)
		
		error_path (string): path to store any adsorbates flagged as
		problematic (defaults to /errors within the directory of the
		starting structure files)
	Returns:
		new_atoms (Atoms object): ASE Atoms object of MOF with adsorbate
		
		new_name (string): name of MOF with adsorbate
	"""
```

`get_adsorbate_grid` requires a PEG in one of two formats: ASCII or .cube.

The ASCII format required by MAI is that generated by RASPA. Specifically, for each point in the PEG, it must contain three columns consisting of the x, y, and z direct coordinates and the computed energy (any values of `?` for the energy column and any additional columns are ignored). An example is provided in `examples/example_MOFs/energy_grids/AHOKIR01-O.grid` for reference. 

`Todo`: Support for the .cube format is under development.

## Setup

### Installing MAI
1. MAI requires [Python](https://www.python.org/) 3.6 or newer. If you do not already have Python installed, the easiest option is to download the [Anaconda](https://www.anaconda.com/download/) distribution.
2. Download or clone the MAI repository and run `pip install .` from the MAI base directory.

### Required Dependencies
MAI requires the following Python packages:
1. [ASE](https://wiki.fysik.dtu.dk/ase/) 3.16.0 or newer
2. [Pymatgen](http://pymatgen.org/) 2018.7.15 or newer.

Both packages can be installed via `pip install ase pymatgen`. 

### External Programs

#### Detecting Open Metal Sites
MAI relies on one of two programs to automatically identify OMSs via `add_adsorbate_oms`.
1. [OpenMetalDetector (OMD)](https://github.com/emmhald/open_metal_detector). To install, download or clone the OMD repository and run `pip install .` frmo the OMD base directory. Refer to the `README` on the OMD GitHub page for details of how to use the program.
2. [Zeo++](http://www.zeoplusplus.org/) version 0.3. To install, just download the source folder to your desired directory. By default, Zeo++'s OMS detection algorithm does not output the positions of each OMS and its coordinating atoms, which are needed for MAI. To address this, copy `network.cc` from `patches/network.cc` in the MAI base directory and replace the corresponding `network.cc` file in the base directory of Zeo++ before installation. The relevant changes can be found starting on line 1165. To use Zeo++'s' OMS detection algorithm, run `/zeo++-0.3/network -omsex filepath`, where `-omsex` requests OMS detection with extended output, and `filepath` is the path to the CIF file of the MOF. This will produce the required `.oms` and `.omsex` files for each MOF.

#### Potential Energy Grids
MAI relies on one of two programs to add adsorbates based on computed PEGs via `add_adsorbate_grid`.
1. [RASPA](https://www.tandfonline.com/doi/full/10.1080/08927022.2015.1010082).
2. [PorousMaterials.jl](https://github.com/SimonEnsemble/PorousMaterials.jl).