Potential Energy Grids
======================
Some adsorbates do not lend themselves well to the geometric approaches laid out thus far. This is particularly the case for adsorbates that are physisorbed near an adsorption site, as opposed to chemisorbed. In these cases, an alternative way of initializing the adsorbate is by mapping a out a potential energy grid of each MOF and putting the adsorbate in a low-energy site within some cutoff radius of the proposed adsorption site.

MAI supports two different formats for potential energy grids. The first is a space-delimited file with four columns of (x,y,z,E), where E is the potential energy and (x,y,z) are the coordinates. Each new line represents a new set of (x,y,z,E) data. This is the file-format produced from RASPA_-generated potential energy grids. MAI also supports the cube_ file format for grids, such as those generated from `PorousMaterials.jl <https://github.com/SimonEnsemble/PorousMaterials.jl>`_. Currently, only single-site CH4 adsorbates are supported, although in principle it should be trivial to consider other adsorbates as well. When use potential energy grids to initialize the position of CH4 adsorbates, the C atom of the CH4 molecule will be placed in the low-energy site, and the four remaining H atoms will be arranged to form the tetrahedral structure of CH4, with one of the H atoms pointed directly toward the adsorption site.

>details go here. get_adsobrate_grid etc.

.. _RASPA: https://www.tandfonline.com/doi/full/10.1080/08927022.2015.1010082
.. _cube: http://paulbourke.net/dataformats/cube/
