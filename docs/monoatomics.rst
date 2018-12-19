Monoatomics
===========

In this example, we'll work through how to add a single atom adsorbate to an open metal site (OMS) in a MOF. The CIF file for the MOF we'll use for this example can be found here: :download:`_static/Cu-BTC.cif`. This MOF is known as Cu-BTC and has the structure shown below:
|Cu-BTC|

The metal (Cu) sites here are shown in orange. As you can see, there are multiple Cu sites per unit cell, and each Cu site is in a paddlewheel-like structure. For this example, we consider the initialization of an O atom adsorbate to a single coordinatively unsaturated Cu site. 

.. |Cu-BTC| image:: _static/Cu-BTC.png
   :align: middle

We'll start with the code that can do the job. Then we'll walk through what it all means.

.. literalinclude:: _static/monoatomic.py

Okay, let's dive right in! MAI requires the calling of a module known as the :class:`mai.adsorbate_constructor` object, which tells MAI what kind of adsorbate you'd like to make. For simple monoatomic species, there are only a few arguments you need to worry about. 

1. The :mod:`ads_species` argument is a string of the element (or molecule, for larger adsorbates) that you want to add to the structure. In this example, we wanted to an oxygen atom, so we set this argument to 'O'. Note that :mod:`ads_species` is case-sensitive.
2. The :mod:`bond_dist` argument is the desired distance between the adsorption site (i.e. the Cu species) and the adsorbate (in Å). Here, we set :mod:`bond_dist` to 1.75 Å.
3. The :mod:`site_idx` keyword argument is an integer representing the ASE :mod:`Atoms` index of the adsorption site (i.e. the Cu species). Later in this tutorial, we'll show how this parameter can be determined automatically, but for now we have manually set it to the 0th index, which corresponds to one of the Cu atoms. To find out the ASE indices, you can read in the CIF file via :mod:`from ase.io import read; read('Cu-BTC.cif')` and inspect_ or visualize_ the :mod:`Atoms` object. Generally, they are the same order indices as you'd find in your favorite CIF viewer (e.g. VESTA_).

That takes care of initializing :class:`mai.adsorbate_constructor` object. Now we can use this object to call a function ot make initialize the adsorbate. This is done via :mod:`adsorbate_constructor.get_adsorbate`. Generally, the only information you'll need to provide is the file path to the CIF file you want to add the adsorbate to. In this case, we set that to just :mod:`Cu-BTC.cif` in our current working directory. The output of calling :mod:`get_adsorbate` is a new ASE :mod:`Atoms` object with the adsorbate initialized.

Now let's see what happens as a result of running this code! The initialized structure is shown below:
|Cu-BTC-O|
Exactly what we'd expect! Note that, at a high level, MAI aims to satisfy two major conditions. The first is that it tries to maximize between the first coordination sphere and the adsorbate. In this case, the geometry is square planar prior to adsorption, so MAI makes a square pyramidal structure. The second condition is that MAI tries to minimize steric interactions when possible. In the case of a paddlewheel structure like Cu-BTC, the monatomic adsorbate could have been initialize in either the + or - direction normal to the planar first coordination sphere. However, only one of those directions is geometrically accessible (the other is pointed inward inbetween the paddlewheel, which would not be reasonable).

That concludes our tutorial with monatomic adsorbates. Join me as we move onto more complicated systems! Up next is diatomics!

.. |Cu-BTC-O| image:: _static/Cu-BTC-O.png
   :align: middle
.. _inspect: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
.. _visualize: https://wiki.fysik.dtu.dk/ase/ase/visualize/visualize.html
.. _VESTA: https://jp-minerals.org/vesta/en/