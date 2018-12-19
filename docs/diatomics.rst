Diatomics
===========

In this example, we'll work through how to add two-atom adsorbates to an open metal site in a MOF. The CIF file for the MOF we'll use for this example can be found here: :download:`Ni2Cl2-BBTA.cif`. This MOF is known Ni2Cl2-BBTA and has the structure shown below:
|Ni2Cl2-BBTA|
The metal (Ni) sites here are shown in silver. This MOF has a honeycomb-like structure with an infinite chain of metals running down the crystallographic c-axis, each one in a square pyramidal geometry. Neat!

.. |Ni2Cl2-BBTA| image:: _static/Ni2Cl2-BBTA.png
   :align: middle

-----------
Homoatomic
-----------
For this example, we will consider the initialization of an O2 molecule to a single coordinatively unsaturated Ni site. O2 can bind in an end-on or side-on mode depending on the structure. We'll consider both for this example. The code to handle this is shown below. 
.. literalinclude:: diatomic_homo.py

Like with the monatomic example, we need to initialize an :class:`mai.adsorbate_constructor` object and then provide it the MOF of interest. In the case of diatomics, we have a few new keywords to introduce. In addition to the arguments described in the monatomic tutorial, we now need to be able to tell MAI what kind of denticity we would like (i.e. end-on or side-on adsorption) and what we want the X1-X2 bond length and M-X1-X2 bond angle to be (if X1-X2 is our diatomic of interest and M is our metal adsorption site). The arguments used here are described below:

1. The :mod:`ads_species` argument is a string of the molecule that you want to add to the structure. In this example, we wanted to add an O2 molecule, so we set this argument to 'O2'. Note that MAI will internally strip any characters following (and including) an underscore, so 'O2_end' and 'O2_side' both get stripped to 'O2'. That being said, the full string for the :mod:`ads_species` argument will be used when writing the filenames of the new CIF files, so using an underscore can be helpful for organizational purposes.
2. The :mod:`bond_dist` argument is the desired distance between the adsorption site (i.e. the Ni species) and the adsorbate (in Å). If the adsorbate is bound in an end-on fashion, this represents the M-X1 distance. If the adsorbate is bound in a side-on fashion, this represents the distance between M and the midpoint between X1 and X2. Here, we set :mod:`bond_dist` to 1.5 Å.
3. The :mod:`site_idx` keyword argument is an integer representing the ASE :mod:`Atoms` index of the adsorption site (i.e. the Ni species).
4. The :mod:`eta` keyword argument is an integer representing the denticity. In other words, :mod:`eta=1` would be an end-on adsorption mode, whereas :mod:`eta=2` would be a side-on adsorption mode. By default, :mod:`eta=1` if unspecified. For this example, we decided to explore both options.
5. The :mod:`d_bond` keyword argument is the desired distance between X1 and X2 (in Å). If not specified, it will default to the value for :mod:`bond_dist`. Here, we decided to set :mod:`d_bond` to 1.2 Å, which is a reasonable O-O bond distance.
6. The :mod:`angle` keyword argument is the angle between the adsorption site and the adsorbate (in degrees). If the adsorbate is bound in an end-on fashion, this represents the M-X1-X2 bond angle. If the adsorbate is bound in a side-on fashion, this represents the angle between M, the midpoint between X1 and X2, and X2. By default, it assumes :mod:`angle=180` if :mod:`eta=1` or :mod:`angle=90` if :mod:`eta=2`. For this example, we use :mod:`angle=120` and :mod:`angle=90`, respectively, which is representative of common O2 binding geometries.

That takes care of initializing the :class:`mai.adsorbate_constructor` object. With this, we provide the object with the path to the MOF, and it will initialize the adsorbate. Now let's see what happens as a result of running this code! The initialized structure is shown below:
|Ni2Cl2-BBTA-O2|
Exactly what we'd expect yet again! You can see that in the first example, O2 is bound end-on, whereas in the second it is bound side-on, as specified in the example script. The bond angles and distances are the same as those specified in the input file. As mentioned in the previous tutorial, MAI attempts to maximize symmetry and minimize sterics. For the symmetry in this case, the originally square pyrmidal geomtry becomes octahedral. When it comes to handling sterics for diatomics, the adsorbate it is spun around at the user-specified bond distance and angle so that it minimizes interactions with the framework atoms.

.. |Ni2Cl2-BBTA-O2| image:: _static/Ni2Cl2-BBTA-O2.png
   :align: middle

------------
Heteroatomic
------------
MAI also supports heteratomic adsorbates. In this example, we'll consider the adsorption of a single CO molecule with the same Ni2Cl2-BBTA MOF. The only thing that changes for heteroatomic adsorbates is that you need to tell MAI which atom is the "connecting atom" (i.e. the atom of the adsorbate bound to the metal adsorption site) if bound in an end-on fashion. By default, MAI will assume that the first atom in :mod:`ads_species` is the connecting atom. Therefore, setting :mod:`ads_species='CO'` or :mod:`ads_species='OC'` would yield M-C-O or M-O-C binding modes, respectively. There is also a :mod:`connect` keyword argument in :class:`mai.adsorbate_constructor`, which is an integer representing the atom in :mod:`ads_species` that should be the connecting atom. For :mod:`ads_species='Co'`, using :mod:`connect=1` or :mod:`connect=2` would once again refer to M-C-O or M-O-C binding modes, respectively.
.. literalinclude:: diatomic_hetero.py 

.. |Ni2Cl2-BBTA-CO| image:: _static/Ni2Cl2-BBTA-CO.png
   :align: middle

That concludes our tutorial for diatomic adsorbates. Now onto triatomics!