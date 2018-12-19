Diatomics
===========

In this example, we'll work through how to add two-atom adsorbate to an open metal site (OMS) in a MOF. The CIF file for the MOF we'll use for this example can be found here: :download:`_static/Ni2Cl2-BBTA.cif`. This MOF is known as Ni2Cl2(BBTA) and has the structure shown below:
|Ni2Cl2-BBTA|
The metal (Ni) sites here are shown in silver. This MOF has a honeycomb-like structure with an infinite chain of metal cations going into the crystallographic c-axis, with each metal species having a five-coordinate, square pyramidal structure. Beautiful! For this example, we will consider the initialization of a single O2 molecule to a single coordinatively unsaturated Ni site.

.. |Ni2Cl2-BBTA| image:: _static/Ni2Cl2-BBTA.png
   :align: middle

We'll start with the code that can do the job. Then we'll walk through what it all means. 
.. literalinclude:: diatomic.py

Before we dive in, if you haven't already checked out the tutorial for monatomic species, please do that first so you can be brought up to speed on the general workflow of MAI. Once you've done that, it'll be easy to see how diatomics are a simple extension of monatomics!

Like with the monatomic example, we need to initialize an :clas:`mai.adsorbate_constructor` object and then provide it the MOF of interest. In the case of diatomics, we have a few new keywords to introduce. In addition to the arguments previously described in the last tutorial, we also now need to be able to tell MAI what kind of denticity we would like (i.e. end-on or side-on adsorption) and what we want the X1-X2 bond length and M-X1-X2 bond angle to be (if X1-X2 is our diatomic of interest and M is our metal adsorption site). The arguments used here are described below:

1. The :mod:`ads_species` argument is a string of the molecule that you want to add to the structure. In this example, we wanted to an O2 molecule, so we set this argument to 'O2'. Note that heteroatomic species are supported, and MAI will default to having the first atom bound to the metal site. So, setting :mod:`ads_species` to 'CO' or 'OC' would lead to M-C-O and M-O-C, respectively.
2. The :mod:`bond_dist` argument is the desired distance between the adsorption site (i.e. the Ni species) and the connecting atom of the adsorbate (in Å). Here, we set :mod:`bond_dist` to 1.5 Å.
3. The :mod:`site_idx` keyword argument is an integer representing the ASE :mod:`Atoms` index of the adsorption site (i.e. the Ni species).
4. The :mod:`eta` keyword argument is an integer representing the denticity (must be 1 or 2 and defaults to 1 if not specified). In other words, :mod:`eta=1` would be an end-on adsorption mode, whereas :mod:`eta=2` would be a side-on adsorption mode. For this example, we decided to explore both options.
5. The :mod:`d_bond` keyword argument is the desired distance between X1 and X2 in the diatomic (in Å). If not specified, it will default to the value for :mod:`bond_dist`. Here, we decided to set :mod:`d_bond` to 1.2 Å, which is a reasonable O-O bond distance.
6. The :mod:`angle` keyword argument is the desired M-X1-X2 bond angle (in degrees). By default, it assumes :mod:`angle=180` if :mod:`eta=1` or :mod:`angle=90` if :mod:`eta=2`. For this example, we use :mod:`angle=120` and :mod:`angle=90`, respectively, which is representative of common O2 binding modes.

That takes care of initializing :clas:`mai.adsorbate_constructor` object. Now we can use this object to call a function ot make initialize the adsorbate with the file path to the MOF's CIF file.

Now let's see what happens as a result of running this code! The initialized structure is shown below:
|Ni2Cl2-BBTA-O2|
Exactly what we'd expect yet again! You can see that in the first example, O2 is bound end-on, whereas the latter it is bound side-on, as specified in the example script. The bond angles and distances are the same as those specified in the input file. As mentioned in the previous tutorial, MAI attempts to minimize sterics. When it comes to diatomics, the adsorbate is spun around the M-X1 axis so that it minimizes interactions with framework atoms.

That concludes our tutorial with diatomic adsorbates. Join me as we move onto more complicated systems! Up next is triatomics!

.. |Ni2Cl2-BBTA-O2| image:: _static/Ni2Cl2-BBTA-O2.png
   :align: middle