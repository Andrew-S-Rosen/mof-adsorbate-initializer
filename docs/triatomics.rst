Triatomics
===========

In this example, we'll work through how to add three-atom adsorbates to an open metal site in a MOF. The CIF file for the MOF we'll use for this example can be found here: :download:`Ni2Cl2-BBTA.cif`. This MOF is known as Ni2Cl2-BBTA and has the structure shown below:
|Ni2Cl2-BBTA|
The metal (Sc) sites here are shown in purple. This MOF has a trimetallic node, with two coordinatively unsaturated metal sites per node. For this example, we will consider the initialization of a N2O molecule to a single coordinatively unsaturated Sc site.

.. |Ni2Cl2-BBTA| image:: _static/Ni2Cl2-BBTA.png
   :align: middle

We'll start with the code that can do the job. Then we'll walk through what it all means. 
.. literalinclude:: _static/triatomic.py

Like with the previous examples, we need to initialize an :class:`mai.adsorbate_constructor` object and then provide it the MOF of interest. In the case of triatomics, we have a few new keywords to introduce. In addition to the arguments previously described in the previous tutorials, we also now need to be able to tell MAI which connecting atom we would like. This was simple for diatomics, but for triatomics we can 
 we want the X2-X3 bond length to be and which connecting what we want the X1-X2-X3 bond angle to be kind of denticity we would like (i.e. end-on or side-on adsorption) and what we want the X1-X2 bond length and M-X1-X2 bond angle to be (if X1-X2 is our diatomic of interest and M is our metal adsorption site). The arguments used here are described below:

1. The ``ads_species`` argument is a string of the molecule that you want to add to the structure. In this example, we wanted to an O2 molecule, so we set this argument to 'O2'. Note that heteroatomic species are supported, and MAI will default to having the first atom bound to the metal site. So, setting ``ads_species`` to 'CO' or 'OC' would lead to M-C-O and M-O-C, respectively.
2. The ``bond_dist`` argument is the desired distance between the adsorption site (i.e. the Ni species) and the connecting atom of the adsorbate (in Å). Here, we set ``bond_dist`` to 1.5 Å.
3. The ``site_idx`` keyword argument is an integer representing the ASE ``Atoms`` index of the adsorption site (i.e. the Ni species).
4. The ``eta`` keyword argument is an integer representing the denticity (must be 1 or 2 and defaults to 1 if not specified). In other words, ``eta=1`` would be an end-on adsorption mode, whereas ``eta=2`` would be a side-on adsorption mode. For this example, we decided to explore both options.
5. The ``d_bond`` keyword argument is the desired distance between X1 and X2 in the diatomic (in Å). If not specified, it will default to the value for ``bond_dist``. Here, we decided to set ``d_bond`` to 1.2 Å, which is a reasonable O-O bond distance.
6. The ``angle`` keyword argument is the desired M-X1-X2 bond angle (in degrees). By default, it assumes ``angle=180`` if ``eta=1`` or ``angle=90`` if ``eta=2``. For this example, we use ``angle=120`` and ``angle=90``, respectively, which is representative of common O2 binding modes.

That takes care of initializing :class:`mai.adsorbate_constructor` object. Now we can use this object to call a function ot make initialize the adsorbate with the file path to the MOF's CIF file.

Now let's see what happens as a result of running this code! The initialized structure is shown below:
|Ni2Cl2-BBTA-N2O|
Exactly what we'd expect yet again! You can see that in the first example, O2 is bound end-on, whereas the latter it is bound side-on, as specified in the example script. The bond angles and distances are the same as those specified in the input file. As mentioned in the previous tutorial, MAI attempts to minimize sterics. When it comes to diatomics, the adsorbate is spun around the M-X1 axis so that it minimizes interactions with framework atoms.

That concludes our tutorial with diatomic adsorbates. Join me as we move onto more complicated systems! Up next is triatomics!

.. |Ni2Cl2-BBTA-N2O| image:: _static/Ni2Cl2-BBTA-N2O.png
   :align: middle