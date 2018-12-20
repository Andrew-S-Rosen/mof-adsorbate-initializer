Triatomics
===========

In this example, we'll work through how to add three-atom adsorbates to open metal sites in MOFs. The CIF for the MOF we'll use for this example can be found :download:`here <_static/Ni2Cl2-BTDD.cif>`. This MOF is known Ni2Cl2(BTDD) and has the structure shown below:

|Ni2Cl2-BTDD|

The metal (Ni) sites here are shown in silver. This MOF has a honeycomb-like structure with square pyramidal Ni cations that run down the crystallographic c-axis. Isn't it beautiful?

.. |Ni2Cl2-BTDD| image:: _static/Ni2Cl2-BTDD.png
   :align: middle

--------------------
Contiguous Adsorbate
--------------------
For this example, we will consider the initialization of a "contiguous" triatomic adsorbate. What I mean by this is that if we treat the adsorbate as an arbitrary molecule X1-X2-X3, then the adsorption process is described by M-X1-X2-X3. In this example, we'll consider the adsorption of N2O to a Ni site, in both an η1-N and η1-O binding mode. The code to handle this is shown below. 

.. literalinclude:: _static/triatomic_linear.py

Like with the previous examples, we need to initialize an :class:`~mai.adsorbate_constructor.adsorbate_constructor` object and then provide it the MOF of interest. In the case of triatomics, we have a few new keywords to introduce. In addition to the arguments described in the previous examples, we can now provide MAI additional geometric parameters if desired. Namely, the new features are now that we can include the X2-X3 bond length and the X1-X2-X3 bond angle. The arguments used here are described below:

1. The ``ads_species``, ``bond_dist``, ``site_idx``, ``d_bond``, and ``angle`` are the same as before.
2. Now, we have the option to add the ``d_bond2`` keyword argument, which specifies the X1-X2 distance. It defaults to ``d_bond2=d_bond1`` if not specified.
3. We can also add the ``angle2`` keyword argument, which specifies the X1-X2-X3 bond angle. It defaults to ``angle2=180`` if not specified. We'll leave this one at the default for this example.

That takes care of initializing the :class:`~mai.adsorbate_constructor.adsorbate_constructor` object. With this, we provide the object with the path to the MOF, and it will initialize the adsorbate. Now let's see what happens as a result of running this code! The initialized structures are shown below:

|Ni2Cl2-BTDD-N2O|

Exactly what we'd expect once more! You can see that in the first example, N2O is bound in an η1-N mode, whereas the second is bound in an η1-O mode, as specified. Feel free to play around with the bond distance and bond angle arguments to get a feel for how MAI works.

.. |Ni2Cl2-BTDD-N2O| image:: _static/Ni2Cl2-BTDD-N2O.png
   :align: middle

---------------------
Noncontiguous Adsorbate
---------------------
The last bit of trickery comes into play when dealing with what I'll call "noncontiguous" adsorbates. These are adsorbates like water, where it is triatomic, but it is not bound in a sequential fashion. As with water, you will have a central atom of the adsorbate bound to the metal. Therefore, ``connect`` must be set to ``connect=2``, and the name of the species must be specified with the connecting atom listed as the second atom in the ``ads_species`` name (e.g. ``ads_species=HOH`` for water bound in an η1-O mode if ``connect=2``). An example code is shown below. The main thing to keep in mind is that now the connecting atom of the adsorbate is X2 instead of X1.

.. literalinclude:: _static/triatomic_central.py

.. |Ni2Cl2-BTDD-CO| image:: _static/Ni2Cl2-BTDD-H2O.png
   :align: middle

That concludes our tutorial for triatomic adsorbates.