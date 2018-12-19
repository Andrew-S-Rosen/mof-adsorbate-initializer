============
Installation
============

Dependencies
============

1. Python_ 3.6 or newer

2. Pymatgen_ 2018.11.30 or newer

3. ASE_ 3.16.0 or newer

Optional:

4. OpenMetalDetector_

.. _Python: http://www.python.org/
.. _Pymatgen: http://pymatgen.org/
.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _OpenMetalDetector: https://github.com/emmhald/open_metal_detector


Installation Instructions
=========================
1. If you don't already have Python installed on your machine, you'll need to install Python 3.6 or newer. I recommend using the Anaconda_ distribution of Python 3 if you don't already have it installed.

2. You will then need to install MAI and all of the required dependencies. The easiest way to do this is to find the :mod:`requirements.txt` file in the MAI base directory. If you have an Anaconda distribution of Python, you can then use :mod:`while read requirement; do conda install --yes $requirement; done < requirements.txt` to install the necessary packages. Otherwise, you can use the Python :mod:`pip` package manager with the command :mod:`pip install -r requirements.txt`.

.. _Anaconda: https://www.anaconda.com/download/