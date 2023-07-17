Energy Parameters
=================

For secondary structure free energy evaluation we usually utilize the set of
thermodynamic **Nearest Neighbor** energy parameters also used in other software,
such as *UNAfold* and *RNAstructure*.

.. toctree::
   :maxdepth: 1
   :caption: Specialized Modules:

   params/salt
   params/io
   params/convert

Available Parameter Sets
------------------------

While the *RNAlib* already contains a compiled-in set of the latest
*Turner 2004 Free Energy Parameters*, we defined a file format that allows to
change these parameters at runtime. The ViennaRNA Package already comes with
a set of parameter files containing

- Turner 1999 RNA parameters
- Mathews 1999 DNA parameters
- Andronescu 2007 RNA parameters
- Mathews 2004 DNA parameters

Energy Parameter API
--------------------

.. doxygengroup:: energy_parameters

