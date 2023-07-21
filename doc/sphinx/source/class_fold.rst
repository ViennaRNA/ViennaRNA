Classified Dynamic Programming Variants
=======================================

Usually, thermodynamic properties using the basic recursions for :doc:`/mfe`,
:doc:`pf_fold`, and so forth, are computed over the entire structure space. However,
sometimes it is desired to partition the structure space *a priori* and compute
the above properties for each of the resulting partitions. This approach directly leads
to *Classified Dynamic Programming*.

.. toctree::
   :maxdepth: 1
   :caption: Specialized Modules:

   classified/kl_neighborhood
   classified/dos

.. doxygengroup:: class_fold
