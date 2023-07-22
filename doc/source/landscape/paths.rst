(Re-)folding Paths, Saddle Points, Energy Barriers, and Local Minima
====================================================================

API for various RNA folding path algorithms.

This part of our API allows for generating RNA secondary structure (re-)folding
paths between two secondary structures or simply starting from a single
structure.

This is most important if an estimate of the refolding energy barrier between
two structures is required, or a structure's corresponding local minimum needs
to be determined, e.g. through a gradient-descent walk.

This part of the interface is further split into the following sections:

- :doc:`/landscape/paths_direct`, and
- :doc:`/landscape/paths_walk`


.. doxygengroup:: paths
    :no-title:
