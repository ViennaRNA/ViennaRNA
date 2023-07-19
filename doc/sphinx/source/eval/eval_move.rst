Energy Evaluation for Atomic Moves
==================================


Functions to evaluate the free energy change of a structure after application of
(a set of) atomic moves

Here, atomic moves are not to be confused with moves of actual physical atoms. Instead,
an atomic move is considered the smallest conformational change a secondary structure
can undergo to form another, distinguishable structure. We currently support the
following moves

- Opening (dissociation) of a single base pair
- Closing (formation) of a single base pair
- Shifting one pairing partner of an existing pair to a different location

.. toctree::
   :maxdepth: 1
   :caption: Contents:

.. doxygengroup:: eval_move
    :no-title:


