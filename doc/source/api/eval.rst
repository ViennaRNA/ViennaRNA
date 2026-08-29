Free Energy Evaluation
======================

.. _api-energy-parameters:

Energy Parameters
-----------------

Functions related to the set of thermodynamic **Nearest Neighbor** energy parameters.

.. doxygengroup:: energy_parameters
    :no-title:

.. _api-energy-parameters-salt:

Salt Corrections
^^^^^^^^^^^^^^^^

All relevant functions to compute salt correction at a given salt concentration and temperature.

.. doxygengroup:: energy_parameters_salt
    :no-title:

Loading / Saving Energy Parameter Sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read and Write energy parameter sets from and to files or strings

.. doxygengroup:: energy_parameters_rw
    :no-title:

.. _api-energy-parameters-convert:

Converting Energy Parameter Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Converting energy parameter files into the latest format.

.. doxygengroup:: energy_parameters_convert
    :no-title:


.. _api-model-details:

Fine-tuning of the Implemented Models
-------------------------------------

Functions and data structures to fine-tune the implemented
secondary structure evaluation model.

.. doxygengroup:: model_details
    :no-title:


.. _api-eval-loops:

Energy Evaluation for Individual Loops
--------------------------------------

General
^^^^^^^

.. doxygengroup:: eval_loops
    :no-title:


Exterior Loops
^^^^^^^^^^^^^^

Functions to evaluate the free energy contributions for exterior (external) loops.

.. doxygengroup:: eval_loops_ext
    :no-title:


Hairpin Loops
^^^^^^^^^^^^^^

Functions to evaluate the free energy contributions for hairpin loops.

.. doxygengroup:: eval_loops_hp
    :no-title:


Internal Loops
^^^^^^^^^^^^^^

Functions to evaluate the free energy contributions for internal loops.

.. doxygengroup:: eval_loops_int
    :no-title:


Multibranch Loops
^^^^^^^^^^^^^^^^^

Functions to evaluate the free energy contributions for mutlibranch loops.

.. doxygengroup:: eval_loops_mb
    :no-title:


.. _api-eval-moves:

Energy Evaluation for Atomic Moves
----------------------------------

Functions to evaluate the free energy change of a structure after application of
(a set of) atomic moves

.. doxygengroup:: eval_move
    :no-title:


.. _api-eval-structures:

Evaluation of Structures
------------------------

Functions to evaluate the free energy of secondary structures under a particular set of
parameters.

.. doxygengroup:: eval
    :no-title:


Deprecated Interface for Free Energy Evaluation
-----------------------------------------------

.. warning::

    Using the functions below is discouraged as they have been marked
    deprecated and will be removed from the library in the (near) future!

.. doxygengroup:: eval_deprecated
    :no-title:
