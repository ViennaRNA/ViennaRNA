Experimental Structure Probing Data
===================================

In the following, you'll find the respective API symbols that allow for the integration
of experimental probing data.

Generic Probing Data API
------------------------

.. doxygengroup:: probing_data
    :no-title:


Strategies for Linear Data
--------------------------

Generic Probing Data Strategy API
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: probing_data_strategy
    :no-title:


The Deigan 2009 Strategy API
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: probing_data_strategy_deigan
    :no-title:


The Zarringhalam 2012 Strategy API
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: probing_data_strategy_zarringhalam
    :no-title:


The Eddy 2014 Strategy API
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: probing_data_strategy_eddy
    :no-title:


SHAPE Reactivity Data
---------------------

Incorporate SHAPE reactivity structure probing data into the folding
recursions by means of soft constraints.

.. doxygengroup:: SHAPE_reactivities
    :no-title:


Generate Soft Constraints from Data
-----------------------------------

Find a vector of perturbation energies that minimizes the discripancies between
predicted and observed pairing probabilities and the amount of neccessary adjustments.

.. doxygengroup:: perturbation
    :no-title:
