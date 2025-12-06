Strategies for Linear Data
==========================

Recent literature lists several methods and strategies that deal with the
problem to convert linear experimental RNA structure probing data into pseudo
energy terms to guide RNA secondary strucure predictions. Such data may come
from SHAPE, DMS, lead, inline probing, or similar techniques. The most commonly
used method may be the one presented in :cite:t:`deigan:2009` where the free
energy evaluation for stacked base pairs adds a pseudo energy contribution
derived from the reactivity values of the probing data. The higher the
reactivity the stronger the stack is penalized by a positive energy.


.. contents:: Table of Contents
    :local:
    :depth: 2


Generic Probing Data Strategy API
---------------------------------

.. doxygengroup:: probing_data_strategy
    :no-title:


The Deigan 2009 Strategy API
----------------------------

.. doxygengroup:: probing_data_strategy_deigan
    :no-title:


The Zarringhalam 2012 Strategy API
----------------------------------

.. doxygengroup:: probing_data_strategy_zarringhalam
    :no-title:


The Eddy 2014 Strategy API
--------------------------

.. doxygengroup:: probing_data_strategy_eddy
    :no-title:
