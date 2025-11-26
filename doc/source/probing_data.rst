Experimental Structure Probing Data
===================================

While RNA secondary structure prediction yields good predictions in general,
the model implemented in the prediction algorithms and its parameters are not
perfect. This may be due to several reasons, such as uncertainties in the parameters
and the simplified assumptions of the model itself. However, prediction performance
can be increased by integrating (experimental) RNA structure probing data, such
as derived from selective 2'-hydroxyl acylation analyzed by primer extension (SHAPE),
dimethyl sulfate (DMS), inline probing, or similar techniques.

Such experimental probing data is usually integrated in the form of small pertubations
in the evaluated energy contributions (:doc:`/grammar/constraints/soft`) that effectively
guide the prediction towards the information gained from the experiment.


.. toctree::
   :maxdepth: 1
   :caption: Specialized Modules:

   probing/strategies
   probing/SHAPE
   probing/perturbation


Structure Probing Data Workflow
-------------------------------

Our implementations follow a particular workflow that basically divides the application
of the structure probing data into three steps:

* Pre-processing of the probing data
* Conversion of the pre-processed data into pseudo energy contributions
* Passing the pseudo energies to the soft-constraints interface

The first two steps are handled by a dedicated **strategy**, that basically is a
single callback function (:c:type:`vrna_probing_strategy_f`). Then, the probing data and
the strategy need to be bundled, e.g. by :c:func:`vrna_probing_data_linear`, to obtain
a data structure of type :c:type:`vrna_probing_data_t`. This data structure is then passed
to the :c:func:`vrna_sc_probing` function that performs the last of the three steps by
calling the strategy to process the data and adding the derived pseudo energy
contributions as soft constraints.

In :doc:`probing/strategies` we already provide a few different strategies that can be
used out-of-the-box, i.e. we provide wrappers that generate the bundled data structure
for a particular strategy. But our interface allows for any user-defined strategy, thus new
strategies can be easily implemented by the user and then passed to our API.

After all the above three steps are done, any prediction (MFE, partition function, etc.)
will acknowledge the probing data and therefore guide the prediction based on it.


In the following, you'll find the respective API symbols that allow for the integration
of experimental probing data. In particular, we implement the most commonly used methods
of how such data can be converted into pseudo energies that can then be turned into
*soft constraints*.


Generic Probing Data API
------------------------

.. doxygengroup:: probing_data
    :no-title:
