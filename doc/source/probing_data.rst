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

In the following, you'll find the respective API symbols that allow for the integration
of experimental probing data. In particular, we implement the most commonly used methods
of how such data can be converted into pseudo energies that can then be turned into
*soft constraints*.

.. toctree::
   :maxdepth: 1
   :caption: Specialized Modules:

   probing/SHAPE
   probing/perturbation

Generic Probing Data API
------------------------

.. doxygengroup:: probing_data
    :no-title:
