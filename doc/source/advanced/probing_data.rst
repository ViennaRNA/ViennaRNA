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
in the evaluated energy contributions (:ref:`advanced-constraints-soft`) that effectively
guide the prediction towards the information gained from the experiment.


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

We already implemented a few different strategies that can be used out-of-the-box, i.e.
we provide wrappers that generate the bundled data structure for a particular strategy.
The available API symbols can be found in :doc:`probing/strategies`.

However, the generic interface allows for any user-defined strategy, thus new
strategies can be easily implemented by the user and then passed to our API.

After all the above three steps are done, any prediction (MFE, partition function, etc.)
will acknowledge the probing data and therefore guide the prediction based on it.


In the following, you'll find the respective API symbols that allow for the integration
of experimental probing data.


Strategies for Linear Data
--------------------------

Recent literature lists several methods and strategies that deal with the
problem to convert linear experimental RNA structure probing data into pseudo
energy terms to guide RNA secondary strucure predictions. Such data may come
from SHAPE, DMS, lead, inline probing, or similar techniques. The most commonly
used method may be the one presented in :cite:t:`deigan:2009` where the free
energy evaluation for stacked base pairs adds a pseudo energy contribution
derived from the reactivity values of the probing data. The higher the
reactivity the stronger the stack is penalized by a positive energy.



SHAPE Reactivity Data
---------------------

Incorporate SHAPE reactivity structure probing data into the folding
recursions by means of soft constraints.

Details for our implementation to incorporate SHAPE reactivity data to guide
secondary structure prediction can be found in :cite:t:`lorenz:2016a`.


Generate Soft Constraints from Data
-----------------------------------

Find a vector of perturbation energies that minimizes the discripancies between
predicted and observed pairing probabilities and the amount of neccessary adjustments.

