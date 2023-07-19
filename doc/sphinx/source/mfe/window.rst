Local (sliding window) MFE Prediction
=====================================

Variations of the local (sliding window) Minimum Free Energy (MFE) prediction algorithm.

We provide implementations for the local (sliding window) MFE prediction algorithm for

* Single sequences,
* Multiple sequence alignments (MSA), and

Note, that our implementation scans an RNA sequence (or MSA) from the 3' to the 5'
end, and reports back locally optimal (consensus) structures, the corresponding free
energy, and the position of the sliding window in global coordinates.

For any particular RNA sequence (or MSA) multiple locally optimal (consensus)
secondary structures may be predicted. Thus, we tried to implement an interface that
allows for an effortless conversion of the corresponding hits into any target data
structure. As a consequence, we provide two distinct ways to retrieve the corresponding
predictions, either

* through directly writing to an open ``FILE`` stream on-the-fly, or
* through a callback function mechanism.

The latter allows one to store the results in any possible target data structure. Our
implementations then pass the results through the user-implemented callback as soon as
the prediction for a particular window is finished.

.. doxygengroup:: mfe_window
    :no-title:
