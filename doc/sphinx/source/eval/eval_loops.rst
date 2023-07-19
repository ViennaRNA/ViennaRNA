Energy Evaluation for Individual Loops
======================================

To assess the free energy contribution of a particular loop :math:`L` within a
secondary structure, two variants are provided

* The *bare* free energy :math:`E_L` (usually in units of deka-calories, i.e.
  multiples of :math:`10 \text{cal} \cdot \text{mol}^{-1}`, and
* The *Boltzmann weight* :math:`q = exp(-\beta E_L)` of the free energy :math:`E_L`
  (with :math:`\beta = \frac{1}{RT}`, gas constant :math:`R` and temperature :math:`T`)

The latter is usually required for partition function computations.

.. contents:: Table of Contents
    :local:


General
-------

.. doxygengroup:: eval_loops
    :no-title:


Exterior Loops
--------------

Functions to evaluate the free energy contributions for exterior (external) loops.

.. doxygengroup:: eval_loops_ext
    :no-title:


Hairpin Loops
-------------

Functions to evaluate the free energy contributions for hairpin loops.

.. doxygengroup:: eval_loops_hp
    :no-title:


Internal Loops
--------------

Functions to evaluate the free energy contributions for internal (interior) loops.

.. doxygengroup:: eval_loops_int
    :no-title:


Multibranch Loops
-----------------

Functions to evaluate the free energy contributions for mutlibranch loops.

.. doxygengroup:: eval_loops_mb
    :no-title:

