Post-transcriptional Base Modifications
=======================================

Energy parameter corrections for modified bases.

Many RNAs are known to be (heavily) modified post-trasnciptionaly. The best
known examples are tRNAs and rRNAs. To-date, more than 150 different modifications
are listed in the MODOMICS database (http://genesilico.pl/modomics/) :cite:p:`boccaletto:2022`.

Many of the modified bases change the pairing behavior compared to their
unmodified version, affecting not only the pairing partner preference, but
also the resulting stability of the loops the base pairs may form.

Here, we provide a simple soft constraints callback implementation to
correct for some well known modified bases where energy parameters are
available for. This mechanism also supports arbitrary new modified base
energy parameters supplied in JSON format
(see :ref:`io/parameters:modified bases` for details).

.. doxygengroup:: modified_bases
    :no-title:
