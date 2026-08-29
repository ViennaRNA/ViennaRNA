===================
Advanced User Guide
===================

Our library is grouped into several modules, each addressing different
aspects of RNA secondary structure related problems. This is an overview
of the concepts and algorithms for which implementations can be found
in this library.

Almost all of them rely on the physics based Nearest Neighbor Model
for RNA secondary structure prediction.

.. toctree::
    :maxdepth: 1
    :caption: Modules:
    :hidden:

    advanced/eval
    advanced/grammar
    advanced/mfe
    landscape
    pf_fold
    subopt_and_representatives
    cofold
    class_fold
    inverse
    probing_data
    ligand_binding
    paired_modules
    modified_bases
    utils

.. toctree::
    :maxdepth: 1
    :caption: General Remarks:
    :hidden:

    wrappers
    vrna_3.0


.. grid:: 1 1 2 2
    :gutter: 1 1 2 3

    .. grid-item-card::

        **Free Energy Evaluation**
        ^^^
        Want to assess the free energy of a particular loop
        or en entire secondary structure? Want to change the
        underlying set of parameters, e.g. according to
        environmental conditions?
        +++
        .. button-ref:: advanced/eval
            :ref-type:  doc
            :color: primary
            :align: center

            To the Energy Evaluation Guide

    .. grid-item-card::

        **RNA Folding Grammar**
        ^^^
        +++
        .. button-ref:: advanced/grammar
            :ref-type:  doc
            :color: primary
            :align: center

            To the Folding Grammar Guide

    .. grid-item-card::

        **Minimum Free Energy (MFE)**
        ^^^
        +++
        .. button-ref:: advanced/mfe
            :ref-type:  doc
            :color: primary
            :align: center

            To the MFE Algorithms Guide

