*******************
Advanced User Guide
*******************

Our library is grouped into several modules, each addressing different
aspects of RNA secondary structure related problems. This is an overview
of the concepts and algorithms for which implementations can be found
in this library.

.. only:: html

    .. grid:: 1 1 2 2
        :gutter: 1 1 2 3

        .. grid-item-card::

            **Free Energy Evaluation**
            ^^^
            Want to assess the **free energy** of a particular loop
            or en entire secondary structure **given a set of
            energy parameters** and environmental conditions?
            This section will help you learn how to do so!
            +++
            .. button-ref:: advanced/eval
                :ref-type:  doc
                :color: primary
                :align: center

                To the Energy Evaluation Guide

        .. grid-item-card::

            **RNA Folding Grammar**
            ^^^
            Require certain **constraints** to be applied to the
            prediction algorithms or want to use or add **auxiliary
            RNA folding grammar rules**? Find soultions how to accomplish
            that here!
            +++
            .. button-ref:: advanced/grammar
                :ref-type:  doc
                :color: primary
                :align: center

                To the Folding Grammar Guide

        .. grid-item-card::

            **Minimum Free Energy (MFE)**
            ^^^
            Uncertain how to predict MFE structures, e.g. for **very
            long RNA sequences**, multiple interacting RNAs, or
            **multiple sequence alignments (MSA)**? There might be
            answers to your questions in this section.
            +++
            .. button-ref:: advanced/mfe
                :ref-type:  doc
                :color: primary
                :align: center

                To the MFE Algorithms Guide

        .. grid-item-card::

            **Partition Function and Equilibrium Probabilities**
            ^^^
            Want to compute **equilibrium probabilities** for certain
            structural features, such as base pairs or helices or
            entire secondary structures? Have a look at the
            partition function and equilibrium probability
            section!
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center

                To the Partition Function Guide


.. toctree::
    :maxdepth: 1
    :caption: Modules:
    :hidden:

    advanced/eval
    advanced/grammar
    advanced/mfe
    advanced/partfunc
    landscape
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


