:html_theme.sidebar_secondary.remove:

#####################
The ViennaRNA Package
#####################

.. only:: html

    .. image:: gfx/vrna.png
        :height: 200px
        :align: right

.. only:: latex

    .. image:: gfx/vrna.png
        :height: 200px
        :align: center

**Predict, analyze, compare and design RNA secondary structures.**

The ViennaRNA Package offers a comprehensive collection of algorithms
dedicated to RNA secondary structures. **Physics based RNA secondary
structure prediction** is at the heart of the package. It allows for
the prediction of **stable structures** and **structural ensembles**,
to compute **thermodynamic probabilities**, such as for base-pairs,
**compare alternative folds**, and the **design of RNA sequences**
that preferentially fold into a particular structure.

For **individual analyses**, **automated workflows**, or the **development
of new methods**, our software provides both, ready-to-use :ref:`programs`
and an API for ``C``, ``Python`` and ``Perl`` (see :doc:`advanced_usage`).


.. only:: html

    .. grid:: 1 1 2 2
        :gutter:  1 1 2 2

        .. grid-item-card::
            :img-top: _static/gfx/index_getting_started.svg
            :text-align: justify

            **Get Started**
            ^^^

            New to the ViennaRNA Package? This guide describes
            how to **download and install** our software, including
            instructions for **customizing the build process** to
            your requirements.
            +++

            .. button-ref:: install
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :click-parent:
                :tooltip: To the Installation Guide

                To the Installation Guide


        .. grid-item-card::
            :img-top: _static/gfx/index_user_guide.svg
            :text-align: justify

            **User Guides and Tutorials**
            ^^^

            * :doc:`usage/programs`
            * :doc:`usage/tutorials`
            * :doc:`usage/io`
            * :doc:`advanced_usage`

            +++

            .. button-ref:: user_guides
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :tooltip: To the User Guides and Tutorials

                To the User Guides


        .. grid-item-card::
            :img-top: _static/gfx/index_api.svg
            :text-align: justify

            **API Reference**
            ^^^

            Want to learn more about specific ``RNAlib`` functions? The
            reference guide contains a **detailed description** of the
            ``RNAlib`` API. It also describes how the C-functions can
            be accessed in a scripting language like Python.

            +++

            .. button-ref:: api
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :click-parent:
                :tooltip: To the API reference

                To the API reference

        .. grid-item-card::
            :img-top: _static/gfx/index_contribute.svg
            :text-align: justify

            **Contribute**
            ^^^

            Found a bug in our programs? Want to extend the
            existing implementation? Then have a look at our
            contributing guidelines to further **improve the
            ViennaRNA Package**.

            +++

            .. button-ref:: contributors
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :outline:
                :click-parent:
                :tooltip: To the Contributors Guide

                To the Contributors Guide


About the ViennaRNA Package
===========================

The core of the **ViennaRNA Package** (:cite:t:`lorenz:2011,hofacker:1994`)
is formed by a collection of algorithms for the **prediction and comparison**
of RNA secondary structures.

**Physics based RNA secondary structure prediction** is the main functionality
in the package. We provide variants of three different kinds of dynamic
programming algorithms for structure prediction, (i) the **minimum free energy**
algorithm of :cite:t:`zuker:1981` which yields a single optimal structure,
(ii) the **partition function** algorithm of :cite:t:`mccaskill:1990` that allows
for **base pair probability** computation or stochastic backtracking of structure
samples, and (iii) the **suboptimal folding** algorithm of :cite:t:`wuchty:1999`
which generates all suboptimal structures within a given energy range of the
optimal energy. For secondary structure comparison, the package contains
several measures of distance (dissimilarities) using either **string alignment**
or **tree-editing** :cite:t:`shapiro:1990`. Finally, we provide an algorithm to
design sequences with a predefined structure, i.e. **inverse folding**.


.. raw:: latex

    \part{Quickstart}

.. toctree::
    :maxdepth: 3
    :caption: Get Started
    :hidden:

    install

.. raw:: latex

    \part{User Guides and Tutorials}

.. toctree::
    :maxdepth: 3
    :caption: Tutorials
    :hidden:

    user_guides

.. raw:: latex

    \part{Advanced Usage}

.. toctree::
    :maxdepth: 3
    :caption: Advanced User Guide
    :hidden:

    advanced_usage
    api

.. raw:: latex

    \part{Miscellaneous}

.. toctree::
    :maxdepth: 3
    :caption: Miscellaneous
    :hidden:

    What's new <news>
    changelog
    bibliography
    Citation <citing>
    FAQ <faq>
    Contributing <contributors>
    license

