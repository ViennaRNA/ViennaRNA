:html_theme.sidebar_secondary.remove:

.. raw:: html

    <style type="text/css">
        h1 {display: none;}
        .bd-main .bd-content .bd-article-container {max-width: 100%;}
        .big-font {font-size: var(--pst-font-size-h4); color: var(--pst-color-primary)}
        .bd-footer-article{display: none;}
    </style>

#####################
The ViennaRNA Package
#####################

.. div:: sd-text-justify

    .. only:: html

        .. image:: gfx/vrna_logo_light.svg
            :class: only-light
            :height: 200px
            :align: right

        .. image:: gfx/vrna_logo_dark.svg
            :class: only-dark
            :height: 200px
            :align: right

    .. only:: latex

        .. image:: gfx/vrna_logo_light.svg
            :height: 200px
            :align: center

    .. rst-class:: font-weight-bold big-font

        Predict, analyze, compare and design RNA secondary structures

    The ViennaRNA Package offers a comprehensive collection of algorithms
    dedicated to RNA secondary structures. **Physics based RNA secondary
    structure prediction** is at the heart of the package. It allows for
    the prediction of **stable structures** and **structural ensembles**,
    to compute **thermodynamic probabilities**, such as for base-pairs,
    to **compare alternative folds**, and to **design RNA sequences**
    that preferentially fold into a particular structure.

    For **individual analyses**, **automated workflows**, or the **development
    of new methods**, our software provides both, ready-to-use :ref:`programs`
    and an API for ``C``, ``Python`` and ``Perl`` (see :doc:`advanced_usage`).


.. only:: html

    .. grid:: 1 1 2 4
        :gutter:  1 2 2 4

        .. grid-item-card::
            :img-top: _static/gfx/index_getting_started.svg
            :text-align: justify

            **Get Started**
            ^^^
            New to the ViennaRNA Package? Find **download
            and install** instructions as well as steps
            to **customize the builds** for your requirements
            here.
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
            ``RNAlib`` API.
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
            Found a bug in our programs? Want to extend an
            existing implementation? Our contributors guide
            shows how to **improve the ViennaRNA Package**.
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
    About <about>
    license

