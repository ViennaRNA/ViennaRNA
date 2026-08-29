:html_theme.sidebar_secondary.remove: true

.. RNA documentation master file, created by
   sphinx-quickstart on Thu Jul 28 23:52:58 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#####################
The ViennaRNA Package
#####################

.. only:: html

    .. image:: gfx/vrna.png
        :height: 200px
        :align: left

.. only:: latex

    .. image:: gfx/vrna.png
        :height: 200px
        :align: center

The core of the **ViennaRNA Package** (:cite:t:`lorenz:2011,hofacker:1994`)
is formed by a collection of algorithms for the **prediction and comparison**
of RNA secondary structures. These implementations can be accessed through
stand-alone programs, such as ``RNAfold``, ``RNAdistance``, etc., which
most often are sufficient for regular users.

For those who wish to develop their own programs we provide ``RNAlib``,
a C-library that can be linked against your own programs or used in your
Python or Perl scripts and pipelines. The latter can be easily accomplished
through our :doc:`/wrappers` for :doc:`Python <api_python>` and Perl 5.

.. only:: html

    .. grid:: 1 1 2 2
        :gutter:  1 1 2 2

        .. grid-item-card::
            :img-top: _static/gfx/index_getting_started.svg

            **Get Started**
            ^^^

            * :doc:`Download and Installation <install>`
            * :doc:`programs`
            * :doc:`using_rnalib`
            * :doc:`io`

            +++

            .. button-ref:: get_started
                :ref-type:  doc
                :color: primary
                :align: center

                To the Install and Usage Guide


        .. grid-item-card::
            :img-top: _static/gfx/index_user_guide.svg

            **Advanced User Guide**
            ^^^

            Our library consists of several modules, each addressing
            different aspects of RNA secondary structure related problems.
            Here you find an overview of the *concepts and algorithms* for
            our implementations.

            +++

            .. button-ref:: advanced_usage
                :ref-type:  doc
                :color: primary
                :align: center

                To the Advanced User Guide


        .. grid-item-card::
            :img-top: _static/gfx/index_api.svg

            **API Reference**
            ^^^

            Want to learn more about specific ``RNAlib`` functions? The
            reference guide contains a detailed description of the
            ``RNAlib`` API. It also describes how the C-functions can
            be accessed in a scripting language like Python.

            +++

            .. button-ref:: api
                :ref-type:  doc
                :color: primary
                :align: center

                To the API reference

        .. grid-item-card::
            :img-top: _static/gfx/index_contribute.svg

            **Contribute**
            ^^^

            Found a bug in our programs? Want to extend the
            existing implementation? Then have a look at our
            contributing guidelines to further improve the
            ViennaRNA Package.

            +++

            .. button-ref:: contributors
                :ref-type:  doc
                :color: primary
                :align: center

                To the Contributors Guide


.. raw:: latex

    \part{Quickstart}

.. toctree::
    :maxdepth: 3
    :caption: Get Started
    :hidden:

    get_started

.. raw:: latex

    \part{User Guides}

.. toctree::
    :maxdepth: 3
    :caption: Tutorials
    :hidden:

    tutorials

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

