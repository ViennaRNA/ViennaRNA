**********
User Guide
**********

Our implementations can be accessed through stand-alone programs, such as
``RNAfold``, ``RNAdistance``, etc., which most often are sufficient for
regular users. For those who wish to develop their own programs we provide
``RNAlib``, a ``C``-library that can be linked against your own programs or used
in your ``Python`` or ``Perl`` scripts and pipelines. The latter can be easily
accomplished through our :ref:`advanced-swig-wrappers` for ``Python 3`` and
``Perl 5``, see also :ref:`advanced-python-interface`.


.. only:: html

    .. grid:: 1 1 2 2
        :gutter: 1 1 2 3

        .. grid-item-card:: The Programs of the ViennaRNA Package
            :link:  programs
            :link-type: ref

            The ViennaRNA Package ships several executable programs for various purposes.
            A more or less elaborate collection of tutorials and manuals on how to use our
            software can be found here.

        .. grid-item-card:: Supported File Formats
            :link:  io-formats
            :link-type: ref

            Here, we introduce the most common notations of sequence and structure data, specifications
            of bioinformatics sequence and structure file formats, and various input and output file
            formats used by our library.

        .. grid-item-card:: Tutorials
            :link:  tutorials
            :link-type: ref

            Find more elaborate instructions how to use our executable programs here.

.. toctree::
    :caption: User Guides and Tutorials
    :maxdepth: 1
    :hidden:

    usage/programs
    usage/manpages
    usage/io
    usage/tutorials

