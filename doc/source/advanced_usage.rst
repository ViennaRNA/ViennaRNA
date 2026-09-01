*******************
Advanced User Guide
*******************

Introduction
============

This guide will show you how to use our implementations in your own
``C``/``C++``, ``Python 3``, or ``Perl 5`` programs and scripting
pipelines without the need to call programs such as ``RNAfold``.
Interfacing directly with ``RNAlib`` allows you to easily **integrate,
combine, and extend our algorithms** for particular tasks and
applications.

.. seealso::

    Minimal instructions for a technical setup to use our library
    can be found in :ref:`advanced-python-interface` and
    :ref:`advanced-linking-against-rnalib`.


Concepts and Algorithms
=======================

Our library is grouped into several modules, each addressing different
aspects of RNA secondary structure related problems. Below is an overview
of the **concepts and algorithms** we've been implementing.

.. only:: html

    .. grid:: 1 1 2 3
        :gutter: 1 2 3 4

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Minimum Free Energy (MFE)**
            ^^^
            Uncertain how to predict MFE structures, e.g. for **very
            long RNA sequences** or **multiple sequence alignments (MSA)**?
            Answers might be in this section.
            +++
            .. button-ref:: advanced/mfe
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the MFE Algorithms Guide

                MFE Algorithms Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Partition Function**
            ^^^
            Want to predict **thermodynamic properties** within
            the entire structure ensemble, such as
            **base pair probabilities**? Learn more in this section.
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Partition Function and Equilibrium Probabilitiy Guide

                Partition Function Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Suboptimal Structures**
            ^^^
            Want to explore or even enumerate the ensemble of structures?
            Find algorithms for **representative and suboptimal**
            structures here.
            +++
            .. button-ref:: advanced/subopt
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Suboptimals and Representative Structures Guide

                Suboptimal Structures Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Free Energy Evaluation**
            ^^^
            Want to assess the **free energy** of a particular loop
            or secondary structure **given a set of energy parameters**?
            This section will help you learn how to do so!
            +++

            .. button-ref:: advanced/eval
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Free Energy Evaluation Guide

                Energy Evaluation Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **RNA Folding Grammar**
            ^^^
            Require **constraints** for the prediction algorithms or
            want to add **auxiliary RNA folding grammar rules**?
            Find solutions how to accomplish that here!
            +++
            .. button-ref:: advanced/grammar
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the RNA Folding Grammar Guide

                Folding Grammar Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **RNA-RNA Interaction**
            ^^^
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the RNA-RNA Interaction Guide

                RNA-RNA Interaction Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Classified Dynamic Programming Variants**
            ^^^
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Classified Dynamic Programming Variants Guide

                Classified Dynamic Programming Variants Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Inverse Folding (Design)**
            ^^^
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Inverse Folding (Design) Guide

                Inverse Folding (Design) Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Experimental Structure Probing Data**
            ^^^
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Experimental Structure Probing Data Guide

                Experimental Structure Probing Data Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Post-transcriptional Base Modifications**
            ^^^
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Post-transcriptional Base Modifications Guide

                Post-transcriptional Base Modifications Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Ligands Binding to RNA Structures**
            ^^^
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Ligands Binding to RNA Structures Guide

                Ligands Binding to RNA Structures Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Structure Modules and Pseudoknots**
            ^^^
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Structure Modules and Pseudoknots Guide

                Structure Modules and Pseudoknots Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **The Secondary Structure Landscape**
            ^^^
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the The RNA Secondary Structure Landscape Guide

                The RNA Secondary Structure Landscape Guide

        .. grid-item-card::
            :text-align: justify

            :octicon:`book;1em;sd-text-primary sd-font-weight-bold` **Utilities**
            ^^^
            +++
            .. button-ref:: advanced/partfunc
                :ref-type:  doc
                :color: primary
                :align: center
                :expand:
                :click-parent:
                :outline:
                :tooltip: To the Utilities Guide

                Utilities Guide



.. toctree::
    :maxdepth: 1
    :caption: Modules:
    :hidden:

    advanced/eval
    advanced/grammar
    advanced/mfe
    advanced/partfunc
    advanced/subopt
    landscape
    cofold
    class_fold
    inverse
    probing_data
    ligand_binding
    paired_modules
    modified_bases
    utils


.. _advanced-linking-against-rnalib:

Linking against RNAlib
======================

To use our implemented algorithms you simply need to link your program
to our ``RNAlib`` C-library that usually comes along with the ViennaRNA
Package installation. If you've installed the ViennaRNA Package as a
pre-build binary package, you probably need the corresponding development
package, e.g. ``viennarna-devel``, or ``viennarna-dev``. The only thing
that is left is to include the ViennaRNA header files into your source
code, e.g.:

.. code:: c

  #include <ViennaRNA/mfe.h>

and start using our fast and efficient algorithm implementations.

.. seealso::

    In the :doc:`/examples/c` section, we list a small set of example code that
    usually is a good starting point for your application.

Compiler and Linker flags
-------------------------

Of course, simply adding the ViennaRNA header files into your source code is
usually not enough. You probably need to tell your compiler where to find the
header files, and sometimes add additional pre-processor directives. Whenever
your installation of *RNAlib* was build with default settings and the header
files were installed into their default location, a simple::

  -I/usr/include

pre-processor/compile flag should suffice. It can even be omitted in this case,
since your compiler should search this directory by default anyway. You only
need to change the path from ``/usr/include`` to the correct location whenever
the header files have been installed into a non-standard directory.

If you've compiled *RNAlib* with some non-default settings then you probably
need to define some additional pre-processor macros:

* ``VRNA_DISABLE_C11_FEATURES`` ... Disable C11/C++11 features.

  .. warning::
  
    Add this directive to your pre-processor/compile flags only if ``RNAlib``
    was build with the ``--disable-c11`` configure option.

  .. admonition:: See also...
  
    :ref:`configuration:disable c11/c++11 features` and :c:func:`vrna_C11_features()`

* ``VRNA_WARN_DEPRECATED`` ... Enable warnings for using deprecated symbols.

  .. note::
  
    Adding this directive enables compiler warnings whenever you use symbols
    in *RNAlib* that are marked *deprecated*.

  .. admonition:: See also...
  
    :ref:`configuration:deprecated symbols` and :doc:`api/deprecated`

* ``USE_FLOAT_PF`` ... Use single precision floating point operations instead
  of double precision in partition function computations.

  .. warning::
  
    Define this macro only if *RNAlib* was build with the ``--enable-floatpf``
    configure option!

  .. admonition:: See also...

    :ref:`configuration:Single precision`


For instance, you might want to add the following definition(s) to your
pre-processor/compile flags::

  -DVRNA_DISABLE_C11_FEATURES

Finally, linking against *RNAlib* is achieved by adding the following linker flag::

  -L/usr/lib -lRNA -flto -fopenmp

Again, the path to the library, ``/usr/lib``, may be omitted if this path is
searched for libraries by default. The second flag tells the linker to include
``libRNA.a``, and the remaining two flags activate :ref:`configuration:link time optimization`
and :ref:`configuration:openmp` support, respectively.

.. note::

  Depending on your linker, the last two flags may differ.

  Depending on your configure time decisions, you can drop one or both of the last flags.

  In case you've compiled *RNAlib* with LTO support (See :ref:`configuration:link time
  optimization`) and you are using a different compiler for your third-party project that
  links against our library, you may add the ``-fno-lto`` flag to disable Link Time
  Optimization.

.. admonition:: See also...

  :ref:`faq:linking fails with lto error`


The pkg-config tool
-------------------

Instead of hard-coding the required compiler and linker flags, you can also let the
``pkg-config`` tool automatically determine the required flags. This tool is usually
packaged for any Linux distribution and should be available for MacOS X and MinGW as
well. We ship a file ``RNAlib2.pc`` which is installed along with the static ``libRNA.a``
C-library and populated with all required compiler and linker flags that correspond to
your configure time decisions.

The compiler flags required for properly building your code that uses *RNAlib* can be
easily obtained via::

  pkg-config --cflags RNAlib2

You get the corresponding linker flags using::

  pkg-config --libs RNAlib2

With this widely accepted standard it is also very easy to integrate *RNAlib* in your
``autotools`` project, just have a look at the ``PKG_CHECK_MODULES`` macro.


.. _advanced-swig-wrappers:

SWIG Wrappers
=============

For an easy integration into scripting languages, we provide an automatically
generated interface to the RNAlib C-library, generated with SWIG. Currently,
we support ``Perl 5`` and ``Python`` as target languages.

.. admonition:: See also...

  :doc:`/api/python`

Function Renaming
-----------------

To provide a namespace-like separation of function symbols from our C library and
third-party code, we use the prefix ``vrna_`` or ``VRNA_`` whenever possible. This,
however, is not necessary for the scripting language interface, as it uses the
separate namespace or package `RNA` anyway. Consequently, symbols that appear to
have the ``vrna_`` or ``VRNA_`` prefix in the C-library have the corresponding prefix
stripped away.

For instance, the C code

.. code:: c

  mfe = vrna_fold(sequence, structure);

translates to

.. code:: perl

  my ($structure, $mfe) = RNA::fold($sequence)

in the Perl 5 interface, and

.. code:: python

  structure, mfe = RNA.fold(sequence)

for Python. Note, that in this example we also make use of the possibilty to
return multiple data at once in the scripting language, while the C library function
uses additional parameters to return multiple data.

Functions that are dedicated to work on specific data structures only,
e.g. the :c:type:`vrna_fold_compound_t`, are usually not exported at all. Instead,
they are attached as object methods of a corresponding class
(see :ref:`wrappers:object oriented interface` for detailed information).


Object Oriented Interface
-------------------------


For data structures, typedefs, and enumerations the ``vrna_`` prefixes are
dropped as well, together with their suffixes ``_s``, ``_t``, and ``_e``, respectively.
Furthermore, data structures are usually transformed into classes and
relevant functions of the C-library are attached as methods.

.. tip::

    Examples on the basic usage of the scripting language interfaces can be
    found in the :doc:`/examples/perl5` and :doc:`/examples/python` sections.

.. seealso::

    The full list of specific differences between the C API and the
    SWIG Wrapper for the supported scripting languages ca be found at
    :doc:/api/swig_wrappers`.


.. _advanced-python-interface:

Python Interface
================

Almost all symbols of the API available in our *RNAlib* C-library is wrapped for
use in Python using ``swig``. That makes our fast and efficient algorithms and
tools available for third-party Python programs and pipelines.

.. tip::

    The Python interface is usually part of the installation of the ViennaRNA Package,
    see also :ref:`install` and :doc:`/install:scripting-language-interfaces`.

    If for any reason your installation does not provide our Python interface or
    in cases where you don't want to install the full ViennaRNA Package but *only*
    the Python bindings to *RNAlib*, you may also install them via Pythons ``pip``:

    .. code:: bash

        python -m pip install viennarna

Usage
-----

To use our Python bindings simply ``import`` the ``RNA`` or ``ViennaRNA`` package
like

.. code:: python

  import RNA

or

.. code:: python

  import ViennaRNA

The ``RNA`` module that provides access to our *RNAlib* C-library can also be imported
directly using

.. code:: python

  from RNA import RNA

or

.. code:: python

  from ViennaRNA import RNA

.. note::

  In previous release of the ViennaRNA Packge, only the ``RNA`` package/module has
  been available.
  Since version 2.6.2 we maintain the `ViennaRNA <https://pypi.org/project/ViennaRNA/>`_
  project at https://pypi.org. The former maintainer additionally introduced the
  ``ViennaRNA`` package which we intend to keep and extend in future releases.


Global Variables
----------------

For the Python interface(s) SWIG places global variables of the C-library
into an additional namespace ``cvar``. For instance, changing the global temperature
variable thus becomes

.. code:: python

  RNA.cvar.temperature = 25


Pythonic interface
------------------

Since our library is written in ``C`` the functions we provide in our API might
seem awkward for users more familiar with Pythons object oriented fashion. Therefore,
we spend some effort on creating a more *pythonic* interface here. In particular, we
tried to group together particular data structures and functions operating on them to
derive classes and objects with corresponding methods attached.

If you browse through our :doc:`reference manual </api>`, many C-functions have additional
:doc:`/api/swig_wrapper` in their description. These descriptions should give an idea how
the function is available in the Python interface. Usually, our ``C`` functions,
data structures, typedefs, and enumerations use the ``vrna_`` prefixes and ``_s``,
``_t``, ``_e`` suffixes. Those decorators are useful in ``C`` but of less use in the
context of Python packages or modules. Therefore, these prefixes and suffixes are
*dropped* from the Python interface.

Object orientation
^^^^^^^^^^^^^^^^^^

Consider the C-function :c:func:`vrna_fold_compound`. This creates a
:c:type:`vrna_fold_compound_t` data structure that is then passed around to various
functions, e.g. to :c:func:`vrna_mfe` from :doc:`advanced/mfe`.

A corresponding C-code may look like this:

.. code:: c

  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>

  #include <ViennaRNA/utils/basic.h>
  #include <ViennaRNA/fold_compound.h>
  #include <ViennaRNA/mfe.h>
  
  int
  main(int  argc,
       char *argv[])
  {
    char *seq, *ss;
    float mfe;
    vrna_fold_compound_t *fc;

    seq = "AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUG";
    ss  = vrna_alloc(sizeof(char) * (strlen(seq) + 1));
    fc  = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
    mfe = vrna_mfe(fc, ss);

    printf("%s\n%s (%6.2f)\n", seq, ss, mfe);

    return EXIT_SUCCESS;
  }

In our Python interface, the :c:type:`vrna_fold_compound_t` data structure becomes the
:py:class:`RNA.fold_compound` class, the :c:func:`vrna_fold_compound` becomes one of
its constructors and the :c:func:`vrna_mfe` function becomes the method
:py:meth:`RNA.fold_compound.mfe()`. So, the Python code would probably translate
to something like

.. code:: python

  import RNA

  seq = "AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUG"
  fc  = RNA.fold_compound(seq)
  (ss, mfe) = fc.mfe()

  print(f"{seq}\n{ss} ({mfe:6.2f})")

.. note::

  The C-function :c:func:`vrna_mfe` actually returns two values, the MFE in units
  of :math:`\text{kcal} \cdot \text{mol}^{-1}` and the corresponding MFE structure.
  The latter is written to the ``ss`` pointer. This is necessary since ``C`` functions
  can at most return one single value. In Python, function and methods may return
  arbitrarily many values instead, and in addition, passing parameters to a function
  or method such that it changes its content is generally discouraged. Therefore,
  our functions that return values through function parameters usually return them
  *regularly* in the Python interface.

Lists and Tuples
^^^^^^^^^^^^^^^^

C-functions in our API that return or receive list-like data usually utilize *pointers*.
Since there are no such things in Python, they would be wrapped as particular kind of
objects that would then be tedious to work with. For the Python interface, we therefore
tried to wrap the majority of these instances to *native* Python types, such as ``list``
or ``tuple``. Therefore, one can usually pass a ``list`` to a function that uses pointers
to array in ``C``, and expect to receive a ``list`` or ``tuple`` from functions that return
pointers to arrays.

Energy Parameters
-----------------

Energy parameters are compiled into our library, so there is usually no necessity
to load them from a file. All parameter files shipped with the ViennaRNA Package
can be loaded by simply calling any of the dedicated functions:

* :py:func:`RNA.params_load_RNA_Turner2004()` (default RNA parameters)
* :py:func:`RNA.params_load_DNA_Mathews2004()` (default DNA parameters)
* :py:func:`RNA.params_load_DNA_Mathews1999()` (old DNA parameters)
* :py:func:`RNA.params_load_RNA_Turner1999()` (old RNA parameters)
* :py:func:`RNA.params_load_RNA_Andronescu2007()` (trained RNA parameters)
* :py:func:`RNA.params_load_RNA_Langdon2018()` (trained RNA parameters)
* :py:func:`RNA.params_load_RNA_misc_special_hairpins()` (special hairpin loop parameters)


.. tip::

    A few more Python code examples can be found :doc:`here </examples/python>`.

.. seealso::

    :doc:`/api/python` contains a full Python API description.


Callback Functions
==================

Many functions in ``RNAlib`` support so-called **callback mechanisms** to propagate
the computation results. In essence, this means that a user defines a function
that takes computation results as input and does whatever should be done with
these results. Then, this function is provided to our algorithms in ``RNAlib``
such that they can call the function and provide corresponding data as soon
as it has been computed.

Why Callbacks?
--------------

Using callback mechanisms, our library enables users not only to retrieve computed
data without the need for parsing complicated data structures, but also allows one
to tweak our implementation to do additional tasks without the requirement of a
re-implementation of basic algorithms.

Our implementation of the callback mechanisms always follows the same scheme:

The user ...

* ... defines a function that complies with the interface we've defined, and
* ... passes a pointer to said function to our implementations

In addition to the specific arguments of our callback interfaces, virtually all callbacks
receive an additional *pass-through-pointer* as their last argument. This enables one
to ...

* ... encapsulate data, and
* ... provide thread-safe operations,

since this pointer is simply passed through by our library functions. It may therefore hold
the address of an arbitrary, user-defined data structure.

Scripting Language Support
--------------------------

Our callback mechanisms also work in a cross-language specific manner. This
means that users can write functions (or methods) in a scripting language,
such as Python, and use them as callbacks for our C-library. Our scripting
language interface will take care of transforming the relevant data structures
from the target language into C and vice-versa. Whenever our algorithms trigger
the callback, they will then be calling the actual scripting language function
and provide the corresponding data directly to them.

.. warning::

    Keep in mind that the translation between the scripting language and C
    involves many extra function calls to prepare and evaluate the corresponding
    data structures. This in turn impacts the runtime of our algorithms that
    can be substantial. For instance providing callbacks for the hard- or
    soft constraints framework from a scripting language can lead to a slow-down
    of up to a factor of 10.



RNAlib API v3.0
===============

With version 2.2 we introduce the new API that will take over the old one
in the future version 3.0. By then, backwards compatibility will be broken, and
third party applications using RNAlib need to be ported. This switch of API became
necessary, since many new features found their way into the RNAlib where a balance
between threadsafety and easy-to-use library functions is hard or even impossible
to establish. Furthermore, many old functions of the library are present as
slightly modified copies of themself to provide a crude way to overload functions.

Therefore, we introduce the new v3.0 API very early in our development stage
such that developers have enough time to migrate to the new functions and interfaces.
We also started to provide encapsulation of the RNAlib functions, data structures,
typedefs, and macros by prefixing them with ``vrna_`` and ``VRNA_`` , respectively.

Header files should also be included using the ``ViennaRNA/`` namespace, e.g.

.. code:: c

  #include <ViennaRNA/fold.h>

instead of just using

.. code:: c

  #include <fold.h>

as is has been required for RNAlib 1.x and 2.x.

This eases the work for programmers of third party applications that would otherwise
need to put much effort into renaming functions and data types in their own
implementations if their names appear in our library. Since we still provide backward
compatibility up to the last version of RNAlib 2.x, this advantage may be fully
exploited only starting from v3.0 which will be released in the future. However, our
plan is to provide the possibility for an early switch-off mechanism of the backward
compatibility in one of our next releases of ViennaRNA Package 2.x.

.. tip::

    Examples on how to use the new v3.0 API can be found in the :doc:`/examples/c`
    section.


Example Code
============

.. toctree::
    :maxdepth: 1

    examples/c
    examples/python
    examples/perl5
