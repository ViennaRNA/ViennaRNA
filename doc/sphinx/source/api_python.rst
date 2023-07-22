Python API
==========

Almost all symbols of the API available in our *RNAlib* C-library is wrapped for
use in Python using ``swig``. That makes our fast and efficient algorithms and
tools available for third-party Python programs and scripting languages.

.. note::

  Our Python API is automatically generated and translated from our C-library documentation.
  If you find anything problematic or want to to help us improve the documentation, do
  not hesitate to contact us or make a PR at our `official github repository <https://github.com/ViennaRNA/ViennaRNA>`_.

Installation
------------

The Python interface is usually part of the installation of the ViennaRNA Package,
see also :doc:`/install` and :ref:`configuration:scripting language interfaces`.

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
*SWIG Wrapper Notes* in their description. These descriptions should give an idea how
the function is available in the Python interface. Usually, our ``C`` functions,
data structures, typedefs, and enumerations use the ``vrna_`` prefixes and ``_s``,
``_t``, ``_e`` suffixes. Those decorators are useful in ``C`` but of less use in the
context of Python packages or modules. Therefore, these prefixes and suffixes are
*dropped* from the Python interface.

Object orientation
^^^^^^^^^^^^^^^^^^

Consider the C-function :c:func:`vrna_fold_compound`. This creates a
:c:type:`vrna_fold_compound_t` data structure that is then passed around to various
functions, e.g. to :c:func:`vrna_mfe` to compute the :doc:`MFE structure </mfe/global>`.
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


Examples
--------

A few more Python code examples can be found :doc:`here </examples/python>`.

The ``RNA`` Python module
-------------------------

.. automodule:: RNA
   :imported-members:
   :members:
   :undoc-members:
   :show-inheritance:
