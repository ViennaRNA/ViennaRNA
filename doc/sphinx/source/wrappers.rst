Scripting Language Interface(s)
===============================

Introduction
------------

For an easy integration into scripting languages, we provide an automatically
generated interface to the RNAlib C-library, generated with SWIG.


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

Global Variables
^^^^^^^^^^^^^^^^

For the Python interface(s) SWIG places global variables of the C-library
into an additional namespace ``cvar``. For instance, changing the global temperature
variable thus becomes

.. code: python

  RNA.cvar.temperature = 25


Object Oriented Interface
-------------------------


For data structures, typedefs, and enumerations the ``vrna_`` prefixes are
dropped as well, together with their suffixes ``_s``, ``_t``, and ``_e``, respectively.
Furthermore, data structures are usually transformed into classes and
relevant functions of the C-library are attached as methods.

Examples
--------

Examples on the basic usage of the scripting language interfaces can be
found in the :doc:`/examples/examples_perl5` and :doc:`/examples/examples_python` sections.

SWIG Wrapper notes
------------------

Special notes on how functions, structures, enums, and macro definitions are
actually wrapped, can be found below

.. doxygenpage:: wrappers
   :content-only:

