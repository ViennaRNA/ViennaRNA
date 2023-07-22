.. RNA documentation master file, created by
   sphinx-quickstart on Thu Jul 28 23:52:58 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

A Library for predicting and comparing RNA secondary structures

Introduction
------------

The core of the ViennaRNA Package (:cite:t:`lorenz:2011,hofacker:1994`)
is formed by a collection of routines
for the prediction and comparison of RNA secondary structures. These
routines can be accessed through stand-alone programs, such as *RNAfold*,
*RNAdistance* etc., which should be sufficient for most users. For those who
wish to develop their own programs we provide a library which can be  
linked to your own code.

This document describes the library and will be primarily useful to
programmers. However, it also contains details about the implementation
that may be of interest to advanced users. The stand-alone programs are
described in separate man pages. The latest version of the package
including source code and html versions of the documentation can be found
at

http://www.tbi.univie.ac.at/RNA

.. toctree::
   :caption: Installation
   :maxdepth: 1

   install
   configuration

.. toctree::
   :caption: Usage:
   :maxdepth: 1

   getting_started
   io
   examples

.. toctree::
   :caption: RNAlib
   :maxdepth: 1

   api
   vrna_3.0
   callbacks
   wrappers
   api_python
   deprecated
   linking

.. toctree::
   :caption: Miscellaneous

   faq
   contributing
   changelog
   bibliography
   license


.. doxygenpage:: index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
