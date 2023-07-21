Configuration
=============

The ViennaRNA Package includes additional executable programs such as

* ``RNAforester``,
* ``Kinfold``,
* ``Kinwalker``,
* ``RNAlocmin``, and
* ``RNAxplorer``.

Furthermore, we include several features in our C-library that may be
activated by default, or have to be explicitly turned on at configure-time.
Below we list a selection of the available configure options that affect
the features included in all executable programs, the *RNAlib* C-library,
and the corresponding scripting language interface(s).


Streaming SIMD Extension
------------------------

Since version 2.3.5 our sources contain code that implements a faster
multibranch loop decomposition in global MFE predictions, as used e.g.
in ``RNAfold``. This implementation makes use of modern processors
*streaming SIMD extension* (SSE) that provide the capability to execute
particular instructions on multiple data simultaneously (*SIMD - single
instruction multiple data*, thanks to W. B. Langdon for providing the
modified code). Consequently, the time required to assess the minimum of
all multibranch loop decompositions is reduced up to about one half
compared to the runtime of the original implementation. This feature is
enabled by default since version 2.4.11 and a dispatcher ensures that
the correct implementation will be selected at runtime. If for any
reason you want to disable this feature at compile-time use the following:

.. code:: bash

  ./configure --disable-simd

Scripting Language Interfaces
-----------------------------

The ViennaRNA Package comes with scripting language interfaces for ``Perl 5``,
``Python`` (provided by `SWIG <https://www.swig.org/>`_), that allow one to
use the implemented algorithms directly without the need of calling an
executable program. The necessary requirements are determined at
configure-time and particular languages may be deactivated automatically
if the requirements are not met.

.. note::

  Building the Python 2 interface is deactivated by default since it reached
  its end-of-life on January 1st, 2020. If for any reason you still want to
  build that interface, you may use the ``--with-python2`` configure option
  to turn it back on.

You may also switch-off particular languages by passing the
``--without-perl`` and/or ``--without-python`` configure options, e.g.:

.. code:: bash

  ./configure --without-perl --without-python

will turn-off the ``Perl 5`` and ``Python 3`` interfaces.

.. tip::

  Disabling the scripting language support all-together can be accomplished using the
  following switch:

  .. code:: bash

    ./configure --without-swig


Cluster Analysis
----------------

The programs ``AnalyseSeqs`` and ``AnalyseDists`` offer some cluster analysis
tools (split decomposition, statistical geometry, neighbor joining, Ward's method)
for sequences and distance data. To also build these programs add ``--with-cluster``
to your configure options.

Kinfold
-------

The ``kinfold`` program can be used to simulate the folding dynamics of an RNA
molecule, and is compiled by default. Use the ``--without-kinfold`` option to
skip compilation and installation of Kinfold.

RNAforester
-----------

The ``RNAforester`` program is used for comparing secondary structures using tree
alignment. Similar to ``kinfold`, use the ``--without-forester`` option to skip
compilation and installation of ``RNAforester``.

Kinwalker
---------

The ``kinwalker`` algorithm performs co-transcriptional folding of RNAs, starting at
a user specified structure (default: open chain) and ending at the minimum free
energy structure. Compilation and installation of this program is deactivated by
default. Use the ``--with-kinwalker`` option to enable building and installation
of ``kinwalker``.

RNAlocmin
---------

The ``RNAlocmin`` program is part of the **Basin Hopping Graph Framework** and reads
secondary structures and searches for local minima by performing a gradient walk
from each of those structures. It then outputs an energetically sorted list of
local minima with their energies and number of hits to particular minimum, which
corresponds to a size of a gradient basin. Additional output consists of barrier
trees and Arhenius rates to compute various kinetic aspects. Compilation and
installation of this program is activated by default. Use the ``--without-rnalocmin``
option to disable building and installation of ``RNAlocmin``.

RNAxplorer
----------

The ``RNAxplorer`` is a multitool, that offers different methods to explore RNA
energy landscapes. In default mode it takes an RNA sequence as input and produces
a sample of RNA secondary structures. The repellant sampling heuristic used in
default mode iteratively penalizes base pairs of local minima of structures that
have been seen too often. This results in a diverse sample set with the most
important low free energy structures.
Compilation and installation of this program is activated by default. Note, that
this tool depends on the ``LAPACK`` library. Use the ``--without-rnaxplorer``
option to disable building and installation of ``RNAxplorer``.


Link Time Optimization
----------------------

To increase the performance of our implementations, the ViennaRNA Package tries
to make use of the *Link Time Optimization* (LTO) feature of modern C-compilers.
If you are experiencing any troubles at make-time or run-time, or the configure
script for some reason detects that your compiler supports this feature although
it doesn't, you can deactivate it using the flag:

.. code:: bash

  ./configure --disable-lto


Note, that ``gcc`` before version 5 is known to produce unreliable LTO code,
especially in combination with :ref:`SIMD <configuration:streaming simd extension>`.
We therefore recommend using a more recent compiler (GCC 5 or above) or to turn off
one of the two features, LTO or SIMD optimized code.

OpenMP
------

To enable concurrent computation of our implementations and in some cases
parallelization of the algorithms we make use of the `OpenMP <https://www.openmp.org/>`_
API. This interface is well understood by most modern compilers. However, in
some cases it might be necessary to deactivate OpenMP support and therefore transform
*RNAlib* into a C-library that is not entirely *thread-safe*. To do so, add the
following configure option:

.. code:: bash

  ./configure --disable-openmp

POSIX threads
-------------

To enable concurrent computation of multiple input data in ``RNAfold``, and for our
implementation of the concurrent unordered insert, ordered output flush data structure
:c:type:`vrna_ostream_t` we make use of POSIX threads (pthread). This should be
supported on all modern platforms and usually does not pose any problems. Unfortunately,
we use a threadpool implementation that is not compatible with Microsoft Windows yet.
Thus, POSIX thread support can not be activated for Windows builds until we have
fixed this problem. If you want to compile ``RNAfold`` and *RNAlib* without POSIX threads
support for any other reasons, add the following configure option:

.. code:: bash

  ./configure --disable-pthreads

SVM Z-score filter
------------------

By default, ``RNALfold`` that comes with the ViennaRNA Package allows for Z-score
filtering of its predicted results using a *Support Vector Machine* (SVM) provided
by the `LIBSVM <https://www.csie.ntu.edu.tw/~cjlin/libsvm/>`_ library. However, this
library is statically linked to our own *RNAlib*. If this introduces any problems for
your own third-party programs that link against *RNAlib*, you can safely switch off
the Z-scoring implementation using:

.. code:: bash

  ./configure --without-svm

GNU Scientific Library
----------------------

The program ``RNApvmin`` computes a pseudo-energy perturbation vector that aims
to minimize the discrepancy of predicted, and observed pairing probabilities. For
that purpose it implements several methods to solve the optimization problem. Many
of them are provided by the `GNU Scientific Library (GSL) <https://www.gnu.org/software/gsl/>`_,
which is why the ``RNApvmin`` program, and the *RNAlib* C-library are required to
be linked against ``libgsl``. If this introduces any problems in your own third-party
programs that link against *RNAlib*, you can turn off a larger portion of available
minimizers in ``RNApvmin`` and linking against ``libgsl`` all-together, using:

.. code:: bash

  ./configure --without-gsl

Multiple-precision Floating-Point Computations
----------------------------------------------

Our :doc:`Non-redundant Boltzmann Sampling </subopt/stochbt>` implementation uses
multi-precision floating-point computations provided by the
`GNU MPFR library <https://www.mpfr.org/>`_ by default. This requires linking against
``libmpfr`` and ``libgmp``. You can switch off this feature using:

.. code:: bash

  ./configure --disable-mpfr


Universal binaries
------------------

If you intend to build the ViennaRNA for Mac OS X such that it runs *on both*, ``x86_64``
and the ``arm64`` (Apple Silicon Processors) architectures, you need to build a
so-called *universal binary*. Note, however, that to accomplish this task, you might
need to deactivate any third-party library dependency as in most cases, only one
architecture will be available at link time. This includes the ``Perl 5`` and ``Python``
interfaces but might also concern also :ref:`MPFR <configuration:multiple-precision floating-point computations>`
and :ref:`GSL support <configuration:gnu scientific library>`, possibly even more.
In order to compile and link the programs, library, and scripting language interfaces
of the ViennaRNA Package for multiple architectures, we've added a new configure switch
that sets up the required changes automatically:

.. code:: bash

  ./configure --enable-universal-binary

.. note::

  With link time optimization turned on, MacOS X's default compiler (``llvm``/``clang``)
  generates an intermediary binary format that can not easily be combined into a
  multi-architecture library. Therefore, the ``--enable-universal-binary`` switch
  turns off :ref:`configuration:link time optimization`!


Disable C11/C++11 features
--------------------------

By default, we use C11/C++11 features in our implementations. This mainly accounts
for unnamed unions/structs within *RNAlib*. The configure script automatically
detects whether or not your compiler understands these features. In case you are
using an older compiler, these features will be deactivated by setting a specific
pre-processor directive. If for some reason you want to deactivate C11/C++11
features despite the capabilities of your compiler, use the following configure
option::

  ./configure --disable-c11

Deprecated symbols
------------------

Since version 2.2 we are in the process of transforming the API of our *RNAlib*.
Hence, several symbols are marked as *deprecated* whenever they have been
replaced by the new API. By default, deprecation warnings at compile time are
deactivated. If you want to get your terminal spammed by tons of deprecation
warnings, enable them using::

  ./configure --enable-warn-deprecated

Single precision
----------------

Calculation of partition functions (via ``RNAfold -p``) uses double precision
floats by default, to avoid overflow errors on longer sequences. If your machine
has little memory and you don't plan to fold sequences over 1,000 bases in
length you can compile the package to do the computations in single precision by
running::

  ./configure --enable-floatpf

.. warning::

  Using this option is discouraged and not necessary on most modern computers.

Help
----

For a complete list of all ``./configure`` options and important environment variables,
type::

  ./configure --help

For more general information on the build process see the *INSTALL* file.


