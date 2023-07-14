Installation and Configuration
==============================

A documentation on how to configure the different features of *RNAlib*, how to install
the ViennaRNA Package, and finally, how to link you own programs against *RNAlib*

.. toctree::
   :maxdepth: 1
   :caption: Contents:


Installing the ViennaRNA Package
--------------------------------

For best portability the ViennaRNA package uses the GNU `autoconf` and `automake` tools.
The instructions below are for installing the ViennaRNA package from source. However,
pre-compiled binaries for various Linux distributions, as well as for Windows users are
available from Download section of the `ViennaRNA homepage <https://www.tbi.univie.ac.at/RNA`_.

Quick-start
^^^^^^^^^^^

Usually you'll just unpack, configure and make. To do this type::

  tar -zxvf ViennaRNA-2.6.2.tar.gz
  cd ViennaRNA-2.6.2
  ./configure
  make
  sudo make install


Installation without root privileges
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you do not have root privileges on your computer, you might want to install the ViennaRNA
Package to a location where you actually have write access to. To do so, you can set the
installation prefix of the `./configure` script like so::

  ./configure --prefix=/home/username/ViennaRNA
  make install

This will install the entire ViennaRNA Package into a new directory ViennaRNA directly
into the users username home directory.

Notes for MacOS X users
^^^^^^^^^^^^^^^^^^^^^^^

Compilation
"""""""""""

Although users will find `/usr/bin/gcc` and `/usr/bin/g++` executables in their directory
tree, these programs are not at all what they pretend to be. Instead of including the GNU
programs, Apple decided to install `clang`/`llvm` in disguise. Unfortunately, the default version
of `clang`/`llvm` does not support `OpenMP` (yet), but only complains at a late stage of
the build process when this support is required. Therefore, it seems necessary to deactivate
OpenMP support by passing the option `--disable-openmp` to the `./configure` script.

Missing EXTERN.h
""""""""""""""""

Furthermore, as far as we are informed, users are discouraged to use the Perl 5 interpreter
that is shipped with Mac OS X.

Instead, one should install a more recent version from another source, e.g. ``homebrew``.
If, however, for any reason you do not want to install your own Perl 5 interpreter but use
the one from Apple, you need to specify its include path to enable building the ViennaRNA
Perl interface. Otherwise, the file ``EXTERN.h`` will be missing at compile time.
To fix this problem, you first need to find out where ``EXTERN.h`` is located::

  sudo find /Library -type f -name EXTERN.h

Then choose the one that corresponds to your default perl interpreter (find out the version
number with ``perl -v | grep version``), simply execute the following before running the
``./configure`` script, e.g.::

  export CPATH=/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/System/Library/Perl/5.18/darwin-thread-multi-2level/CORE

if your default perl is v5.18 running on MacOSX10.15. Change the paths according to your
current setup. After that, running ``./configure`` and compilation should run fine.

.. admonition:: See also...

  Related question at stackoverflow: https://stackoverflow.com/q/52682304/18609162


Universal binaries
^^^^^^^^^^^^^^^^^^

Additionally, if you intend to build the ViennaRNA such that it runs on both, `x86_64`
and the `arm64` (such as for the M1 processors in recent MacBooks), architectures, you
need to build a so-called universal binary. Note, however, that to accomplish this task,
you might need to deactivate any third-party library dependency as in most cases, only
one architecture will be available at link time. This includes the Perl 5 and Python
interfaces but also MPFR and GSL support, possibly even more. In order to compile and
link the programs, library, and scripting language interfaces of the ViennaRNA Package
for multiple architectures, we've added a new configure switch that sets up the required
changes automatically::

  ./configure --enable-universal-binary

.. note::

  With link time optimization turned on, MacOS X's default compiler (llvm/clang)
  generates an intermediary binary format that can not easily be combined into a
  multi-architecture library. Therefore, the `--enable-universal-binary` switch
  turns off link time optimization!


RNAlib configuration
--------------------

The ViennaRNA Package includes additional executable programs such as `RNAforester`,
`Kinfold`, and `Kinwalker`. Furthermore, we include several features in our C-library
that may be activated by default, or have to be explicitly turned on at configure-time.
Below we list a selection of the available configure options that affect the features
included in all executable programs, the RNAlib C-library, and the corresponding scripting
language interface(s).


Streaming SIMD Extension
^^^^^^^^^^^^^^^^^^^^^^^^

Since version 2.3.5 our sources contain code that implements a faster multibranch loop
decomposition in global MFE predictions, as used e.g. in RNAfold. This implementation
makes use of modern processors streaming SIMD extension (SSE) that provide the capability
to execute particular instructions on multiple data simultaneously (SIMD - single instruction
multiple data, thanks to W. B. Langdon for providing the modified code). Consequently, the
time required to assess the minimum of all multibranch loop decompositions is reduced up
to about one half compared to the runtime of the original implementation. This feature is
enabled by default since version 2.4.11 and a dispatcher ensures that the correct
implementation will be selected at runtime. If for any reason you want to disable this
feature at compile-time use the following configure flag::

  ./configure --disable-simd

Scripting Language Interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ViennaRNA Package comes with scripting language interfaces for Perl 5,
Python (provided by swig), that allow one to use the implemented algorithms
directly without the need of calling an executable program. The necessary
requirements are determined at configure time and particular languages may be
deactivated automatically, Building the Python 2 interface is deactivated by
default since it reached its end-of-life on January 1st, 2020. If for any reason
you still want to build that interface, you may use the ``--with-python2`` configure
option to turn it back on. You may also switch-off particular languages by
passing the ``--without-perl`` and/or ``--without-python`` configure options, e.g.::

  ./configure --without-perl --without-python

will turn-off the Perl 5 and Python 3 interfaces.

.. tip::

  Disabling the scripting language support all-together can be accomplished using the
  following switch::

    ./configure --without-swig


Cluster Analysis
^^^^^^^^^^^^^^^^

The programs ``AnalyseSeqs`` and ``AnalyseDists`` offer some cluster analysis
tools (split decomposition, statistical geometry, neighbor joining, Ward's method)
for sequences and distance data. To also build these programs add ``--with-cluster``
to your configure options.

Kinfold
^^^^^^^

The ``kinfold`` program can be used to simulate the folding dynamics of an RNA
molecule, and is compiled by default. Use the ``--without-kinfold`` option to
skip compilation and installation of Kinfold.

RNAforester
^^^^^^^^^^^

The ``RNAforester`` program is used for comparing secondary structures using tree
alignment. Similar to ``kinfold`, use the ``--without-forester`` option to skip
compilation and installation of ``RNAforester``.

Kinwalker
^^^^^^^^^

The ``kinwalker`` algorithm performs co-transcriptional folding of RNAs, starting at
a user specified structure (default: open chain) and ending at the minimum free
energy structure. Compilation and installation of this program is deactivated by
default. Use the ``--with-kinwalker`` option to enable building and installation
of ``kinwalker``.

RNAlocmin
^^^^^^^^^

The ``RNAlocmin`` program is part of the **Basin Hopping Graph Framework** and reads
secondary structures and searches for local minima by performing a gradient walk
from each of those structures. It then outputs an energetically sorted list of
local minima with their energies and number of hits to particular minimum, which
corresponds to a size of a gradient basin. Additional output consists of barrier
trees and Arhenius rates to compute various kinetic aspects. Compilation and
installation of this program is activated by default. Use the ``--without-rnalocmin``
option to disable building and installation of ``RNAlocmin``.

RNAxplorer
^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^

To increase the performance of our implementations, the ViennaRNA Package tries
to make use of the Link Time Optimization (LTO) feature of modern C-compilers.
If you are experiencing any troubles at make-time or run-time, or the configure
script for some reason detects that your compiler supports this feature although
it doesn't, you can deactivate it using the flag::

  ./configure --disable-lto


Note, that ``gcc`` before version 5 is known to produce unreliable LTO code,
especially in combination with :ref:`SIMD <start/install:streaming simd extension>`.
We therefore recommend using a more recent compiler (GCC 5 or above) or to turn off
one of the two features, LTO or SIMD optimized code.

OpenMP
^^^^^^

To enable concurrent computation of our implementations and in some cases
parallelization of the algorithms we make use of the OpenMP API. This interface
is well understood by most modern compilers. However, in some cases it might be
necessary to deactivate OpenMP support and therefore transform *RNAlib* into a
C-library that is not entirely *thread-safe*. To do so, add the following configure
option::

  ./configure --disable-openmp

POSIX threads
^^^^^^^^^^^^^

To enable concurrent computation of multiple input data in RNAfold, and for our
implementation of the concurrent unordered insert, ordered output flush data structure
:c:type:`vrna_ostream_t` we make use of POSIX threads (pthread). This should be
supported on all modern platforms and usually does not pose any problems. Unfortunately,
we use a threadpool implementation that is not compatible with Microsoft Windows yet.
Thus, POSIX thread support can not be activated for Windows builds until we have
fixed this problem. If you want to compile RNAfold and RNAlib without POSIX threads
support for any other reasons, add the following configure option::

  ./configure --disable-pthreads

SVM Z-score filter
^^^^^^^^^^^^^^^^^^

By default, ``RNALfold`` that comes with the ViennaRNA Package allows for Z-score
filtering of its predicted results using a support vector machine (SVM). However,
the library we use to implement this feature (``libsvm``) is statically linked to
our own RNAlib. If this introduces any problems for your own third-party programs
that link against RNAlib, you can safely switch off the Z-scoring implementation
using::

  ./configure --without-svm

GNU Scientific Library
^^^^^^^^^^^^^^^^^^^^^^

The program ``RNApvmin`` computes a pseudo-energy perturbation vector that aims
to minimize the discrepancy of predicted, and observed pairing probabilities. For
that purpose it implements several methods to solve the optimization problem. Many
of them are provided by the GNU Scientific Library, which is why the ``RNApvmin``
program, and the RNAlib C-library are required to be linked against ``libgsl``.
If this introduces any problems in your own third-party programs that link against
RNAlib, you can turn off a larger portion of available minimizers in ``RNApvmin``
and linking against ``libgsl`` all-together, using the switch::

  ./configure --without-gsl

Disable C11/C++11 features
^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, we use C11/C++11 features in our implementations. This mainly accounts
for unnamed unions/structs within *RNAlib*. The configure script automatically
detects whether or not your compiler understands these features. In case you are
using an older compiler, these features will be deactivated by setting a specific
pre-processor directive. If for some reason you want to deactivate C11/C++11
features despite the capabilities of your compiler, use the following configure
option::

  ./configure --disable-c11

Deprecated symbols
^^^^^^^^^^^^^^^^^^

Since version 2.2 we are in the process of transforming the API of our *RNAlib*.
Hence, several symbols are marked as *deprecated* whenever they have been
replaced by the new API. By default, deprecation warnings at compile time are
deactivated. If you want to get your terminal spammed by tons of deprecation
warnings, enable them using::

  ./configure --enable-warn-deprecated

Single precision
^^^^^^^^^^^^^^^^

Calculation of partition functions (via ``RNAfold -p``) uses double precision
floats by default, to avoid overflow errors on longer sequences. If your machine
has little memory and you don't plan to fold sequences over 1,000 bases in
length you can compile the package to do the computations in single precision by
running::

  ./configure --enable-floatpf

.. warning::

  Using this option is discouraged and not necessary on most modern computers.

Help
^^^^

For a complete list of all ``./configure`` options and important environment variables,
type::

  ./configure --help

For more general information on the build process see the *INSTALL* file.


Linking against RNAlib
----------------------

In order to use our implemented algorithms you simply need to link your program
to our *RNAlib* C-library that usually comes along with the ViennaRNA Package
installation. If you've installed the ViennaRNA Package as a pre-build binary package,
you probably need the corresponding development package, e.g. ``viennarna-devel``, or
``viennarna-dev``. The only thing that is left is to include the ViennaRNA header
files into your source code, e.g.:

.. code:: c

  #include <ViennaRNA/mfe.h>

and start using our fast and efficient algorithm implementations.

.. admonition:: See also...

  In the :ref:`examples:c` and :ref:`examples:newAPI` sections, we list a small
  set of example code that usually is a good starting point for your application.

Compiler and Linker flags
^^^^^^^^^^^^^^^^^^^^^^^^^

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
  
    Add this directive to your pre-processor/compile flags only if *RNAlib*
    was build with the ``--disable-c11`` configure option.

  .. admonition:: See also...
  
    :ref:`start/install:disable c11/c++11 features` and :c:func:`vrna_C11_features()`

* ``VRNA_WARN_DEPRECATED`` ... Enable warnings for using deprecated symbols.

  .. note::
  
    Adding this directive enables compiler warnings whenever you use symbols
    in *RNAlib* that are marked *deprecated*.

  .. admonition:: See also...
  
    :ref:`start/install:deprecatd symbols` and :ref:`deprecated`

* ``USE_FLOAT_PF`` ... Use single precision floating point operations instead
  of double precision in partition function computations.

  .. warning::
  
    Define this macro only if *RNAlib* was build with the ``--enable-floatpf``
    configure option!

  .. admonition:: See also...

    :ref:`start/install:Single precision`


For instance, you might want to add the following definition(s) to your
pre-processor/compile flags::

  -DVRNA_DISABLE_C11_FEATURES

Finally, linking against *RNAlib* is achieved by adding the following linker flag::

  -L/usr/lib -lRNA -flto -fopenmp

Again, the path to the library, ``/usr/lib``, may be omitted if this path is
searched for libraries by default. The second flag tells the linker to include
``libRNA.a``, and the remaining two flags activate :ref:`start/install:link time optimization`
and :ref:`start/install:openmp` support, respectively.

.. note::

  Depending on your linker, the last two flags may differ.

  Depending on your configure time decisions, you can drop one or both of the last flags.

  In case you've compiled *RNAlib* with LTO support (See :ref:`start/install:link time
  optimization`) and you are using a different compiler for your third-party project that
  links against our library, you may add the ``-fnolto`` flag to disable Link Time
  Optimization.

The pkg-config tool
^^^^^^^^^^^^^^^^^^^

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
