Linking against RNAlib
======================

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

  In the :doc:`/examples/c` section, we list a small set of example code that
  usually is a good starting point for your application.

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
  
    :ref:`configuration:disable c11/c++11 features` and :c:func:`vrna_C11_features()`

* ``VRNA_WARN_DEPRECATED`` ... Enable warnings for using deprecated symbols.

  .. note::
  
    Adding this directive enables compiler warnings whenever you use symbols
    in *RNAlib* that are marked *deprecated*.

  .. admonition:: See also...
  
    :ref:`configuration:deprecated symbols` and :doc:`deprecated`

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
