RNAlib API v3.0
===============

Introduction
------------

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

Major changes
-------------

...


Porting to the new API
----------------------

...


Examples
--------

Examples on how to use the new v3.0 API can be found in the :doc:`/examples/c`
section.
