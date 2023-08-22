Frequently Asked Questions
==========================

Missing EXTERN.h
----------------

When compiling from source at the Mac OS X platform, users often encounter
an error message stating that the file ``EXTERN.h`` is missing for compilation
of the ``Perl 5`` wrapper. This is a known problem and due to the fact that
users are discouraged to use the ``Perl 5`` interpreter that is shipped with Mac OS X.

Instead, one should install a more recent version from another source, e.g. ``homebrew``.
If, however, for any reason you do not want to install your own ``Perl 5`` interpreter
but use the one from Apple, you need to specify its include path to enable building the
ViennaRNA Perl interface. Otherwise, the file ``EXTERN.h`` will be missing at compile time.
To fix this problem, you first need to find out where ``EXTERN.h`` is located:

.. code:: bash

  sudo find /Library -type f -name EXTERN.h

Then choose the one that corresponds to your default perl interpreter (find out the version
number with ``perl -v | grep version``), simply execute the following before running the
``./configure`` script, e.g.::

  export CPATH=/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/System/Library/Perl/5.18/darwin-thread-multi-2level/CORE

if your default perl is v5.18 running on MacOSX10.15. Change the paths according to your
current setup. After that, running ``./configure`` and compilation should run fine.

.. admonition:: See also...

  Related question at stackoverflow: https://stackoverflow.com/q/52682304/18609162


Linking fails with LTO error
----------------------------

By default, *RNAlib* is compiled with :ref:`configuration:link time optimization`.
This may introduce problems upon linking a third-party program that was either compiled
with a different compiler or compiler version. As a work-around solution, we include the
``-fno-lto`` linker flag in the output of the :ref:`pkg-config <linking:the pkg-config tool>`.
This tells the linker to not perform link time optimization even though LTO code is
included in the library. Usually, this should not affect the runtime of the algorithms
too much.

However, some linkers may not support the ``-fno-lto`` flag and fail at the linker stage.
In addition, if *RNAlib* has been compiled with ``clang``, it may not include the non-LTO
code required for linking without LTO. To resolve this issue, you may need to deactivate
link time optimization while building *RNAlib*.

.. admonition:: See also...

   :ref:`configuration:link time optimization`

