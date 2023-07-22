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
