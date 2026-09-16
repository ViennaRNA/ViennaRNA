RNA Secondary Structures
========================

Functions to create, parse, convert, manipulate, and compare secondary
structure representations.

Dot-Bracket Notation
--------------------


WUSS Notation
-------------

Washington University Secondary Structure (WUSS) notation
related functions


Pair Tables
-----------

Pair Table Representation of Secondary Structures

.. code:: python-console

    >>> import RNA
    ... # create a pair table from dot-bracket structure
    ... pt = RNA.ptable("(((.(((...))))))")
    ... # print each base pair
    ... for i, v in enumerate(pt[1:], start = 1):
    ...     if pt[i] > i:
    ...         print(f"base pair ({i},{pt[i]})")
    ...         
    base pair (1,16)
    base pair (2,15)
    base pair (3,14)
    base pair (5,13)
    base pair (6,12)
    base pair (7,11)


Pair Lists
----------

Pair List Representation of Secondary Structures


Abstract Shapes
---------------

Abstract Shapes Representation of Secondary Structures


Helix Lists
-----------

Helix List Representation of Secondary Structures


Tree Representation
-------------------

Tree Representation of Secondary Structures


Distance measures
-----------------

Distance measures between Secondary Structures


Benchmarking
------------

Benchmarking Secondary Structure Predictions


Other Structure Utilities
-------------------------


