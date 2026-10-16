Combinatorics Algorithms
========================

For some applications, e.g. :doc:`/advanced/rip`, algorithms from the
field of **combinatorics** are essential. Within our library, we implement
a handful of such algorithms.

Enumerating Necklaces
---------------------

One of such algorithms is the so-called **fast algorithm to generate necklaces with
fixed content** by :cite:t:`sawada:2003`. Given a list of entity counts (the elements
of the necklace) it efficiently enumerates all non-cyclic permutations of the content,
i.e. the complete list of necklaces. This means each necklace consists of the sum of
entity counts of the input list.

Here is a small ``Python`` example that demonstrates the usage of our implementation:

.. code:: python-console

    >>> import RNA
    ... 
    ... # 4 different elements, each of count 1
    ... for n in RNA.enumerate_necklaces((1, 1, 1, 1)):
    ...     print(n)
    ...     
    (0, 3, 2, 1)
    (0, 3, 1, 2)
    (0, 2, 3, 1)
    (0, 2, 1, 3)
    (0, 1, 3, 2)
    (0, 1, 2, 3)

The result can be interpreted as follows: the first necklace consists of entity 0,
entity 3, entity 2, and entity 1, in that order. The second consists of entities
0, 3, 1, and 2, and so on.


Rotational Symmetry
-------------------

Another problem that appears in the context of **circular RNAs** as well as
RNA-RNA interaction is the detection of rotational symmetry. We implement functions
to detect and compute the order of rotational symmetry for sequences of objects
and secondary structures.

A small ``Python`` example to use our function(s) may look as follows:

.. code:: python-console

    >>> import RNA
    ... 
    ... # a sequence of unsigned integer objects with rotational symmetry of 2
    ... print(len(RNA.rotational_symmetry([0, 1, 1, 0, 1, 1])))
    ... 
    2

This also works for ``strings``

.. code:: python-console

    >>> print(len(RNA.rotational_symmetry("ABBABBABB")))
    ... 
    3

and :ref:`dot-bracket-notation` for secondary structures:

.. code:: python-console

    >>> fc = RNA.fold_compound("GCGCGCGCGCGC", RNA.md(circ=1))
    ... print(len(fc.rotational_symmetry_db("(((...)))...")))
    ... 
    2

The latter example shows a circular RNA sequence with rotational symmetry
of 2. Here, shifting the structure by 6 nucleotides we would yield another
structure that is indistinguishable from the first one.


Multichoose
-----------

Sometimes, it is also useful to enumerate all combinations of subsets numbers of
a particular range, e.g. *n choose k* or *n multichoose k*. The latter is implemented
as part of our library as it serves well for enumerating complexes of RNA strands
in the context of RNA-RNA interactions. In particular, it compiles a list of
**k-combinations**, i.e. a list of multisets of size ``k`` from a set of integer
values in the range from ``0`` to ``n - ``. Here is a small ``Python`` snippet
demonstrating its usage:

.. code:: python-console

    >>> import RNA
    ... 
    ... print(RNA.n_multichoose_k(3, 2))
    ... 
    ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


Boustrophedon Scheme
--------------------

Interestingly, the complexity of some problems that loop of a sequence of numbers
are dependent on the order in which this sequence is traversed. Here, alternating
the direction of processing from the left and right can turn out beneficial. We
use this strategy in :ref:`advanced-boltzmann-sampling` and more details can be
found in :cite:t:`ponty:2008`.

In ``Python`` you can obtain a list of alternating numbers to loop over like this:

.. code:: python-console

    >>> import RNA
    ... 
    ... for i in RNA.boustrophedon(11, 20)[1:]:
    ...     print(i)
    ...
    11
    20
    12
    19
    13
    18
    14
    17
    15
    16

Note here, that we sliced-off the 0th element from the list that is returned by
:c:func:`vrna_boustrophedon`, since the first element simply stores the length
of the list.

Combinatorics API
-----------------

.. seealso::

    See the :doc:`/api/utils/combinatorics` section for a full overview of our
    implemented functions.
