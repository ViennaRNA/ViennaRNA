Command Files
=============

The RNAlib and many programs of the ViennaRNA Package can parse and apply data from
so-called *command files*. These commands may refer to structure constraints or even
extensions of the RNA folding grammar (such as :ref:`grammar:unstructured domains`).

Commands are given as a line of whitespace delimited data fields. The syntax we use
extends the constraint definitions used in the `mfold <http://mfold.rna.albany.edu/?q=mfold>`_
or `UNAfold <http://mfold.rna.albany.edu/?q=DINAMelt/software>`_ software, where
each line begins with a command character followed by a set of positions.

However, we introduce several new commands, and allow for an optional loop type context
specifier in form of a sequence of characters, and an orientation flag that enables one
to force a nucleotide to pair upstream, or downstream.

Constraint commands
-------------------

The following set of commands is recognized:

* ``F`` ... Force
* ``P`` ... Prohibit
* ``C`` ... Conflicts/Context dependency
* ``A`` ... Allow (for non-canonical pairs)
* ``E`` ... Soft constraints for unpaired position(s), or base pair(s)

RNA folding grammar exensions
-----------------------------

* ``UD`` ... Add ligand binding using the :ref:`grammar:unstructured domains` feature

Specification of the loop type context
--------------------------------------

The optional loop type context specifier ``[LOOP]`` may be a combination of the following:

* ``E`` ... Exterior loop
* ``H`` ... Hairpin loop
* ``I`` ... Internal/Interior loop
* ``M`` ... Multibranch loop
* ``A`` ... All loops

For structure constraints, we additionally allow one to address base pairs enclosed
by a particular kind of loop, which results in the specifier ``[WHERE]`` which consists
of ``[LOOP]`` plus the following character:

* ``i`` ... enclosed pair of an Interior loop
* ``m`` ... enclosed pair of a Multibranch loop

If no ``[LOOP]`` or ``[WHERE]`` flags are set, all contexts are considered
(equivalent to ``A`` ).

Controlling the orientation of base pairing
-------------------------------------------

For particular nucleotides that are forced to pair, the following ``[ORIENTATION]`` flags
may be used:

* ``U`` ... Upstream
* ``D`` ... Downstream

If no ``[ORIENTATION]`` flag is set, both directions are considered.

Sequence coordinates
--------------------

Sequence positions of nucleotides/base pairs are 1-based and consist of three
positions :math:`i`, :math:`j`, and :math:`k`. Alternativly, four positions may
be provided as a pair of two position ranges :math:`[i:j]`, and :math:`[k:l]`
using the ``-`` sign as delimiter within each range, i.e. ``i-j``, and ``k-l``.

Valid constraint commands
-------------------------

Below are resulting general cases that are considered *valid* constraints:

* **"Forcing a range of nucleotide positions to be paired"**::

    F i 0 k [WHERE] [ORIENTATION]

  Description:

  Enforces the set of :math:`k` consecutive nucleotides starting at
  position :math:`i` to be paired. The optional loop type specifier ``[WHERE]``
  allows to force them to appear as closing/enclosed pairs of certain types of
  loops.

* **"Forcing a set of consecutive base pairs to form"**:::

    F i j k [WHERE]

  Description:

  Enforces the base pairs :math:`(i,j), \ldots, (i+(k-1), j-(k-1))` to form.
  The optional loop type specifier ``[WHERE]`` allows to specify in which loop
  context the base pair must appear.

* **"Prohibiting a range of nucleotide positions to be paired"**::

    P i 0 k [WHERE]

  Description:

  Prohibit a set of :math:`k` consecutive nucleotides to participate
  in base pairing, i.e. make these positions unpaired. The optional loop type
  specifier ``[WHERE]`` allows to force the nucleotides to appear within the
  loop of specific types.

* **"Probibiting a set of consecutive base pairs to form"**::

      P i j k [WHERE]

  Description:

  Probibit the base pairs :math:`(i,j), \ldots, (i+(k-1), j-(k-1))` to form.
  The optional loop type specifier ``[WHERE]`` allows to specify the type of
  loop they are disallowed to be the closing or an enclosed pair of.

* **"Prohibiting two ranges of nucleotides to pair with each other"**::

      P i-j k-l [WHERE]

  Description:

  Prohibit any nucleotide :math:`p \in [i:j]` to pair with any other nucleotide
  :math:`q \in [k:l]`. The optional loop type specifier ``[WHERE]`` allows to
  specify the type of loop they are disallowed to be the closing or an enclosed pair of.

* **"Enforce a loop context for a range of nucleotide positions"**::

      C i 0 k [WHERE]

  Description:

  This command enforces nucleotides to be unpaired similar to *prohibiting* nucleotides to be paired,
  as described above. It too marks the corresponding nucleotides to be unpaired, however,
  the ``[WHERE]`` flag can be used to enforce specfic loop types the nucleotides must appear in.

* **"Remove pairs that conflict with a set of consecutive base pairs"**::

    C i j k

  Description:

  Remove all base pairs that conflict with a set of consecutive base pairs
  :math:`(i,j), \ldots, (i+(k-1), j-(k-1))`. Two base pairs :math:`(i,j)` and
  :math:`(p,q)` conflict with each other if :math:`i < p < j < q`, or
  :math:`p < i < q < j`.

* **"Allow a set of consecutive (non-canonical) base pairs to form"**::

    A i j k [WHERE]

  Description:

  This command enables the formation of the consecutive base pairs
  :math:`(i,j), \ldots, (i+(k-1), j-(k-1))`, no matter if they are *canonical*,
  or *non-canonical*. In contrast to the above ``F`` and ``W`` commands, which remove
  conflicting base pairs, the ``A`` command does not. Therefore, it may be used to
  allow *non-canoncial* base pair interactions. Since the RNAlib does not contain
  free energy contributions :math:`E_{ij}` for non-canonical base pairs :math:`(i,j)`,
  they are scored as the *maximum* of similar, known contributions. In terms of a
  *Nussinov* like scoring function the free energy of non-canonical base pairs is
  therefore estimated as

  .. math::

    E_{ij} = \min \left[  \max_{(i,k) \in \{GC, CG, AU, UA, GU, UG\}} E_{ik}, \max_{(k,j) \in \{GC, CG, AU, UA, GU, UG\}} E_{kj} \right].

  The optional loop type specifier ``[WHERE]`` allows to specify in which loop
  context the base pair may appear.

* **"Apply pseudo free energy to a range of unpaired nucleotide positions"**::

    E i 0 k e

  Description:

  Use this command to apply a pseudo free energy of :math:`e` to the set of :math:`k`
  consecutive nucleotides, starting at position :math:`i`. The pseudo free energy is
  applied only if these nucleotides are considered unpaired in the recursions, or
  evaluations, and is expected to be given in units of :math:`\text{kcal} \cdot \text{mol}^{-1}`.

* **"Apply pseudo free energy to a set of consecutive base pairs"**::

    E i j k e

  Description:

  Use this command to apply a pseudo free energy of :math:`e` to the set of base pairs
  :math:`(i,j), \ldots, (i+(k-1), j-(k-1))`. Energies are expected to be given in
  units of :math:`\text{kcal} \cdot \text{mol}^{-1}`.


Valid domain extensions commands
--------------------------------

* **"Add ligand binding to unpaired motif (a.k.a. unstructured domains)"**::

    UD m e [LOOP]

  Description:

  Add ligand binding to unpaired sequence motif :math:`m` (given in IUPAC format,
  capital letters) with binding energy :math:`e` in particular loop type(s).

  Example::

    UD  AAA   -5.0    A

  The above example applies a binding free energy of :math:`-5\,\text{kcal} \cdot \text{mol}^{-1}`
  for a motif ``AAA`` that may be present in all loop types.


