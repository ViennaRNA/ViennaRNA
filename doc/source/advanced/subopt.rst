Suboptimals and Representative Structures
=========================================

Suboptimal Structures sensu Zuker
---------------------------------

The algorithm to compute optimal secondary structures that contain
a particular base pair has been published by :cite:t:`zuker:1989`
and is based on ideas for predicting structures for circular RNAs,
in particular viroids, as presented in :cite:t:`steger:1984`.


Suboptimal Structures within an Energy Band around the MFE
----------------------------------------------------------

.. code:: python

    import RNA

    sequence = "GGGGAAAACCCC"

    subopt_data = { 'counter' : 1, 'sequence' : sequence }

    # Print a subopt result as FASTA record
    def print_subopt_result(structure, energy, data):
        if not structure == None:
            print(f">subopt {data['counter']:d}")
            print(f"{data['sequence']}\n{structure} [{energy:6.2f}]")
            # increase structure counter
            data['counter'] += 1

    # Get model details with unique ML decomposition
    md = RNA.md(uniq_ML = 1)

    # Create a 'fold_compound' for our sequence
    fc = RNA.fold_compound(sequence, md)

    # Enumerate all structures 500 dacal/mol = 5 kcal/mol arround
    # the MFE and print each structure using the function above
    fc.subopt_cb(500, print_subopt_result, subopt_data);


Boltzmann Sampling
------------------

Below follow a few ``Python`` examples that demonstrate the usage of the
Boltzmann sampling implementations. To prepare the Boltzmann sampling
implementation, we first need to compute the partition function and thereby
fill the dynamic programming matrices required for the stochastic backtracking:

.. code:: python

    import RNA

    sequence = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGU \
                GGAAUCCGUACUUAUACCGCCUGUGCGGACUACUAUCCUGACCACAUAGU"

    """
    First we prepare a fold_compound object and run the usual
    steps to compute the partition function
    """

    # create model details
    md = RNA.md()

    # activate unique multibranch loop decomposition
    md.uniq_ML = 1

    # create fold compound object
    fc = RNA.fold_compound(sequence, md)

    # compute MFE
    (ss, mfe) = fc.mfe()

    # rescale Boltzmann factors according to MFE
    fc.exp_params_rescale(mfe)

    # compute partition function to fill DP matrices
    fc.pf()

    """
    Now we are ready to perform Boltzmann sampling
    """

Let us now step by step go through the different ways to create samples
of (sub)structures from the ensemble.

We begin with sub-structures of a particular length that start at the 5'
end, i.e. pos 1, of the input sequence:

.. tab-set::

    .. tab-item:: Return value based

        .. code:: python

            # 1. backtrace a single sub-structure of length 10
            print(fc.pbacktrack5(10))

            # 2. backtrace a single sub-structure of length 50
            print(fc.pbacktrack5(50))

            # 3. backtrace multiple sub-structures of length 10 at once
            for s in fc.pbacktrack5(20, 10):
                print(s)

            # 4. backtrace multiple sub-structures of length 50 at once
            for s in fc.pbacktrack5(100, 50):
                print(s)

    .. tab-item:: Callback based

        .. code:: python

            def store_structure(structure, data):
                """
                A simple callback function that stores
                a structure sample into a list
                """
                if structure:
                    data.append(structure)

            # Draw a sample of 100 sub-structures of length 50
            ss  = []
            i   = fc.pbacktrack5(100, 50, store_structure, ss)

            for s in ss:
                print(s)

Of course one can sample substructures from anywhere within the input
sequence as well:

.. tab-set::

    .. tab-item:: Return value based

        .. code:: python

            # 1. backtrace a single substructure from the subsequence [10:50]
            print(fc.pbacktrack_sub(10, 50))

            # 2. backtrace multiple substructures from the subsequence [10:50]
            for s in fc.pbacktrack_sub(100, 10, 50):
                print(s)

            # 3. backtrace multiple substructures from the subsequence [10:50] non-redundantly
            for s in fc.pbacktrack_sub(100,
                                       10,
                                       50,
                                       RNA.PBACKTRACK_NON_REDUNDANT):
                print(s)

    .. tab-item:: Callback based

        .. code:: python

            def store_structure(structure, data):
                """
                A simple callback function that stores
                a structure sample into a list
                """
                if structure:
                    data.append(structure)

            # list container to store our structures
            ss  = []

            # 1. backtrace multiple substructures from the subsequence [10:50]
            i   = fc.pbacktrack_sub(100, 10, 50, store_structure, ss)

            # 2. backtrace multiple substructures from the subsequence [10:50] non-redundantly
            i   = fc.pbacktrack_sub(100,
                                    10,
                                    50,
                                    store_structure,
                                    ss,
                                    RNA.PBACKTRACK_NON_REDUNDANT)

            for s in ss:
                print(s)


Finally, let us draw full-length structures:

.. tab-set::

    .. tab-item:: Return value based

        .. code:: python

            # 1. backtrace a single structure (full length)
            print(fc.pbacktrack())

            # 2. backtrace multiple structures at once
            for s in fc.pbacktrack(100):
                print(s)

            # 3. backtrace multiple structures non-redundantly
            for s in fc.pbacktrack(100, RNA.PBACKTRACK_NON_REDUNDANT):
                print(s)

    .. tab-item:: Callback based

        .. code:: python

            def store_structure(structure, data):
                """
                A simple callback function that stores
                a structure sample into a list
                """
                if structure:
                    data.append(structure)

            # list container to store our structures
            ss  = []

            # 1. backtrace multiple full-length structures in callback mode
            i   = fc.pbacktrack(100, store_structure, ss)

            # 2. non-redundantly backtrace multiple full-length structures in callback mode
            i   = fc.pbacktrack(100,
                                store_structure,
                                ss,
                                RNA.PBACKTRACK_NON_REDUNDANT)

            for s in ss:
                print(s)

The **non-redundant** mode for Boltzmann sampling (activated by 
:code:`RNA.PBACKTRACK_NON_REDUNDANT`) removes each sample that has been drawn
from the ensemble. Thus, no structure will be sampled more than once.

Whenever you drew a few structures in non-redundant mode and want to continue
sampling from the remaining sub-ensemble that is left, you can use the **resume**
capability for this sampling process:

.. tab-set::

    .. tab-item:: Return value based

        .. code:: python

            # perform 15 iterations and draw 500 non-redundant structures each time
            # (with resume option)
            num_samples = 500
            iterations  = 15
            d           = None # pbacktrack memory object
            s_list      = []

            for i in range(0, iterations):
                d, ss   = fc.pbacktrack(num_samples, d, RNA.PBACKTRACK_NON_REDUNDANT)
                s_list.extend(ss)

            for s in s_list:
                print(s)

    .. tab-item:: Callback based

        .. code:: python

            def store_structure(structure, data):
                """
                A simple callback function that stores
                a structure sample into a list
                """
                if structure:
                    data.append(structure)

            # Non-redundantly backtrace multiple full length structures
            # in callback mode with resume option
            num_samples = 500
            iterations  = 15
            d           = None # pbacktrack memory object
            s_list      = []

            for i in range(0, iterations):
                d, i = fc.pbacktrack(num_samples,
                                     store_structure,
                                     s_list,
                                     d,
                                     RNA.PBACKTRACK_NON_REDUNDANT)

            for s in ss:
                print(s)



Maximum Expected Accuracy (MEA) Structure
-----------------------------------------

Another useful structure representative one can determine from base pairing
probabilities :math:`p_{ij}` is the structure that exhibits the *maximum
expected accuracy* (MEA). By assuming the base pair probability is a good
measure of correctnes of a pair :math:`(i,j)`, the expected accuracy of a
structure is

.. math::

    EA(s) = \sum_{(i,j) \in s} 2 \gamma p_{ij} + \sum_{i \notin s} p^u_i

with unpaired probabilities :math:`p^u_i = 1 - \sum_{j \neq i} p_{ij}` and
weighting factor :math:`\gamma` that allows us to weight *paired* against
*unpaired* positions. The MEA algorithm attempts to find the structure that
maximizes :math:`EA`. **Higher values** of :math:`\gamma` result in more base
pairs of lower probability and thus **higher sensitivity**. **Low values** of
:math:`\gamma` result in structures containing only highly likely pairs, i.e.
**high specificity**.


Example MEA structure Prediction

.. tab-set::

    .. tab-item:: C

        .. literalinclude:: ../../../examples/mea_structure.c
            :encoding: latin-1
            :language: c
            :emphasize-lines: 33-34

    .. tab-item:: Python

        .. literalinclude:: ../../../examples/Python/mea_structure.py
            :encoding: latin-1
            :language: python
            :emphasize-lines: 22

    .. tab-item:: Perl

        .. literalinclude:: ../../../examples/Perl5/mea_structure.pl
            :encoding: latin-1
            :language: perl
            :emphasize-lines: 22

The examples above should all print the following::

    AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGA
    .....((((..(((...)))..))))((((.(....).)))).((.(((.....))).)).. ( -8.30 MEA = 47.47)



The Centroid Structure
----------------------

The centroid is the structure with the minimal average distance to all other structures

.. math::

    <d(s)> = \sum_{(i,j) \in s} (1-p_{ij}) + \sum_{(i,j) \notin s} p_{ij}

Thus, the centroid is simply the structure containing all pairs with :math:`p_ij > 0.5`


Example Centroid structure Prediction

.. tab-set::

    .. tab-item:: C

        .. literalinclude:: ../../../examples/centroid_structure.c
            :encoding: latin-1
            :language: c
            :emphasize-lines: 33-34

    .. tab-item:: Python

        .. literalinclude:: ../../../examples/Python/centroid_structure.py
            :encoding: latin-1
            :language: python
            :emphasize-lines: 22

    .. tab-item:: Perl

        .. literalinclude:: ../../../examples/Perl5/centroid_structure.pl
            :encoding: latin-1
            :language: perl
            :emphasize-lines: 22

The examples above should all print the following::

    AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGA
    .....((((..((.....))..))))((((........)))).((.(((.....))).)).. ( -7.20 d = 8.96)

