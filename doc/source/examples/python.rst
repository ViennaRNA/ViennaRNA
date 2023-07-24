Python Examples
===============

MFE Prediction (flat interface)
-------------------------------

.. literalinclude:: ../../../examples/Python/helloworld_mfe.py
   :encoding: latin-1

MFE Prediction (object oriented interface)
------------------------------------------

.. literalinclude:: ../../../examples/Python/oo_example1.py
   :encoding: latin-1

Suboptimal Structure Prediction
-------------------------------

.. literalinclude:: ../../../examples/Python/subopt.py
   :encoding: latin-1

Boltzmann Sampling
------------------

a.k.a. Probabilistic Backtracing

.. literalinclude:: ../../../examples/Python/boltzmann_sampling.py
   :encoding: latin-1

RNAfold -p MEA equivalent
-------------------------

.. literalinclude:: ../../../examples/Python/RNAfold_MEA.py
   :encoding: latin-1

MFE Consensus Structure Prediction
----------------------------------

.. literalinclude:: ../../../examples/Python/helloworld_mfe_comparative.py
   :encoding: latin-1

MFE Prediction (deviating from default settings)
------------------------------------------------

.. literalinclude:: ../../../examples/Python/helloworld_nondefault.py
   :encoding: latin-1

Fun with Soft Constraints
-------------------------

.. literalinclude:: ../../../examples/Python/maximum_matching.py
   :encoding: latin-1

Parsing Alignments
------------------

Reading the first entry from a STOCKHOLM 1.0 formatted MSA file ``msa.stk`` may look like this:

.. code:: python

  num, names, aln, id, ss = RNA.file_msa_read("msa.stk")

Similarly, if the file contains more than one alignment, one can use the :py:func:`RNA.file_msa_read_record()`
function to subsequently read each alignment separately:

.. code:: python

  with open("msa.stk") as f:
      while True:
          num, names, aln, id, ss = RNA.file_msa_read_record(f)
          if num < 0:
              break
          elif num == 0:
              print("empty alignment")
          else:
              print(names, aln)

After successfully reading the first record, the variable ``num`` contains the number of
sequences in the alignment (the actual return value of the C-function), while the variables
``names``, ``aln``, ``id``, and ``ss`` are lists of the sequence names and aligned sequences,
as well as strings holding the alignment ID and the structure as stated in the ``SS_cons`` line,
respectively.

.. note::

  The last two return values may be empty strings in case the alignment does not provide
  the required data.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

