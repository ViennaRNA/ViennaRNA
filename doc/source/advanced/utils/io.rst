Files and I/O
=============

Nucleic Acid Sequences and Structures
-------------------------------------


Multiple Sequence Alignments
----------------------------

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



Command Files
-------------


Other I/O Utils
---------------

