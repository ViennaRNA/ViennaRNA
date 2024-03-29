Multiple Sequence Alignments (MSA)
==================================

ClustalW format
---------------

The *ClustalW* format is a relatively simple text file containing a single
multiple sequence alignment of DNA, RNA, or protein sequences. It was first used
as an output format for the *clustalw* programs, but nowadays it may also be
generated by various other sequence alignment tools. The specification is straight
forward:

* The first line starts with the words::

    CLUSTAL W

  or::

    CLUSTALW

* After the above header there is at least one empty line
* Finally, one or more blocks of sequence data are following, where each
  block is separated by at least one empty line.

Each line in a blocks of sequence data consists of the sequence name followed by
the sequence symbols, separated by at least one whitespace character. Usually,
the length of a sequence in one block does not exceed 60 symbols. Optionally,
an additional whitespace separated cumulative residue count may follow the sequence
symbols. Optionally, a block may be followed by a line depicting the degree of
conservation of the respective alignment columns.

.. note::

  Sequence names and the sequences must not contain whitespace characters!
  Allowed gap symbols are the hyphen (``-``), and dot (``.``).

.. warning::

  Please note that many programs that output this format tend to truncate
  the sequence names to a limited number of characters, for instance the
  first 15 characters. This can destroy the uniqueness of identifiers in
  your MSA.

Here is an example alignment in ClustalW format:

.. include:: ../../../examples/files/alignment_clustal.aln
   :literal:


Stockholm 1.0 format
--------------------

Here is an example alignment in Stockholm 1.0 format:

.. include:: ../../../examples/files/alignment_stockholm.stk
   :literal:


.. admonition:: See also...

  :ref:`io/rna_structures:wuss notation` for legal characters and
  their interpretation in the consensus secondary structure line ``SS_cons``.


FASTA (Pearson) format
----------------------

.. note::

  Sequence names must not contain whitespace characters. Otherwise, the parts after
  the first whitespace will be dropped. The only allowed gap character is the hyphen
  (``-``).

Here is an example alignment in FASTA format:

.. include:: ../../../examples/files/alignment_fasta.fa
   :literal:


MAF format
----------

The multiple alignment format (MAF) is usually used to store multiple alignments on DNA level
between entire genomes. It consists of independent blocks of aligned sequences which are
annotated by their genomic location. Consequently, an MAF formatted MSA file may contain
multiple records. MAF files start with a line::

  ##maf

which is optionally extended by whitespace delimited key=value pairs. Lines starting with
the character (``#``) are considered comments and usually ignored.

A MAF block starts with character (``a``) at the beginning of a line, optionally followed
by whitespace delimited ``key=value`` pairs. The next lines start with character (``s``) and
contain sequence information of the form::

  s src start size strand srcSize sequence

where:

* *src* is the name of the sequence source
* *start* is the start of the aligned region within the source (0-based)
* *size* is the length of the aligned region without gap characters
* *strand* is either (``+``) or (``-``), depicting the location of the aligned
  region relative to the source
* *srcSize* is the size of the entire sequence source, e.g. the full chromosome
* *sequence* is the aligned sequence including gaps depicted by the hyphen (``-``)

Here is an example alignment in MAF format (bluntly taken from the
`UCSC Genome browser website <https://genome.ucsc.edu/FAQ/FAQformat.html#format5>`_):

.. include:: ../../../examples/files/alignment_maf.maf
   :literal:
