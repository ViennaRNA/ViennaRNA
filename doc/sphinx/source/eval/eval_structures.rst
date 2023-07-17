Evaluation of Structures
========================

Several different functions to evaluate the free energy of a full secondary structure
under a particular set of parameters and the Nearest Neighbor Energy model are available
in *RNAlib*.

For most of them, two different forms of representations for the secondary structure may
be used:

* The Dot-Bracket string
* A pair table representation

Furthermore, the evaluation functions are divided into **basic** and **simplified** variants,
where **basic** functions require the use of a :c:type:`vrna_fold_compound_t` data structure
holding the sequence string, and model configuration (settings and parameters).

The **simplified** functions, on the other hand, provide often used default model settings
that may be called directly with only sequence and structure data.

Finally, **verbose** variants exist for some functions that allow one to print the
(individual) free energy contributions to some ``FILE`` stream.

.. doxygengroup:: eval


