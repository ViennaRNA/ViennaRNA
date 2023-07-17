Ligands Binding to RNA Structures
=================================

In our library, we provide two different ways to incorporate binding of small molecules
and proteins to specific RNA structures:

- :doc:`/ligands/unstructured_domains`, and
- :doc:`/ligands/constraints`

The first approach is implemented as an actual extension of the folding grammar. It adds
auxiliary derivation rules for each case when consecutive unpaired nucleotides are evaluated.
Therefore, this model is applicable to ligand binding to any loop context.

The second approach, on the other hand, uses the soft-constraints feature to change the energy
evaluation of hairpin- or interior-loops. Hence, it can only be appleid when a ligand binds to
a hairpin-like, or interior-loop like motif.

.. toctree::
   :maxdepth: 1
   :caption: Specialized Modules:

   ligands/unstructured_domains
   ligands/constraints
