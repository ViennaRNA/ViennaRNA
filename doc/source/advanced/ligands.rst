Ligands Binding to RNA Structures
=================================

In our library, we provide two different ways to incorporate binding of small molecules
and proteins to specific RNA structures. The first approach is implemented as an actual
extension of the folding grammar. It adds auxiliary derivation rules for each case when
consecutive unpaired nucleotides are evaluated. Therefore, this model is applicable to
ligand binding to any loop context.

The second approach, on the other hand, uses the soft-constraints feature to change the
energy evaluation of hairpin- or internal-loops. Hence, it can only be appleid when a
ligand binds to a hairpin-like, or internal-loop like motif.


Unstructured Domains
--------------------

Sometime, certain ligands, like single strand binding (SSB) proteins, compete with intramolecular
base pairing of the RNA. In situations, where the dissociation constant of the ligand is known and
the ligand binds to a consecutive stretch of single-stranded nucleotides we can use the :ref:`advanced-domains-up`
functionality to extend the RNA folding grammar. This module provides a convenience default implementation
that covers most of the application scenarios.

The function :c:func:`vrna_ud_add_motif()` attaches a ligands sequence motif and corresponding
binding free energy to the list of known ligand motifs within the ``domains_up`` attribute of
:c:type:`vrna_fold_compound_t`. The first call to this function initializes the :doc:`/grammar/domains_up`
feature with our default implementation. Subsequent calls of secondary structure prediction
algorithms with the modified :c:type:`vrna_fold_compound_t` then directly include the competition
of the ligand with regules base pairing. Since we utilize the unstructured domain extension,
The ligand binding model can be removed again using the :c:func:`vrna_ud_remove()` function.


Hairpin-/Internal Loop Motifs
-----------------------------

Ligand binding to specific hairpin/internal loop like motifs using
the :ref:`advanced-constraints-soft` feature.

