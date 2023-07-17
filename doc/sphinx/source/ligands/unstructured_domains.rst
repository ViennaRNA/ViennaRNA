Ligands Binding to Unstructured Domains
=======================================

Add ligand binding to loop regions using the :doc:`/grammar/domains_up` feature.

Sometime, certain ligands, like single strand binding (SSB) proteins, compete with intramolecular
base pairing of the RNA. In situations, where the dissociation constant of the ligand is known and
the ligand binds to a consecutive stretch of single-stranded nucleotides we can use the @ref domains_up
functionality to extend the RNA folding grammar. This module provides a convenience default implementation
that covers most of the application scenarios.

The function :c:func:`vrna_ud_add_motif()` attaches a ligands sequence motif and corresponding
binding free energy to the list of known ligand motifs within the ``domains_up`` attribute of
:c:type:`vrna_fold_compound_t`. The first call to this function initializes the :doc:`/grammar/domains_up`
feature with our default implementation. Subsequent calls of secondary structure prediction
algorithms with the modified :c:type:`vrna_fold_compound_t` then directly include the competition
of the ligand with regules base pairing. Since we utilize the unstructured domain extension,
The ligand binding model can be removed again using the :c:func:`vrna_ud_remove()` function.

.. doxygengroup:: ligands_up
