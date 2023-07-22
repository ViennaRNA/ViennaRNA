Energy Parameters
=================

Modified Bases
--------------

The functions :c:func:`vrna_sc_mod()`, :c:func:`vrna_sc_mod_json()` and alike
implement an energy correction framework to account for modified bases in the
secondary structure predictions. To supply these functions with the energy
parameters and general specifications of the base modification, the following
``JSON`` data format may be used:

JSON data must consist of a header section ``modified_bases`` This header
is an object with the mandatory keys:

* ``name`` specifying a name of the modified base
* ``unmodified`` that consists of a single upper-case letter of the unmodified
  version of this base,
* the ``one_letter_code`` key to specify which letter is used for the modified
  bases in the subsequent energy parameters, and
* an array of `pairing_partners``

The latter must be uppercase characters. An optional ``sources`` key may contain
an array of related publications, e.g. those the parameters have been derived from.

Next to the header may follow additional keys to specify the actual energy
contributions of the modified base in various loop contexts. All energy
contributions must be specified in free energies :math:`\Delta G` in units of
:math:`\text{kcal} \cdot \text{mol}^{-1}`. To allow for rescaling of the free
energies at temperatures that differ from the default (:math:`37^\circ C`),
enthalpy parameters :math:`\Delta H` may be specified as well. Those, however
are optional. The keys for free energy (at :math:`37^\circ C`) and enthalpy
parameters have the suffixes ``_energies`` and ``_enthalpies``, respectively.

The parser and underlying framework currently supports the following
loop contexts:

* **base pair stacks** (via the ``stacking`` key prefix).

  This key must point to an object with one key value pair for each
  stacking interaction data is provided for. Here, the key consists
  of four upper-case characters denoting the interacting bases, where
  the the first two represent one strand in 5' to 3' direction and the
  last two the opposite strand in 3' to 5' direction. The values
  are energies in :math:`kcal \cdot mol^{-1}`.
* **terminal mismatches** (via the ``mismatch`` key prefix).

  This key points to an object with key value pairs for each mismatch
  energy parameter that is available. Keys are 4 characters long
  nucleotide one-letter codes as used in base pair stacks above.
  The second and fourth character denote the two unpaired mismatching
  bases, while the other two represent the closing base pair.

* **dangling ends** (via the ``dangle5`` and ``dangle3`` key prefixes).

  The object behind these keys, again, consists of key value pairs
  for each dangling end energy parameter. Keys are 3 characters long
  where the first two represent the two nucleotides that form the
  base pair, and the third is the unpaired base that either stacks
  on the 3' or 5' end of the enclosed part of the base pair.

* **terminal pairs** (via the ``terminal`` key prefix).

  Terminal base pairs, such as AU or GU, sometimes receive an
  additional energy penalty. The object behind this key may
  list energy parameters to apply whenever particular base
  pairs occur at the end of a helix. Each of those parameters
  is specified as key value pair, where the key consists of
  two upper-case characters denoting the terminal base pair.

Below is a JSON template specifying most of the possible input
parameters. Actual energy parameter files can be found in the
source code tarball within the ``misc/`` subdirectory.

.. literalinclude:: ../../../misc/rna_mod_template_parameters.json
   :language: json

An actual example of real-world data may look like

.. literalinclude:: ../../../misc/rna_mod_pseudouridine_parameters.json
   :language: json


