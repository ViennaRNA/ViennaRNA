#!/usr/bin/python
#

import RNA

seq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUA"

# create fold_compound data structure (required for all subsequently applied  algorithms)
fc = RNA.fold_compound(seq)

# compute MFE and MFE structure
(mfe_struct, mfe) = fc.mfe()

# rescale Boltzmann factors for partition function computation
fc.exp_params_rescale(mfe)

# compute partition function
(pp, pf) = fc.pf()

# compute centroid structure
(centroid_struct, dist) = fc.centroid()

# compute free energy of centroid structure
centroid_en = fc.eval_structure(centroid_struct)

# compute MEA structure
(MEA_struct, MEA) = fc.MEA()

# compute free energy of MEA structure
MEA_en = fc.eval_structure(MEA_struct)

# print everything like RNAfold -p --MEA
print("{}\n{} ({:6.2f})".format(seq, mfe_struct, mfe))
print("{} [{:6.2f}]".format(pp, pf))
print("{} {{{:6.2f} d={:.2f}}}".format(centroid_struct, centroid_en, dist))
print("{} {{{:6.2f} MEA={:.2f}}}".format(MEA_struct, MEA_en, MEA))
print(" frequency of mfe structure in ensemble {:g}; ensemble diversity {:-6.2f}".format(fc.pr_structure(mfe_struct), fc.mean_bp_distance()))
