#!/usr/bin/python
#

import RNA

seq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGA"

# create fold_compound data structure
# (required for all subsequently applied  algorithms)
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

# print only the centroid structure, its free energy
# in kcal/mol and the distance to the ensemble
print(f"{seq}\n{centroid_struct} ({centroid_en:6.2f} d={dist:.2f})")
