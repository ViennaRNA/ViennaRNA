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

# compute MEA structure
(mea_struct, mea) = fc.MEA(1.0)

# compute free energy of MEA structure
mea_en = fc.eval_structure(mea_struct)

# print only the MEA structure, its free energy
# in kcal/mol and the maximum exprected accuracy
print(f"{seq}\n{mea_struct} ({mea_en:6.2f} d={mea:.2f})")
