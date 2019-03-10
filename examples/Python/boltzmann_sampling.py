import RNA

sequence = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUACUUAUACCGCCUGUGCGGACUACUAUCCUGACCACAUAGU"

"""
A simple callback function that stores
a structure sample into a list
"""
def store_structure(s, data):
    if s:
        data.append(s)


"""
First we prepare a fold_compound object
"""

# create model details
md = RNA.md()

# activate unique multibranch loop decomposition
md.uniq_ML = 1

# create fold compound object
fc = RNA.fold_compound(sequence, md)

# compute MFE
(ss, mfe) = fc.mfe()

# rescale Boltzmann factors according to MFE
fc.exp_params_rescale(mfe)

# compute partition function to fill DP matrices
fc.pf()

"""
Now we are ready to perform Boltzmann sampling
"""

# 1. backtrace a single sub-structure of length 10
print "%s" % fc.pbacktrack5(10)


# 2. backtrace a single sub-structure of length 50
print "%s" % fc.pbacktrack5(50)

# 3. backtrace multiple sub-structures of length 10 at once
for s in fc.pbacktrack5(20, 10):
    print "%s" % s

# 4. backtrace multiple sub-structures of length 50 at once
for s in fc.pbacktrack5(100, 50):
    print "%s" % s

# 5. backtrace a single structure (full length)
print "%s" % fc.pbacktrack()

# 6. backtrace multiple structures at once
for s in fc.pbacktrack(100):
    print "%s" % s

# 7. backtrace multiple structures non-redundantly
for s in fc.pbacktrack(100, RNA.PBACKTRACK_NON_REDUNDANT):
    print "%s" % s

# 8. backtrace multiple structures non-redundantly (with resume option)
num_samples = 500
iterations  = 15
d           = None # pbacktrack memory object
s_list      = list()

for i in range(0, iterations):
    d, ss   = fc.pbacktrack(num_samples, d, RNA.PBACKTRACK_NON_REDUNDANT)
    s_list  = s_list + list(ss)

for s in s_list:
    print "%s" % s

# 9. backtrace multiple sub-structures of length 50 in callback mode
ss  = list()
i   = fc.pbacktrack5(100, 50, store_structure, ss)

for s in ss:
    print "%s" % s

# 10. backtrace multiple full-length structures in callback mode
ss  = list()
i   = fc.pbacktrack(100, store_structure, ss)

for s in ss:
    print "%s" % s

# 11. non-redundantly backtrace multiple full-length structures in callback mode
ss  = list()
i   = fc.pbacktrack(100, store_structure, ss, RNA.PBACKTRACK_NON_REDUNDANT)

for s in ss:
    print "%s" % s

# 12. non-redundantly backtrace multiple full length structures
# in callback mode with resume option
ss = list()
d  = None # pbacktrack memory object

for i in range(0, iterations):
    d, i = fc.pbacktrack(num_samples, store_structure, ss, d, RNA.PBACKTRACK_NON_REDUNDANT)

for s in ss:
    print "%s" % s
