import RNA

# The RNA sequence
seq = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA"

# create a new model details structure with
# non-default settings
md = RNA.md(
        circ        = 1,    # Circular RNA sequence
        temperature = 25.0, # 25 Deg Celcius
        dangles     = 1     # Dangle Model 1
      )

# create a fold compound
fc = RNA.fold_compound(seq, md)

# predict Minmum Free Energy and corresponding secondary structure
(ss, mfe) = fc.mfe()

# print sequence, structure and MFE
print("{}\n{} [ {:6.2f} ]".format(seq, ss, mfe))
