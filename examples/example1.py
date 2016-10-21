import RNA

sequence = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCAA"

# compute minimum free energy (mfe) and corresponding structure
(structure, mfe) = RNA.fold(sequence)

# print output
print "%s\n%s [ %6.2f ]" % (structure, mfe)
