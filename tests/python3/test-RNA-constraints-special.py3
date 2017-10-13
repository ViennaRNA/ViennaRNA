import RNApath

RNApath.addSwigInterfacePath(3)

import RNA
import unittest

seq_con     = "CCCAAAAGGGCCCAAAAGGG"
str_con     = "..........(((....)))"
str_con_def = "(((....)))(((....)))"

datadir = RNApath.getDataDirPath()


class constraintsTest(unittest.TestCase):

    def test_sc_add_hi_motif(self):
        print("test_sc_add_hi_motif")
        fc= RNA.fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG")
        #struct =          ".(((((..((((((((((((((((((....(((((((............)))))))........)))))))))))))...)))))))))).............."
        #structWithMotif=      "((((...((((((((......)))))...)))...))))...((.((((((((((.((((((.((.((((....)))).))..)))))).)))))))))))).."
        ret = fc.sc_add_hi_motif("GAUACCAG&CCCUUGGCAGC","(...((((&)...)))...)",-9.22)
        (ss,mfe) = fc.mfe()
        print(ret,"\t",ss, "[ %6.2f ]" % mfe)
        self.assertEqual(ret,1)


    def test_theophylline_ligand_binding_interface(self):
        print("test_theophylline_ligand_binding_interface")
        RNA.noLonelyPairs = 0
        fc = RNA.fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG")
        (ss, mfe) = fc.mfe()
        print("%s [ %6.2f ]" % (ss, mfe))

        fc.sc_add_hi_motif("GAUACCAG&CCCUUGGCAGC", "(...((((&)...)))...)", -9.22)
        (ss, mfe) = fc.mfe()
        print("%s [ %6.2f ]" % (ss, mfe))

        fc.sc_remove()

        fc.sc_add_hi_motif("GAAAAAU", "(.....)", -19)
        (ss, mfe) = fc.mfe()
        print("%s [ %6.2f ]" % (ss, mfe))


if __name__ == '__main__':
    unittest.main()
