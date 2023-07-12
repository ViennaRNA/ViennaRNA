import unittest

if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA


motif_seq_theo = "GAUACCAG&CCCUUGGCAGC"
motif_str_theo = "(...((((&)...)))...)"
e_binding_theo = -9.22

str_hp    = "((((...((((.....(.(.....).).))))...))))...((.((((((((((.((((((.((.((((....)))).))..)))))).)))))))))))).."
e_hp      = -35.60
str_theo  = "((((...((((((((......)))))...)))...))))...((.((((((((((.((((((.((.((((....)))).))..)))))).)))))))))))).."
e_theo    = -33.82

class constraintsTest(unittest.TestCase):

    def test_sc_add_hi_motif(self):
        """Soft constraints - Motif - Hairpin"""
        fc = RNA.fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG")
        fc.sc_add_hi_motif("GAAAAAU", "(.....)", -19)
        (ss, mfe) = fc.mfe()
        print("%s [ %6.2f ]" % (ss, mfe))


    def test_theophylline_ligand_binding_interface(self):
        """Soft constraints - Motif - Theophylline aptamer"""
        fc  = RNA.fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG")
        ret = fc.sc_add_hi_motif(motif_seq_theo, motif_str_theo, e_binding_theo)
        (ss, mfe) = fc.mfe()
        print("%s [ %6.2f ]" % (ss, mfe))
        self.assertEqual(ret,1)
        self.assertEqual(ss, str_theo)
        self.assertAlmostEqual(mfe, e_theo, delta=1e-3)


if __name__ == '__main__':
    unittest.main(testRunner=taprunner.TAPTestRunner())
