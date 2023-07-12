import unittest


if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA


seq          = "AGGAAACCUUAAUUGGUUA"
structpk     = ".((...))(([[..))]]."
seq2         = "AAAAAAAA"
struct2      = "(......)"
struct3      = "..(..).."

class ensemble_defectTest(unittest.TestCase):
    def test_ensemble_defect_sanity(self):
        """Ensembe defect (sanity check)"""
        fc = RNA.fold_compound(seq2)
        fc.pf()

        ed = fc.ensemble_defect(struct2) * len(seq2)
        self.assertEqual(ed, 2.0)

        ed = fc.ensemble_defect(struct3) * len(seq2)
        self.assertEqual(ed, 2.0)

    def test_ensemble_defect(self):
        """Ensembe defect"""
        fc = RNA.fold_compound(seq)
        fc.pf()

        ed = fc.ensemble_defect(structpk)
        self.assertEqual(ed, 0.614080983833787)

        pt = RNA.ptable(structpk)
        ed = fc.ensemble_defect(pt)
        self.assertEqual(ed, 0.614080983833787)

    def test_ensemble_defect_pk(self):
        """Ensemble defect (pk)"""
        fc = RNA.fold_compound(seq)
        fc.pf()

        ed = fc.ensemble_defect(structpk, RNA.BRACKETS_ANY)
        self.assertEqual(ed, 0.7279184335061499)

        pt = RNA.ptable(structpk, RNA.BRACKETS_ANY)
        ed = fc.ensemble_defect(pt)
        self.assertEqual(ed, 0.7279184335061499)


if __name__ == '__main__':
    unittest.main(testRunner=taprunner.TAPTestRunner())
