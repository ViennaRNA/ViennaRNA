import RNApath

RNApath.addSwigInterfacePath(3)


import RNA
import unittest


seq          = "AGGAAACCUUAAUUGGUUA"
structpk     = ".((...))(([[..))]]."

class ensemble_defectTest(unittest.TestCase):
    def test_ensemble_defect(self):
        print("test_ensemble_defect")
        fc = RNA.fold_compound(seq)
        fc.pf()

        ed = fc.ensemble_defect(structpk)
        self.assertEqual(ed, 0.6140797258673892)

        pt = RNA.ptable(structpk)
        ed = fc.ensemble_defect(pt)
        self.assertEqual(ed, 0.6140797258673892)

        ed = fc.ensemble_defect(structpk, RNA.BRACKETS_ANY)
        self.assertEqual(ed, 0.7279171755397522)

        pt = RNA.ptable(structpk, RNA.BRACKETS_ANY)
        ed = fc.ensemble_defect(pt)
        self.assertEqual(ed, 0.7279171755397522)


if __name__ == '__main__':
    unittest.main()
