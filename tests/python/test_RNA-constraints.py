import os
import unittest

if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA


seq_con     = "CCCAAAAGGGCCCAAAAGGG"
str_con     = "..........(((....)))"
str_con_def = "(((....)))(((....)))"
seq_long    = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUC"

def mfe_window_callback(start, end, structure, energy, data=None):
    data.append({ 'structure': structure, 'start': start, 'end' : end, 'energy' : energy})


class constraintsTest(unittest.TestCase):
    DATADIR = os.environ.get('VRNA_TEST_DATA', os.sep.join(["tests", "data"]))

    def test_constraints_add_hc(self):
        """Add hard  contraints from file"""
        fc = RNA.fold_compound(seq_con)
        (ss, mfe) = fc.mfe()
        print(f"{ss} [{mfe:6.2f}]")
        self.assertEqual(ss,str_con_def)

        fc.constraints_add(os.sep.join([self.DATADIR, "hc.txt"]))
        (ss2, mfe2) = fc.mfe()
        print(f"{ss2} [{mfe2:6.2f}]")
        self.assertEqual(ss2, str_con)


    def test_constraints_add(self):
        """Add soft contraints from file"""
        fc = RNA.fold_compound(seq_con)
        (ss, mfe) = fc.mfe()
        print(f"{ss} [{mfe:6.2f}]")
        self.assertEqual(ss,str_con_def)

        fc.constraints_add(os.sep.join([self.DATADIR, "sc.txt"]))
        (ss, mfeNew) = fc.mfe()
        print(f"{ss} [{mfeNew:6.2f}]")
        self.assertEqual(f"{mfe:6.2f}", f"{(mfeNew + 5):6.2f}")


    def test_hc_add_up(self):
        """Add hard constraints - unpaired"""
        #seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        #str_con_def=    "(((....)))(((....)))"
        #str_con=    "..........(((....)))"
        fc = RNA.fold_compound(seq_con)
        fc.hc_add_up(1,RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
        (ss,mfe) = fc.mfe()
        print(ss, "[ %6.2f ]" % mfe)
        self.assertEqual(ss,".((....)).(((....)))")


    def test_hc_add_bp_nonspecific(self):
        """Add hard constraints - base pairing unspecific"""
        #GGGCCCCCCCCCCCCCCCCC
        #(((......)))........
        fc=RNA.fold_compound("GGGCCCCCCCCCCCCCCCCC")
        fc.hc_add_bp_nonspecific(20, -1, RNA.CONSTRAINT_CONTEXT_ENFORCE | RNA.CONSTRAINT_CONTEXT_ALL_LOOPS) # force the last base to pair with some bases upstream
        (ss,mfe) = fc.mfe()
        print(ss, "[ %6.2f ]" % mfe)
        self.assertEqual(ss,"(((..............)))")


    def test_hc_add_bp(self):
        """Add hard constraints - base pairs"""
        seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        str_con_def=    "(((....)))(((....)))"
        fc=RNA.fold_compound(seq_con)
        fc.hc_add_bp(1,20,RNA.CONSTRAINT_CONTEXT_ENFORCE | RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
        (ss,mfe) = fc.mfe()
        print(ss, "[ %6.2f ]" % mfe)
        self.assertEqual(ss,"(((..............)))")


    def test_hc_add_from_db(self):
        """Add hard constraints from dot-bracket string"""
        #seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        #str_con_def=    "(((....)))(((....)))"
        #hc.txt=    "xxx................."
        #str_con=    "..........(((....)))"
        fc = RNA.fold_compound(seq_con)
        fc.hc_add_from_db("xxx.................")
        (ss,mfe) = fc.mfe()
        print(ss, "[ %6.2f ]" % mfe)
        self.assertEqual(ss,str_con)


    def test_hc_mfe_window_bp(self):
        """Hard constraints - base pairs in sliding-window MFE prediction"""
        fc = RNA.fold_compound(seq_long, None, RNA.OPTION_WINDOW)
        fc.hc_add_bp(1, 10, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_ENFORCE);
        fc.hc_add_bp(101, 110, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_ENFORCE);
        data = []
        mfe = fc.mfe_window_cb(mfe_window_callback, data)
        for hit in data:
            if (hit['start'] <= 101) and (hit['end'] >= 110):
                # must contain base pair (101,110)
                pt = RNA.ptable(hit['structure'])
                self.assertTrue(pt[101 - hit['start'] + 1] == (110 - hit['start'] + 1))

            if (hit['start'] == 1) and (hit['end'] >= 10):
                # must contain base pair (101,110)
                pt = RNA.ptable(hit['structure'])
                self.assertTrue(pt[1] == 10)


if __name__ == '__main__':
    constraintsTest.DATADIR = RNApath.getDataDirPath()
    unittest.main(testRunner=taprunner.TAPTestRunner())
