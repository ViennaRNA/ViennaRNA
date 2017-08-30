import RNApath

RNApath.addSwigInterfacePath()

import RNA
import unittest

seq_con     = "CCCAAAAGGGCCCAAAAGGG"
str_con     = "..........(((....)))"
str_con_def = "(((....)))(((....)))"
seq_long    = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUC"

datadir = RNApath.getDataDirPath()

def mfe_window_callback(start, end, structure, energy, data=None):
    data.append({ 'structure': structure, 'start': start, 'end' : end, 'energy' : energy})


class constraintsTest(unittest.TestCase):

    def test_constraints_add(self):
        #seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        #str_con_def=    "(((....)))(((....)))"
        #hc.txt=    "P 1 0 2"
        #str_con=    "..........(((....)))"

        hc_file = datadir + "hc.txt"
        print "test_constraints_add"
        fc = RNA.fold_compound(seq_con)
        fc.constraints_add(hc_file)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f ]" % mfe
        self.assertEqual(ss,str_con)

        fc.hc_init()
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f ]" % mfe
        self.assertEqual(ss,str_con_def)

        #sc.txt = E 3 8 1 -5
        sc_file = datadir + "sc.txt"
        fc.sc_init()
        fc.constraints_add(sc_file)
        (ss,mfeNew) = fc.mfe()
        print ss, "[ %6.2f ]" % mfeNew
        self.assertEqual("%6.2f" %mfe, "%6.2f" % (mfeNew +5))


    def test_hc_add_up(self):
        #seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        #str_con_def=    "(((....)))(((....)))"
        #str_con=    "..........(((....)))"
        print "test_hc_add_up"
        fc = RNA.fold_compound(seq_con)
        fc.hc_add_up(1,RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f ]" % mfe
        self.assertEqual(ss,".((....)).(((....)))")


    def test_hc_add_bp_nonspecific(self):
        print "test_hc_add_bp_nonspecific"
        #GGGCCCCCCCCCCCCCCCCC
        #(((......)))........
        fc=RNA.fold_compound("GGGCCCCCCCCCCCCCCCCC")
        fc.hc_add_bp_nonspecific(20,-1); # force the last base to pair with some bases upstream
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f ]" % mfe
        self.assertEqual(ss,"(((..............)))")


    def test_hc_add_bp(self):
        print "test_hc_add_bp"
        seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        str_con_def=    "(((....)))(((....)))"
        fc=RNA.fold_compound(seq_con)
        fc.hc_add_bp(1,20,RNA.CONSTRAINT_CONTEXT_ENFORCE | RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f ]" % mfe
        self.assertEqual(ss,"(((..............)))")


    def test_hc_add_from_db(self):
        print "test_hc_add_from_db"
        #seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        #str_con_def=    "(((....)))(((....)))"
        #hc.txt=    "xxx................."
        #str_con=    "..........(((....)))"
        fc = RNA.fold_compound(seq_con)
        fc.hc_add_from_db("xxx.................")
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f ]" % mfe
        self.assertEqual(ss,str_con)


    def test_hc_mfe_window_bp(self):
        print "test test_hc_mfe_window_bp"
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
    unittest.main()
