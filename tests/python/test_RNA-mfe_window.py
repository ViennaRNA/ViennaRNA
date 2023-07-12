import unittest

if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA


seq1 = "CGCAGGGAUACCCGCG"
longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"
s1="CCCCAAAACGGG"
s2="CCCGAAAAGGGG"
s3="CCCCAAAAGGGG"
ali = [s1,s2,s3]


def mfe_window_callback(start, end, structure, energy, data=None):
    data.append({ 'structure': structure, 'start': start, 'end' : end, 'energy' : energy})


class mfe_window_functionTest(unittest.TestCase):

    def test_mfe_window(self):
        """MFE prediction - fold_compound.mfe_window()"""
        fc= RNA.fold_compound(seq1, None, RNA.OPTION_MFE | RNA.OPTION_WINDOW)
        (mfe) = fc.mfe_window()
        print("[ %6.2f ]" % mfe)
        self.assertEqual("%6.2f" % mfe, "%6.2f" % -5.60)


    def test_Lfold_cb(self):
        """MFE prediction - RNA.Lfold_cb"""
        data = []
        mfe = RNA.Lfold_cb(seq1, 150, mfe_window_callback, data)
        self.assertTrue(len(data) == 2)
        self.assertEqual(data[0]['structure'], "((....)).")
        self.assertEqual("%6.2f" % data[0]['energy'], "%6.2f" % -0.80)
        self.assertEqual(data[1]['structure'], "(((.(((...))))))")
        self.assertEqual("%6.2f" % data[1]['energy'], "%6.2f" % -5.60)


    def test_mfe_window_cb(self):
        """MFE prediction - fold_compound.mfe_window_cb"""
        fc= RNA.fold_compound(seq1, None, RNA.OPTION_MFE | RNA.OPTION_WINDOW)
        data = []
        mfe = fc.mfe_window_cb(mfe_window_callback, data)
        self.assertTrue(len(data) == 2)
        self.assertEqual(data[0]['structure'], "((....)).")
        self.assertEqual("%6.2f" % data[0]['energy'], "%6.2f" % -0.80)
        self.assertEqual(data[1]['structure'], "(((.(((...))))))")
        self.assertEqual("%6.2f" % data[1]['energy'], "%6.2f" % -5.60)


    def test_aliLfold_cb(self):
        """MFE prediction - sequence aignments - RNA.aliLfold_cb"""
        data = []
        mfe = RNA.aliLfold_cb(ali, 150, mfe_window_callback, data)
        self.assertTrue(len(data) == 2)
        self.assertEqual(data[0]['structure'], "(((.....)))")
        self.assertEqual("%6.2f" % data[0]['energy'], "%6.2f" % -1.30)
        self.assertEqual(data[1]['structure'], "(((......)))")
        self.assertEqual("%6.2f" % data[1]['energy'], "%6.2f" % -2.70)


    def test_mfe_window_cb(self):
        """MFE prediction - sequence alignments - fold_compound.mfe_window_cb"""
        fc= RNA.fold_compound(ali, None, RNA.OPTION_MFE | RNA.OPTION_WINDOW)
        data = []
        mfe = fc.mfe_window_cb(mfe_window_callback, data)
        self.assertTrue(len(data) == 2)
        self.assertEqual(data[0]['structure'], "(((.....)))")
        self.assertEqual("%6.2f" % data[0]['energy'], "%6.2f" % -1.30)
        self.assertEqual(data[1]['structure'], "(((......)))")
        self.assertEqual("%6.2f" % data[1]['energy'], "%6.2f" % -2.70)



if __name__ == '__main__':
    unittest.main(testRunner=taprunner.TAPTestRunner())
