import RNApath

RNApath.addSwigInterfacePath()


import RNA
import unittest

seq1 = "CGCAGGGAUACCCGCG"


def mfe_window_callback(start, end, structure, energy, data=None):
    data.append({ 'structure': structure, 'start': start, 'end' : end, 'energy' : energy})


class mfe_eval_functionTest(unittest.TestCase):

    def test_mfe_window(self):
        print "test_mfe_window\n"
        fc= RNA.fold_compound(seq1, None, RNA.OPTION_MFE | RNA.OPTION_WINDOW)
        (mfe) = fc.mfe_window()
        print "[ %6.2f" %mfe ,"]\n"
        self.assertEqual("%6.2f" % mfe, "%6.2f" % -5.60)


    def test_Lfold_cb(self):
        print "test_Lfold_cb\n"
        data = []
        mfe = RNA.Lfold_cb(seq1, 150, mfe_window_callback, data)
        self.assertTrue(len(data) == 2)
        self.assertEqual(data[0]['structure'], "(((...))).")
        self.assertEqual("%6.2f" % data[0]['energy'], "%6.2f" % -2.60)
        self.assertEqual(data[1]['structure'], "(((.(((...))))))")
        self.assertEqual("%6.2f" % data[1]['energy'], "%6.2f" % -5.60)


    def test_mfe_window_cb(self):
        print "test_mfe_window_cb\n"
        fc= RNA.fold_compound(seq1, None, RNA.OPTION_MFE | RNA.OPTION_WINDOW)
        data = []
        mfe = fc.mfe_window_cb(mfe_window_callback, data)
        self.assertTrue(len(data) == 2)
        self.assertEqual(data[0]['structure'], "(((...))).")
        self.assertEqual("%6.2f" % data[0]['energy'], "%6.2f" % -2.60)
        self.assertEqual(data[1]['structure'], "(((.(((...))))))")
        self.assertEqual("%6.2f" % data[1]['energy'], "%6.2f" % -5.60)



if __name__ == '__main__':
    unittest.main()
