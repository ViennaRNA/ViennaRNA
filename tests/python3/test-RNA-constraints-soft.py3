import RNApath

RNApath.addSwigInterfacePath(3)

import RNA
import unittest

seq_con     = "CCCAAAAGGGCCCAAAAGGG"
str_con     = "..........(((....)))"
str_con_def = "(((....)))(((....)))"

datadir = RNApath.getDataDirPath()


class constraintsTest(unittest.TestCase):

    def test_sc_set_up(self):
        print("test_sc_set_up")

        #        "1234567890
        seq_sc  =      "CCCAAAAGGG"
        fc = RNA.fold_compound(seq_sc)
        (ss,mfe) = fc.mfe()
        print(ss, "[ %6.2f" %mfe ,"]\n")
        self.assertEqual(ss,"(((....)))")

        fc.sc_init()

        m= [0.0,-5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];  #E 1 0 1 -5 ,  position 1 gets -5 if unpaired , vector starts with 0 and not 1

        fc.sc_set_up(m)
        (ss,mfeNew) = fc.mfe()
        print(ss, "[ %6.2f" %mfeNew ,"]\n")
        self.assertEqual("%6.2f" %mfeNew,"%6.2f" % -5.70)

    def test_sc_set_bp(self):
        print("test_sc_set_bp")

        #add energy of -5 to basepair 1-9 if formed, prefed structure should now be ((.....))., with a energy of -4.90, #matrix is also beginning with position 0
        m = [[0 for x in range(11)] for y in range(11)]
        m[1][9] = -5.0; # base 1-9 should get -5.0 if basepair
        m[9][1] = -5.0; # base 1-9 should get -5.0 if basepair

        seq_sc = "CCCAAAAGGG"
        fc = RNA.fold_compound(seq_sc)
        fc.sc_set_bp(m)
        (ss,mfeNew) = fc.mfe()
        print(ss, "[ %6.4f" %mfeNew ,"]\n")
        self.assertEqual("%6.2f" %mfeNew,"%6.2f" % -4.90)



if __name__ == '__main__':
    unittest.main()
