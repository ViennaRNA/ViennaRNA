import RNApath

RNApath.addSwigInterfacePath()


import RNA
import unittest

seq1      = "CGCAGGGAUACCCGCG"
struct1   = "(((.(((...))))))"

struct11  = "(((.((.....)))))"
seq2      = "GCGCCCAUAGGGACGC"
struct2   = "((((((...))).)))"
seq3      = "GCGCACAUAGUGACGC"
struct3   = "(..(((...)))...)"
align     = [seq1,seq2,seq3]





class FoldCompoundTest(unittest.TestCase):

    def test_create_fold_compound_Single(self):
        print "test_create_fold_compound_Single\n"

        fc = RNA.fold_compound(seq1)
        self.assertEqual(fc.type(),0)

    def test_create_fold_compound_Align(self):
        print "test_create_fold_compound_Align\n"
        fc= RNA.fold_compound(align)
        self.assertEqual(fc.type(),1)

    def test_create_fold_compound_2D(self):
        print "test_create_fold_compound_2D\n"
        fc= RNA.fold_compound(seq1,seq2,seq3)
        self.assertTrue(fc)


    ###centroid.h
    def test_centroid(self):
        print "test_centroid\n"
        fc=RNA.fold_compound(align)
        fc.pf()
        (sc,dist) = fc.centroid()
        print  sc,"\tDistance of :  %6.2f" %dist ,"\n"
        self.assertTrue(sc and dist)


    ## partition function from here
    def test_pf(self):
        print "test_pf"
        fc= RNA.fold_compound(seq1)
        (ss,gfe) = fc.pf()
        print ss, "[ %6.2f" %gfe ,"]\n"
        self.assertTrue(ss)
        bp_dis = fc.mean_bp_distance()
        print seq1 ,"\t meanBPDistance : ", bp_dis,"\n"
        self.assertTrue(bp_dis)


    # hairpin_loops.h from here
    def test_eval_hp_loop(self):
        print "test_eval_hp_loop"
        seq1  =      "GCAAAAGG"
        struct1=    ".(....)."

        fc=RNA.fold_compound(seq1)
        #ehair = fc.eval_hp_loop(2,7)
        ehair = fc.E_hp_loop(2,7)
        print seq1, " 2,7  = [ %6.2f" %ehair ,"] \n"
        self.assertEqual("%6.2f" %ehair,"%6.2f" % +410)

        #Exterior loop evalution not working
        #eExt = fc.E_ext_hp_loop(0,0)
        #print seq1, " 2,6  = [ %6.2f" %eExt ,"] \n"
        #self.assertEqual("%6.2f" %eExt,"%6.2f" % -140)

        #length = 8
        #External loop                           :  -140
        #Hairpin  loop (  2,  7) CG              :   410
        #GCAAAAGG
        #.(....).
        #energy =   2.70


    #def test_exp_E_hp_loop(self):
        #print "test_exp_E_hp_loop"
        #seq1  =      "GCAAAAGG"
        #struct1=    ".(....)."

        #fc=RNA.fold_compound(seq1)
        #fc.pf()
        #ehair = fc.exp_E_hp_loop(2,7)
        #print seq1, " 2,7  = [ %6.2f" %ehair ,"] \n"
        #self.assertEqual("%6.2f" %ehair,"%6.2f" % +410)


    #interior_loops.h from here
    def test_E_int_loop(self):
        print "test_E_int_loop"
        #    "123456789012"
        seq1 =  "AGACAAAAGACA"
        struct1=".(.(....).)."
        fc=RNA.fold_compound(seq1,None,RNA.OPTION_MFE)
        e = fc.E_int_loop(2,11)
        print seq1, " 2,7  = [ %6.2f" %e ,"] \n"
        self.assertEqual("%6.2f" %e,"%6.2f" % +80)


if __name__ == '__main__':
    unittest.main()



