import os
import unittest
from struct import *
import locale

locale.setlocale(locale.LC_ALL, 'C')

if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA


seq1      = "CGCAGGGAUACCCGCG"
struct1   = "(((.(((...))))))"

struct11  = "(((.((.....)))))"
seq2      = "GCGCCCAUAGGGACGC"
struct2   = "((((((...))).)))"
seq3      = "GCGCACAUAGUGACGC"
struct3   = "(..(((...)))...)"
align     = [seq1,seq2,seq3]
(struct, mfe) = RNA.fold(seq1)


class GeneralTests(unittest.TestCase):
    """Some general tests"""
    VERSION = os.environ.get('VRNA_VERSION', '')

    def test_version(self):
        """Version number"""
        if self.VERSION != '':
            self.assertEqual(RNA.__version__, self.VERSION)


    def test_hammingDistance(self):
        """Compute hamming distance between sequences and structures"""
        self.assertEqual(RNA.hamming(seq1,seq2),16)
        self.assertEqual(RNA.bp_distance(struct1,struct2),6)


    def test_temperature(self):
        """Default global Temperature"""
        self.assertEqual(RNA.cvar.temperature,37)


    def test_foldASequence(self):
        """Simple MFE prediction, structure, energy"""
        # new better interface
        (structure, mfe) = RNA.fold(seq1)
        self.assertEqual(structure,struct1)
        # check energy
        self.assertTrue(abs(RNA.energy_of_struct(seq1, struct1) - mfe) < 0.0001)


    def test_constrained_folding(self):
        """Simple constrained MFE folding"""
        RNA.cvar.fold_constrained = 1
        (structure,cmfe) = RNA.fold(seq1,"....xx....xx....")
        self.assertEqual(structure,'(((..........)))')
        self.assertTrue(abs(RNA.energy_of_struct(seq1, structure) - cmfe) < 0.0001)
        RNA.cvar.fold_constrained = 0


    def test_tree_distance(self):
        """Tree edit distance"""
        xstruc = RNA.expand_Full(struct1)
        T1 = RNA.make_tree(xstruc)
        xstruc = RNA.expand_Full(struct2)
        T2 = RNA.make_tree(xstruc)
        RNA.edit_backtrack = 1
        tree_dist = RNA.tree_edit_distance(T1, T2)
        # print RNA::get_aligned_line(0), RNA::get_aligned_line(1),"\n"
        self.assertEqual(tree_dist,2)


    def test_cofold_andMore(self):
        """Cofolding of two sequences, MFE and Partition function, base pairing probabilities"""
        RNA.cvar.cut_point = len(seq1)+1
        (costruct,comfe) = RNA.cofold(seq1 + seq2)
        self.assertEqual(costruct,'(((.(((...))))))((((((...))).)))')
        cmfe = RNA.energy_of_struct(seq1 + seq2, costruct)
        self.assertTrue(abs(comfe - cmfe) < 1e-5)
        (x,ac,bc,fcab,cf) = RNA.co_pf_fold(seq1 + seq2)
        self.assertTrue(cf<comfe)
        self.assertTrue(comfe-cf<1.3)

        (x,usel1, usel2, fcaa, usel3)= RNA.co_pf_fold(seq1 + seq1)
        RNA.cvar.cut_point = len(seq2)+1
        (x,usel1, usel2, fcbb, usel3)= RNA.co_pf_fold(seq2 + seq2)
        (AB,AA,BB,A,B) = RNA.get_concentrations(fcab, fcaa, fcbb,ac, bc, 1e-5, 1e-5)
        AB/=2e-5
        AA/=2e-5
        BB/=2e-5
        A/=2e-5
        B/=2e-5

        ret = (abs(AB-0.0)+abs(AA-0.00578)+abs(BB-0.01100)+abs(A-0.48843)+abs(B-0.47801))
        self.assertTrue(ret<0.0001)
        RNA.cvar.cut_point=-1
        #pf_fo ld
        RNA.cvar.do_backtrack = 1
        s,f = RNA.pf_fold(seq1)
        self.assertTrue(f < mfe)
        self.assertTrue(mfe-f < 0.8)

        p1 = RNA.get_pr(2,15)
        ii = RNA.intP_getitem(RNA.cvar.iindx, 2)

        if(RNA.pf_float_precision() != 0) :
             p2 = RNA.floatP_getitem(RNA.cvar.pr, ii-15)
        else:
             p2 = RNA.doubleP_getitem(RNA.cvar.pr, ii-15)
        self.assertTrue(p1 < 0.999)
        self.assertTrue(abs(p1-p2) < 1.2e-7)


    def test_parse_structure(self):
        """Structure parsing and splitting into components"""
        RNA.parse_structure(struct1)
        self.assertEqual(RNA.cvar.loops,2)
        self.assertEqual(RNA.cvar.pairs,6)
        self.assertEqual(RNA.cvar.unpaired,4)
        self.assertEqual(RNA.intP_getitem(RNA.cvar.loop_degree,1),2)


    def test_rna_plots(self):
        """Creating structure layout and base pair probability EPS plots"""
        RNA.PS_rna_plot(seq1, struct1, "test_ss.ps")
        anote = "2 15 1 gmark\n" + "3 cmark\n"
        RNA.PS_rna_plot_a(seq1, struct1, "test_ss_a.ps", None, anote)
        RNA.PS_dot_plot(seq1, "test_dp.ps")
        RNA.ssv_rna_plot(seq1, struct, "test.coord")
        print("please check the two postscript files test_ss.ps and test_dp.ps")
        RNA.write_parameter_file("test.par")


    def test_different_symbol_set(self):
        """Reduced symbol set"""
        RNA.cvar.symbolset = "GC"
        start = RNA.random_string(len(struct1), "GC")
        cost = 100
        trial = 0
        max_tries = 100
        while cost != 0 and trial < max_tries:
            (sinv, cost) = RNA.inverse_fold(start, struct1)
            trial += 1
        (ss, en) = RNA.fold(sinv)
        self.assertEqual(ss, struct1)
        RNA.cvar.symbolset = "ACGU"


    def test_eos_dimer(self):
        """Energy of dimer structure"""
        RNA.cvar.cut_point = 3
        e =  RNA.energy_of_struct("GCGC", "(())")
        RNA.cvar.cut_point = -1
        self.assertEqual(int(e*100+0.5), 70)


    def test_duplexfold(self):
        """Duplex fold structure prediction"""
        duplex = RNA.duplexfold(seq1, seq2)
        self.assertEqual(duplex.structure, ".(((.....(((.&)))))).")


    def test_alifold(self):
        """Comparative MFE prediction (alifold)"""
        align = ["GCCAUCCGAGGGAAAGGUU", "GAUCGACAGCGUCU-AUCG", "CCGUCUUUAUGAGUCCGGC"]

        (css, cen) = RNA.alifold(align)
        self.assertEqual(css,"(((.(((...)))..))).")
        self.assertEqual(RNA.aln_consensus_mis(align), "SMBHBHYDRBGDVWmVKBB")
        RNA.free_alifold_arrays()


    def test_moveSets(self):
        """Move set - gradient descent (energy, structure)"""
        RNA.cvar.cut_point=-1
        struct1_move = "(..............)"
        #move_standard( sequence, start structure, move_type(GRADIENT, FIRST, ADAPTIVE), verbosity, shifts, noLP)
        (s,energy) = RNA.move_standard(seq1, struct1_move, 0, 1, 0, 0)
        print("energy = ",energy," s = ", s);
        self.assertEqual(s, "................")

    def test_moveSets2(self):
        """Move set - First (structure)"""
        struct1_move = "(..............)"
        (s,energy) =  RNA.move_standard(seq1, struct1_move, 1, 1, 0, 0)
        self.assertEqual(s, "(((.((....)).)))")


    def test_simplexycoordinates(self):
        """Simple XY coordinates"""
        coords = RNA.simple_xy_coordinates(struct1)

        for c in (coords):
            print(c.X, ",", c.Y)


    def test_model_details_defaults(self):
        """Model details data structure - default values"""
        # check model details structure
        md = RNA.md();
        self.assertEqual(int(md.dangles), 2)
        self.assertEqual(md.temperature, 37.0)


    def test_model_details_glob(self):
        """Model details data structure - default values from global settings"""
        RNA.cvar.dangles     = 0
        RNA.cvar.temperature = 40.1
        md = RNA.md()
        self.assertEqual(int(md.dangles), 0)
        self.assertEqual(int(RNA.cvar.dangles), 0)
        self.assertEqual(md.temperature, 40.1)

        #reset globals to default
        RNA.cvar.dangles = 2
        RNA.cvar.temperature= 37.0

    def test_params(self):
        """Energy Parameter structure"""
        # check parameter structures
        md = RNA.md()
        md.temperature = 45
        params = RNA.param()

        self.assertEqual(params.temperature,37.0)
        params = RNA.param(md)
        self.assertEqual(params.temperature,45)

    def test_exp_params(self):
        """Energy Parameter structure - Boltzmann factors"""
        md = RNA.md()
        md.temperature = 42.1
        pf_params = RNA.exp_param()
        self.assertEqual(pf_params.temperature,37.0)
        pf_params = RNA.exp_param(md)
        self.assertEqual(pf_params.temperature,42.1)


class FoldCompoundTest(unittest.TestCase):
    def test_create_fold_compound_Single(self):
        """Creating fold_compound - Single sequence"""
        fc = RNA.fold_compound(seq1)
        self.assertEqual(fc.type, RNA.FC_TYPE_SINGLE)


    def test_create_fold_compound_Align(self):
        """Creating fold_compound - Multiple sequence alignment"""
        fc= RNA.fold_compound(align)
        self.assertEqual(fc.type, RNA.FC_TYPE_COMPARATIVE)


    def test_create_fold_compound_2D(self):
        """Creating fold_compound - Single sequence 2D fold"""
        fc= RNA.fold_compound(seq1,seq2,seq3)
        self.assertTrue(fc)


    ###centroid.h
    def test_centroid(self):
        """fold_compound method - centroid structure"""
        fc=RNA.fold_compound(align)
        fc.pf()
        (sc,dist) = fc.centroid()
        print(sc, "\tDistance of :  %6.2f" % dist)
        self.assertTrue(sc and dist)


    ## partition function from here
    def test_pf(self):
        """fold_compound method - Partition function, mean base pair distance"""
        fc= RNA.fold_compound(seq1)
        (ss,gfe) = fc.pf()
        print(ss, "[ %6.2f ]" % gfe)
        self.assertTrue(ss)
        bp_dis = fc.mean_bp_distance()
        print(seq1, "\t meanBPDistance : ", bp_dis)
        self.assertTrue(bp_dis)


    def test_pf_dimer(self):
        """fold_compound method - Partition function for dimer"""
        fc = RNA.fold_compound(seq1 + "&" + seq2)
        (costruct, comfe) = fc.mfe_dimer()
        self.assertEqual(costruct, "(((.(((...))))))((((((...))).)))")
        cmfe = fc.eval_structure(costruct)
        self.assertTrue(abs(comfe-cmfe) < 1e-5)

        (x,ac,bc,fcab,cf) = fc.pf_dimer()
        self.assertTrue((cf < comfe) and (comfe - cf < 1.3))


if __name__ == '__main__':
    GeneralTests.VERSION = RNApath.VERSION_NUMBER
    unittest.main(testRunner=taprunner.TAPTestRunner())

