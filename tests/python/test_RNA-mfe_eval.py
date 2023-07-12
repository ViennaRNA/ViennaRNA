import unittest

if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA


seq1          = "CGCAGGGAUACCCGCG"
struct1       = "(((.(((...))))))"
seq1Dimer     = "CGCAGGGA&ACCCGCG"
struct1Dimer  = "(((.(((..))))))"
struct11      = "(((.((.....)))))"
##         1234567890
struct1_pt    = [len(struct1),16,15,14,0,13,12,11,0,0,0,7,6,5,3,2,1]


class mfe_eval_functionTest(unittest.TestCase):

    def test_mfe(self):
        """MFE prediction - single sequence"""
        fc= RNA.fold_compound(seq1)
        (ss,mfe) = fc.mfe()
        print(ss, "[ %6.2f ]" % mfe)
        self.assertEqual(ss,struct1)


    def test_mfe_Dimer(self):
        """MFE prediction - dimer"""
        fc=RNA.fold_compound(seq1Dimer)
        (ss,mfe) = fc.mfe_dimer()
        print(ss, "[ %6.2f ]" % mfe)
        self.assertEqual(ss,struct1Dimer)


    def test_eval_structure(self):
        """Structure energy evaluation - dot-bracket string"""
        fc = RNA.fold_compound(seq1)
        energy= fc.eval_structure(struct1)
        self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60)
        print("\n", struct1, "[ %6.2f ]" % energy)


    def test_eval_structure_pt(self):
        """Structure energy evaluation - pair table"""
        fc=RNA.fold_compound(seq1)
        energy= fc.eval_structure_pt(struct1_pt) /100.; #/100 for dcal

        self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60)
        print( struct1, "[ %6.2f ]" % energy)


    # Testing with filehandler and with stdout
    def test_eval_structure_verbose(self):
        """Structure energy evaluation - dot-bracket string - verbose output"""
        fc = RNA.fold_compound(seq1)
        filename= "test_RNA-mfe_eval.py.out"
        try:
            with open(filename, 'w') as f:
                print(filename ," is opened for writing")
                energy = fc.eval_structure_verbose(struct1, f)
                energy2 = fc.eval_structure_verbose(struct1)

                self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60)
                print( struct1, "[ %6.2f ]" % energy)
        except IOError:
            print("Could not open ", filename)


    def test_eval_structure_pt_verbose(self):
        """Structure energy evaluation - pair table - verbose output"""
        filename= "test_RNA-mfe_eval.py.out"
        try:
            with open(filename, 'w') as f:
                print(filename, " is opened for writing")
                fc=RNA.fold_compound(seq1)

                energy = fc.eval_structure_pt_verbose(struct1_pt, f)/100.;  # / 100 for dcal
                energy2 = fc.eval_structure_pt_verbose(struct1_pt)/100.;  # / 100 for dcal

                self.assertEqual("%6.2f" % energy, "%6.2f" % -5.6)
                print( struct1, "[ %6.2f ]" % energy)

        except IOError:
            print("Could not open ", filename)


    def test_eval_covar_structure(self):
        """Structure energy evaluation - Covariance energy contribution"""
        s1="CCCCAAAACGGG"
        s2="CCCGAAAAGGGG"
        s3="CCCCAAAAGGGG"
        ali = [s1,s2,s3]
        covarStructure = "((((....))))"

        fc = RNA.fold_compound(ali)
        pseudoEScore=fc.eval_covar_structure(covarStructure)
        print(covarStructure, "[ %6.2f ]" % pseudoEScore)
        self.assertTrue(pseudoEScore)


    def test_eval_loop_pt(self):
        """Loop energy evaluation"""
        fc= RNA.fold_compound(seq1)
        energy= fc.eval_loop_pt(6, struct1_pt) / 100.; #/100 for dcal
        print("[ %6.2f ]" % energy)
        self.assertEqual("%6.2f" % energy,"%6.2f" % -3.3)


    def test_eval_move_del(self):
        """Move energy evaluation - base pair removal"""
        fc = RNA.fold_compound(seq1)
        energy = fc.eval_move(struct1, -7, -11);  # remove basepair (7,11) ,  energy change should be 2.10
        self.assertEqual("%6.2f" % energy, "%6.2f" % 2.10)
        print("\n", struct1, " moveset (-7,-11) --> [ %6.2f ]" % energy)


    def test_eval_move_ins(self):
        """Move energy evaluation - base pair insertion"""
        fc = RNA.fold_compound(seq1)
        energy = fc.eval_move(struct11, 7, 11);  # add basepair (7,11) ,  energy change should be -2.10
        self.assertEqual("%6.2f" % energy, "%6.2f" % -2.10)
        print("\n", struct11, " moveset (7,11) --> [ %6.2f ]" % energy)


    def test_eval_move_pt_del(self):
        """Move energy evaluation - base pair removal - pair table"""
        fc = RNA.fold_compound(seq1)
        energy = fc.eval_move_pt(struct1_pt, -7, -11) / 100.;  # remove basepair (7,11) ,  energy change should be 2.10
        self.assertEqual("%6.2f" % energy, "%6.2f" % 2.10)
        print("\n", struct1, " moveset (-7,-11) --> [ %6.2f ]" % energy)


if __name__ == '__main__':
    unittest.main(testRunner=taprunner.TAPTestRunner())
