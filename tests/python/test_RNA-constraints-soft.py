import unittest

if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA


seq_con     = "CCCAAAAGGGCCCAAAAGGG"
short_seq   = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUA"
seq_long    = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUC"
str_con     = "..........(((....)))"
str_con_def = "(((....)))(((....)))"

def mfe_window_callback(start, end, structure, energy, data=None):
    data.append({ 'structure': structure, 'start': start, 'end' : end, 'energy' : energy})


class constraintsTest(unittest.TestCase):

    def test_sc_set_up(self):
        """Soft constraints - unpaired nucleotides"""
        seq_sc  =      "CCCAAAAGGG"
        fc = RNA.fold_compound(seq_sc)
        (ss,mfe) = fc.mfe()
        print(ss, "[ %6.2f ]" % mfe)
        self.assertEqual(ss,"(((....)))")

        fc.sc_init()

        m= [0.0,-5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];  #E 1 0 1 -5 ,  position 1 gets -5 if unpaired , vector starts with 0 and not 1

        fc.sc_set_up(m)
        (ss,mfeNew) = fc.mfe()
        print(ss, "[ %6.2f ]" % mfeNew)
        self.assertEqual("%6.2f" %mfeNew,"%6.2f" % -5.70)


    def test_sc_set_bp(self):
        """Soft constraints - base pairs"""
        #add energy of -5 to basepair 1-9 if formed, prefed structure should now be ((.....))., with a energy of -4.90, #matrix is also beginning with position 0
        m = [[0 for x in range(11)] for y in range(11)]
        m[1][9] = -5.0; # base 1-9 should get -5.0 if basepair
        m[9][1] = -5.0; # base 1-9 should get -5.0 if basepair

        seq_sc = "CCCAAAAGGG"
        fc = RNA.fold_compound(seq_sc)
        fc.sc_set_bp(m)
        (ss,mfeNew) = fc.mfe()
        print(ss, "[ %6.2f ]" % mfeNew)
        self.assertEqual("%6.2f" %mfeNew,"%6.2f" % -4.90)


    def test_sc_shift(self):
        """
        Soft constraints - shift energy landscape

        Compute partition function and base pair probabilities both, constrained
        and unconstrained, where the constraint simply shifts the free energy base
        line by -1 kcal/mol per nucleotide.
        When comparing both results, equilibrium probabilities must not have changed,
        except for free energy of the ensemble!
        """
        fc = RNA.fold_compound(short_seq)
        # unconstrained partition function
        ss, dG = fc.pf()

        bpp = fc.bpp()

        # add constraints
        for i in range(1, len(short_seq) + 1):
          fc.sc_add_up(i, -1.0)

        for i in range(1, len(short_seq)):
            for j in range(i + 1, len(short_seq) + 1):
                fc.sc_add_bp(i, j, -2)

        # constrained partition function
        ss2, dG2 = fc.pf()
        bpp2 = fc.bpp()

        # check if ensemble free energies are actually different
        self.assertTrue(dG != dG2)

        # check nucleotide pairing propensities
        self.assertEqual(ss, ss2)

        # check pairing probabilities
        for i in range(1, len(short_seq)):
            for j in range(i + 1, len(short_seq) + 1):
                self.assertEqual("%1.8f" % bpp[i][j], "%1.8f" % bpp2[i][j])


class constraintsMFEWindowTest(unittest.TestCase):
    def test_sc_mfe_window_add_bp(self):
        """Soft constraints - base pairs sliding window MFE"""
        fc = RNA.fold_compound(seq_long, None, RNA.OPTION_WINDOW)

        # add twice -10.0 kcal/mol on base pair (55,60)
        fc.sc_add_bp(55, 60, -10.0, RNA.OPTION_WINDOW)
        fc.sc_add_bp(55, 60, -10.0, RNA.OPTION_WINDOW)

        # allow pair (55,60) to actually form (might be non-canonical)
        fc.hc_add_bp(55, 60, RNA.CONSTRAINT_CONTEXT_NO_REMOVE | RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

        data = []
        mfe = fc.mfe_window_cb(mfe_window_callback, data)

        for hit in data:
            if (hit['start'] <= 55) and (hit['end'] >= 60):
                # compose actual dot bracket string (including potential 5' dangle nucleotide
                if hit['start'] > 1:
                    d = 1
                    ss = "." + hit['structure']
                else:
                    d = 0
                    ss = hit['structure']

                # get corresponding subsequence
                s = seq_long[hit['start'] - 1 - d: hit['end']]

                # re-evaluate free energy of subsequence/hit
                e = RNA.energy_of_struct(s, ss)

                # energy difference between both must be -20.0 kcal/mol (if the constrained base pair is present)
                self.assertEqual("%6.2f" % hit['energy'], "%6.2f" % (e - 20.0))


    def test_sc_mfe_window_add_up(self):
        """Soft constraints - unpaired nucleotides sliding window MFE"""
        fc = RNA.fold_compound(seq_long, None, RNA.OPTION_WINDOW)

        # add twice -5.0 kcal/mol per unpaired nucleotide in segment [55,60]
        for i in range(50, 61):
            fc.sc_add_up(i, -5.0, RNA.OPTION_WINDOW)

        # force segment [50,60] to stay unpaired
        for i in range(50,61):
            fc.hc_add_up(i, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

        data = []
        mfe = fc.mfe_window_cb(mfe_window_callback, data)

        for hit in data:
            if (hit['start'] <= 55) and (hit['end'] >= 60):
                # count number of unpiared nucleotides with bonus
                c = 0
                for i in range(50, 61):
                    if hit['structure'][i - hit['start']] == '.':
                        c = c + 1

                # compose actual dot bracket string (including potential 5' dangle nucleotide
                if hit['start'] > 1:
                    d = 1
                    ss = "." + hit['structure']
                else:
                    d = 0
                    ss = hit['structure']

                # get corresponding subsequence
                s = seq_long[hit['start'] - 1 - d: hit['end']]

                # re-evaluate free energy of subsequence/hit
                e = RNA.energy_of_struct(s, ss)

                # energy difference between both must be c * -5.0 kcal/mol
                self.assertEqual("%6.2f" % hit['energy'], "%6.2f" % (e + (c * -5.0)))



if __name__ == '__main__':
    unittest.main(testRunner=taprunner.TAPTestRunner())
