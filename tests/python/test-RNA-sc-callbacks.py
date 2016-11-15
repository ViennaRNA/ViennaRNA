import RNApath

RNApath.addSwigInterfacePath()

import RNA
import unittest

seq1 = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCC"
struct1 = ".(((..((.((((((((((....).)))).)).))))).)))(((....))).."
seq2 = "GGGGAAAACCCC"
struct2 = "((((....))))"
struct3 = "(((((..)))))"

RNA.cvar.dangles=0

c = { 'what' : "theheck" }

def MaximumMatching(i,j,k,l,d,data=None):
    e = 0

    if d == RNA.DECOMP_PAIR_HP:
        e = e - data['fc'].eval_hp_loop(i,j)
        e = e - 100

    elif d == RNA.DECOMP_PAIR_IL:
        e = e - data['fc'].eval_int_loop(i, j, k, l)
        e = e - 100

    elif d == RNA.DECOMP_PAIR_ML:
        u = j - i - 1 - (k - l + 1)
        e = e - data['params'].MLclosing
        e = e - data['params'].MLintern[0]
        e = e - u * data['params'].MLbase

        e = e - 100

    elif d == RNA.DECOMP_ML_ML_STEM:
        u = j - i + 1 - (j - l + 1) - (k - i + 1)
        e = e - data['params'].MLintern[0]
        e = e - u * data['params'].MLbase

    elif d == RNA.DECOMP_ML_STEM:
        u = j - i + 1 - (k - l + 1)
        e = e - data['params'].MLintern[0]
        e = e - u * data['params'].MLbase

    elif d == RNA.DECOMP_ML_ML:
        u = j - i + 1 - (k - l + 1)
        e = e - u * data['params'].MLbase

    elif d == RNA.DECOMP_ML_UP:
        u = j - i + 1
        e = e - u * data['params'].MLbase

    elif d == RNA.DECOMP_EXT_STEM:
        e = e - data['fc'].E_ext_loop(k, l)

    elif d == RNA.DECOMP_EXT_STEM_EXT:
        e = e - data['fc'].E_ext_loop(i, k)

    elif d == RNA.DECOMP_EXT_EXT_STEM:
        e = e - data['fc'].E_ext_loop(l, j)

    elif d == RNA.DECOMP_EXT_EXT_STEM1:
        e = e - data['fc'].E_ext_loop(l, j-1)

    return e


def blubb(i,j,k,l,d,data=None):
    if d == RNA.DECOMP_PAIR_HP:
        """
        Add -10 kcal/mol to any hairpin
        """
        return -1000

    return 0


"""
The backtracking callback must return a list of base pairs
Here, the base pairs may be given in one of the three ways
shown below:
"""

"""
1. We create a list of dictionaries with 'i' and 'j'
keys that specify the coordinates of the base pair (i,j)
"""
def bt_hp_dict(i,j,k,l,d,data=None):
    if d == RNA.DECOMP_PAIR_HP:
        return [{ 'i' : i+1, 'j' : j-1 }]

    return None


"""
2. We create a list of tuples (i,j)
"""
def bt_hp_tuples(i,j,k,l,d,data=None):
    if d == RNA.DECOMP_PAIR_HP:
        return [(i+1, j-1)]

    return None


"""
3. We create a list of RNA::basepair objects
"""
def bt_hp_basepair(i,j,k,l,d,data=None):
    if d == RNA.DECOMP_PAIR_HP:
        bp = RNA.basepair()
        bp.i = i+1
        bp.j = j-1

        return [ bp ]

    return None


class mfe_eval_functionTest(unittest.TestCase):

    def test_maximum_matching(self):
        print "Revert MFE to MaximumMatching"
        mm_data = { 'fc': RNA.fold_compound(seq1),
                    'params': RNA.param()
                  }
        fc = RNA.fold_compound(seq1)
        fc.sc_add_data(mm_data, None)
        fc.sc_add_f(MaximumMatching)
        (s, mm) = fc.mfe()
        print "%s\n%s (max BPs: %d)\n" %  (seq1, s, -mm)
        self.assertEqual(s, struct1)
        self.assertTrue(-mm == 18)

        mm_data = { 'fc': RNA.fold_compound(seq2),
                    'params': RNA.param()
                  }
        fc = RNA.fold_compound(seq2)
        fc.sc_add_data(mm_data, None)
        fc.sc_add_f(MaximumMatching)
        (s, mm) = fc.mfe()
        print "%s\n%s (max BPs: %d)\n" %  (seq2, s, -mm)
        self.assertEqual(s, struct2)
        self.assertTrue(-mm == 4)


    def test_backtrack_hp_dict(self):
        print "Hairpin backtracking (Dictionary)"
        fc = RNA.fold_compound(seq2)
        fc.sc_add_bt(bt_hp_dict)
        (s, mfe) = fc.mfe()
        print "%s %6.2f\n" %  (s, mfe)
        self.assertEqual(s, struct3)


    def test_backtrack_hp_tuples(self):
        print "Hairpin backtracking (Tuples)"
        fc = RNA.fold_compound(seq2)
        fc.sc_add_bt(bt_hp_tuples)
        (s, mfe) = fc.mfe()
        print "%s %6.2f\n" %  (s, mfe)
        self.assertEqual(s, struct3)


    def test_backtrack_hp_basepair(self):
        print "Hairpin backtracking (basepair)"
        fc = RNA.fold_compound(seq2)
        fc.sc_add_bt(bt_hp_basepair)
        (s, mfe) = fc.mfe()
        print "%s %6.2f\n" %  (s, mfe)
        self.assertEqual(s, struct3)


if __name__ == '__main__':
    unittest.main();
