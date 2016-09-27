import RNApath

RNApath.addSwigInterfacePath()

import RNA
import unittest

sequence = "GGGGAAAACCCC"
#sequence = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCC"

RNA.cvar.dangles=0

a = RNA.fold_compound(sequence)
aa = RNA.fold_compound(sequence)

mm_data = { 'fold_compound': aa,
            'params': RNA.param()
          }


b = { 'test': "something" }
c = { 'what' : "theheck" }

def bla(d, data=None):
    if d == RNA.STATUS_MFE_PRE:
        print "about to start MFE recursions\n"
    if d == RNA.STATUS_MFE_POST:
        print "finished MFE recursions\n"
    print data


def blubb(i,j,k,l,d,data=None):
    if d == RNA.DECOMP_PAIR_HP:
        """
        Add -10 kcal/mol to any hairpin
        """
        return -1000

    return 0

def MaximumMatching(i,j,k,l,d,data=None):
    e = 0

    if d == RNA.DECOMP_PAIR_HP:
        e = e - data['fold_compound'].eval_hp_loop(i,j)
        e = e - 100

    elif d == RNA.DECOMP_PAIR_IL:
        e = e - data['fold_compound'].eval_int_loop(i, j, k, l)
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
        e = e - data['fold_compound'].E_ext_loop(k, l)

    elif d == RNA.DECOMP_EXT_STEM_EXT:
        e = e - data['fold_compound'].E_ext_loop(i, k)

    elif d == RNA.DECOMP_EXT_EXT_STEM:
        e = e - data['fold_compound'].E_ext_loop(l, j)

    elif d == RNA.DECOMP_EXT_EXT_STEM1:
        e = e - data['fold_compound'].E_ext_loop(l, j-1)

    return e


def bt(i,j,k,l,d,data=None):
    """
    The backtracking callback must return a list of base pairs
    Here, the base pairs may be given in one of the three ways
    shown below:
    """

    if d == RNA.DECOMP_PAIR_HP:
        """
        1. We create a list of dictionaries with 'i' and 'j'
        keys that specify the coordinates of the base pair (i,j)
        """
        bp = { 'i' : i+1, 'j' : j-1 }

        """
        2. We create a list of tuples (i,j)
        """
        bp = (i+1, j-1)

        """
        3. We create a list of RNA::basepair objects
        """
        bp = RNA.basepair()
        bp.i = i+1
        bp.j = j-1

        return [ bp ]

    return None

a.add_auxdata(b, None)
a.add_callback(bla)
a.sc_add_data(c, None)
a.sc_add_f(blubb)
a.sc_add_bt(bt)
(s, mfe) = a.mfe()
print "%s %6.2f\n" %  (s, mfe)


a.sc_remove()

a.sc_add_data(mm_data, None)
a.sc_add_f(MaximumMatching)
(s, mfe) = a.mfe()
print "MM"
print "%s %6.2f\n" %  (s, mfe)

def print_subopt_result(structure, energy, data=None):
    if not structure == None:
        print "%s [%6.2f]" % (structure, energy)


RNA.cvar.uniq_ML = 1
a = RNA.fold_compound("GGGGAAAACCCC")
a.subopt_cb(500, print_subopt_result);


if __name__ == '__main__':
    unittest.main();
