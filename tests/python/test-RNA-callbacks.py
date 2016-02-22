import RNApath

RNApath.addSwigInterfacePath()

import RNA
import unittest

a = RNA.fold_compound("GGGGAAAACCCC")

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

if __name__ == '__main__':
    unittest.main();
