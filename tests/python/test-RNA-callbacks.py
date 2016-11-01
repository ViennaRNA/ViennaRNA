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
(s, mfe) = a.mfe()
print "%s %6.2f\n" %  (s, mfe)


a.sc_remove()

def print_subopt_result(structure, energy, data=None):
    if not structure == None:
        print "%s [%6.2f]" % (structure, energy)


RNA.cvar.uniq_ML = 1
a = RNA.fold_compound("GGGGAAAACCCC")
a.subopt_cb(500, print_subopt_result);


if __name__ == '__main__':
    unittest.main();
