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
    print d
    print data
    if d == RNA.DECOMP_PAIR_HP:
        """
        add an extra basepair inside any hairpin
        """
        return [(i+1,j-1)]
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
