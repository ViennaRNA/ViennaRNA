import RNApath

RNApath.addSwigInterfacePath(3)

import RNA
import unittest

a = RNA.fold_compound("GGGGAAAACCCC")

b = { 'test': "something" }
c = { 'what' : "theheck" }

def bla(d, data=None):
    if d == RNA.STATUS_MFE_PRE:
        print("about to start MFE recursions\n")
    if d == RNA.STATUS_MFE_POST:
        print("finished MFE recursions\n")
    print(data)


def blubb(i,j,k,l,d,data=None):
    if d == RNA.DECOMP_PAIR_HP:
        """
        Add -10 kcal/mol to any hairpin
        """
        return -1000

    return 0

a.add_auxdata(b, None)
a.add_callback(bla)
a.sc_add_data(c, None)
a.sc_add_f(blubb)
(s, mfe) = a.mfe()
print(s, "[%6.2f]" % mfe)


def print_subopt_result(structure, energy, data=None):
    if not structure == None:
        print(structure, "[%6.2f]" % energy)


RNA.cvar.uniq_ML = 1
a = RNA.fold_compound("GGGGAAAACCCC")
a.subopt_cb(500, print_subopt_result);


if __name__ == '__main__':
    unittest.main()
