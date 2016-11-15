import RNApath

RNApath.addSwigInterfacePath()

import RNA
import unittest

sequence = "CGCAGGGAUACCCGCG"

def print_subopt_result(structure, energy, data=None):
    if not structure == None:
        print "%s [%6.2f]" % (structure, energy)


class GeneralTests(unittest.TestCase):
    def test_subopt0(self):
        print "test_subopt (pointer mode)\n"

        RNA.cvar.subopt_sorted = 1
        RNA.cvar.noLonelyPairs = 1
        solution = RNA.subopt(sequence, None, 500, None)

        print "%d suboptimals" % solution.size()

        for x in range(0,solution.size()):
        # the last access should produce a "value out of range" warning
            if solution.get(x).structure :
                print "%s %6.2f" % (solution.get(x).structure,solution.get(x).energy)
        RNA.cvar.noLonelyPairs = 0


    def test_subopt1(self):
        print "test_subopt (list mode)\n"
        ## test native array output of subopt()
        RNA.cvar.subopt_sorted = 1
        RNA.cvar.noLonelyPairs = 1
        solution = RNA.subopt(sequence, 500)
        print "%d suboptimals" % len(solution)
        for s in solution:
            print "%s %6.2f" % (s.structure,s.energy)
        RNA.cvar.noLonelyPairs = 0


    def test_subopt2(self):
        print "test_subopt (as fold_compound method)\n"

        a = RNA.fold_compound(sequence)
        solution = a.subopt(500)
        print "%d suboptimals" % len(solution)
        for s in solution:
            print "%s %6.2f" % (s.structure,s.energy)


    def test_subopt3(self):
        print "test_subopt_cb (as fold_compound method)\n"
        a = RNA.fold_compound(sequence)
        a.subopt_cb(500, print_subopt_result);


    def test_zuker_subopt(self):
        print "test Zuker subopt (pointer mode)"
        solution = RNA.zukersubopt(sequence)
        for x in range(0,solution.size()):
        # the last access should produce a "value out of range" warning
            if(solution.get(x).structure) :
                print "%s %6.2f" % (solution.get(x).structure, solution.get(x).energy)


    def test_subopt_zuker(self):
        print "testing subopt_zuker() method"
        fc = RNA.fold_compound(sequence)
        solution = fc.subopt_zuker()
        print sequence
        for s in solution:
            print "%s [ %6.2f ]" % (s.structure, s.energy)


if __name__ == '__main__':
    unittest.main();
