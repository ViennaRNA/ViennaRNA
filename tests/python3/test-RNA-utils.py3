import RNApath

RNApath.addSwigInterfacePath(3)


import RNA
import unittest


seq1          = "CGCAGGGAUACCCGCG"
struct1       = "(((.(((...))))))"
struct1pk     = "(.(.(([....)))])"
struct2       = "(..............)"

class GeneralTests(unittest.TestCase):
    def test_pairtable(self):
        print("test_pairtable")
        pairTable = RNA.ptable(struct1)
        correctPairTable = (16, 16, 15, 14, 0, 13, 12, 11, 0, 0, 0, 7, 6, 5, 3, 2, 1)
        self.assertEqual(pairTable,correctPairTable)

        pairTable_pk = RNA.ptable_pk(struct1pk)
        correctPairTable_pk = (16, 16, 0, 14, 0, 13, 12, 15, 0, 0, 0, 0, 6, 5, 3, 7, 1)
        self.assertEqual(pairTable_pk,correctPairTable_pk)

        #convert pairtable back to struct1
        self.assertEqual(struct1,RNA.db_from_ptable(pairTable))


    def test_basePairDistance(self):
        print("test_basePairDistance")
        d = RNA.bp_distance("(((.(((...))))))","(((..........)))")
        self.assertEqual(d,3)


    def test_basePairDistancePairTable(self):
        print("test_basePairDistancePairTable")
        pt1 = RNA.ptable("(((.(((...))))))")
        pt2 = RNA.ptable("(((..........)))")

        d = RNA.bp_distance(pt1,pt2)

        self.assertEqual(d,3)


    def test_basePairDistancePseudoknot(self):
        print("test_basePairDistancePairTable")
        d1 = RNA.bp_distance(struct1, struct1pk)
        d2 = RNA.bp_distance(struct1, struct1pk, RNA.BRACKETS_ANY)

        self.assertEqual(d1, 2)
        self.assertEqual(d2, 3)


    def test_plists(self):
        print("test_plists")
        plist = RNA.plist(struct1,0.6)
        print(plist)


    def test_filename_sanitize_simple(self):
        print("test_filename_sanitize_simple")
        fn = "bla/bla??_foo\\bar\"r<u>m:ble"
        fs = RNA.filename_sanitize(fn)
        self.assertEqual(fs, "blabla_foobarrumble")
        fs = RNA.filename_sanitize(fn, '-')
        self.assertEqual(fs, "bla-bla--_foo-bar-r-u-m-ble")


    def test_filename_sanitize_special(self):
        print("test_filename_sanitize_special_names")
        fn = "??"
        fs = RNA.filename_sanitize(fn)
        self.assertEqual(fs, "")
        fs = RNA.filename_sanitize(fn, '.')
        self.assertEqual(fs, "")


    def test_filename_sanitize_long(self):
        print("test_filename_sanitize_long_names")
        fn = "%s%s%sDEFGHIJ.svg" % ("A" * 120, "B" * 120, "C" * 10)
        fs = RNA.filename_sanitize(fn)
        self.assertEqual(fs, "%s%s%sD.svg" % ("A" * 120, "B" * 120, "C" * 10))
        fn = "A.%s%s%sDEFGHIsvg" % ("A" * 120, "B" * 120, "C" * 10)
        fs = RNA.filename_sanitize(fn)
        self.assertEqual(fs, "A.%s%s%sDEF" % ("A" * 120, "B" * 120, "C" * 10))



if __name__ == '__main__':
    unittest.main()
