import RNApath

RNApath.addSwigInterfacePath()


import RNA
import unittest


seq1          = "CGCAGGGAUACCCGCG"
struct1       = "(((.(((...))))))"


class GeneralTests(unittest.TestCase):
    def test_pairtable(self):
        print("test_pairtable\n")
        pairTable = RNA.ptable(struct1)
        correctPairTable = (16, 16, 15, 14, 0, 13, 12, 11, 0, 0, 0, 7, 6, 5, 3, 2, 1);
        print pairTable[0];
        self.assertEqual(pairTable,correctPairTable);
        
if __name__ == '__main__':
    unittest.main()
