import RNApath

RNApath.addSwigInterfacePath(3)


import RNA
import unittest


class CombinatoricsTests(unittest.TestCase):
    def test_enumerate_necklaces1(self):
        print("test_enumerate_necklaces (4 sequences, 1 strand per sequences)\n")
        content   = (1, 1, 1, 1)
        perm_real = ( (0, 3, 2, 1),
                      (0, 3, 1, 2),
                      (0, 2, 3, 1),
                      (0, 2, 1, 3),
                      (0, 1, 3, 2),
                      (0, 1, 2, 3))
        permutations = RNA.enumerate_necklaces(content)
        self.assertEqual(set(perm_real), set(permutations))


    def test_enumerate_necklaces2(self):
        print("test_enumerate_necklaces (3 sequences, different number of strands per sequence)\n")
        content = (3, 1, 2)
        perm_real = ( (1, 0, 0, 0, 2, 2),
                      (1, 0, 0, 2, 0, 2),
                      (1, 0, 0, 2, 2, 0),
                      (1, 0, 2, 0, 0, 2),
                      (1, 0, 2, 0, 2, 0),
                      (1, 0, 2, 2, 0, 0),
                      (1, 2, 0, 0, 0, 2),
                      (1, 2, 0, 0, 2, 0),
                      (1, 2, 0, 2, 0, 0),
                      (1, 2, 2, 0, 0, 0))
        permutations = RNA.enumerate_necklaces(content)
        self.assertEqual(set(perm_real), set(permutations))


if __name__ == '__main__':
    unittest.main()
