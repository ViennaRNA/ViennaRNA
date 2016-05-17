import RNApath

RNApath.addSwigInterfacePath()


import RNA
import unittest


seq1          = "CGCAGGGAUACCCGCG"
struct1       = "(((.(((...))))))"
struct1pk     = "(.(.(([....)))])"
struct2       = "(..............)"

class GeneralTests(unittest.TestCase):
    def test_pairtable(self):
        print("test_pairtable\n")
        pairTable = RNA.ptable(struct1)
        correctPairTable = (16, 16, 15, 14, 0, 13, 12, 11, 0, 0, 0, 7, 6, 5, 3, 2, 1)
        self.assertEqual(pairTable,correctPairTable)
        
        pairTable_pk = RNA.ptable_pk(struct1pk)
        correctPairTable_pk = (16, 16, 0, 14, 0, 13, 12, 15, 0, 0, 0, 0, 6, 5, 3, 7, 1)
        self.assertEqual(pairTable_pk,correctPairTable_pk)
        
        #convert pairtable back to struct1
        self.assertEqual(struct1,RNA.db_from_ptable(pairTable))
    
    def test_basePairDistance(self):
      print("test_basePairDistance\n")
      d = RNA.bp_distance("(((.(((...))))))","(((..........)))")
      self.assertEqual(d,3)
    
   # def test_plists(self):
   #   print("test_plists\n")
      #plist = RNA.plist(struct1,0.6)
      #print plist
    
    #def test_findpath(self):
    #  print("test_findpath\n")
    #  fc = RNA.fold_compound(seq1)
    #  energy = fc.path_findpath_saddle(struct1, struct2, 10)
    #  print(energy);
    #  
    #  l = fc.path_findpath(struct1,struct2,10)
    #  print l
      
if __name__ == '__main__':
    unittest.main()
