import RNApath

RNApath.addSwigInterfacePath(3)


import RNA
import unittest

sq1  =  "CGCAGGGAUACCCGCG";
struct1="(((.(((...))))))";
struct11 = "(((.((.....)))))";
#sq2  ="GCGCCCAUAGGGACGC";
#struct2="((((((...))).)))";
#sq4  ="GCGCACAUAGUGACGC";
#struct4="(..(((...)))...)";

class FoldCompoundTest(unittest.TestCase):
	def test_eval_structure(self):		
		fc = RNA.fold_compound(sq1);
		energy= fc.eval_structure(struct1);
		self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60);
		print ("\n", struct1, "%6.2f" % energy);
		
	def test_eval_structure_verbose(self):
		fc = RNA.fold_compound(sq1);
		f=open("output_test.txt",'w');
		energy = fc.eval_structure_verbose(struct1,f);
		self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60);
		print ("\n", struct1, "%6.2f" % energy);
		f.close();
	def test_eval_move_del(self):
		fc = RNA.fold_compound(sq1);
		energy = fc.eval_move(struct1,-7,-11);  # remove basepair (6,11) ,  energy change should be 2.10
		self.assertEqual("%6.2f" % energy, "%6.2f" % 2.10);
		print ("\n", struct1, " moveset (-7,-11) --> %6.2f" % energy);
		
	def test_eval_move_ins(self):
		fc = RNA.fold_compound(sq1);
		energy = fc.eval_move(struct11,7,11);  # add basepair (6,11) ,  energy change should be -2.10
		self.assertEqual("%6.2f" % energy, "%6.2f" % -2.10);
		print ("\n", struct11, " moveset (7,11) --> %6.2f" % energy);	
		
		
		
if __name__ == '__main__':
	unittest.main();
