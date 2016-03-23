

import RNA
import unittest

seq1  =  	"CGCAGGGAUACCCGCG";
struct1=	"(((.(((...))))))";
seq1Dimer = 	"CGCAGGGA&ACCCGCG";
struct1Dimer=	"(((.(((..))))))";

struct11 = 	"(((.((.....)))))";
seq2  =		"GCGCCCAUAGGGACGC";
struct2=	"((((((...))).)))";
seq3  =		"GCGCACAUAGUGACGC";
struct3=	"(..(((...)))...)";
align=[seq1,seq2,seq3];

filename="output_test.txt";

##		 1234567890
seq_pt =	"CCCAAAAGGG"
struct_pt =	"(((....)))"

class FoldCompoundTest(unittest.TestCase):
	
	def test_create_fold_compound_Single(self):
		print "test_create_fold_compound_Single\n";
		
		fc = RNA.fold_compound(seq1);
		self.assertEqual(fc.type(),0);
		
	def test_create_fold_compound_Align(self):
		print "test_create_fold_compound_Align\n";
		fc= RNA.fold_compound(align);
		self.assertEqual(fc.type(),1);
		
	def test_create_fold_compound_2D(self):
		print "test_create_fold_compound_2D\n";
		fc= RNA.fold_compound(seq1,seq2,seq3);
		self.assertTrue(fc);
	
	def test_mfe(self):
		print "test_mfe\n";
		fc= RNA.fold_compound(seq1);
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ss,struct1);
	
	def test_mfe_Dimer(self):
		print "test_mfe_Dimer\n";
		fc=RNA.fold_compound(seq1Dimer);
		(ss,mfe) = fc.mfe_dimer();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ss,struct1Dimer);      
	def test_mfe_window(self):
		print "test_mfe_window\n";
		fc= RNA.fold_compound(seq1);
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ss,struct1);
	
	def test_eval_structure(self):	
		print "test_eval_structure\n";
		fc = RNA.fold_compound(seq1);
		energy= fc.eval_structure(struct1);
		self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60);
		print "\n", struct1, "[%6.2f" % energy,"]\n";
	
	# not working because unequal length of seq and structure OR converting from int to short is not valid 
	def test_eval_structure_pt(self):
		print "test_eval_structure_pt\n";
		fc=RNA.fold_compound(seq_pt);
		pt = [len(struct_pt),10,9,8,0,0,0,0,3,2,1];  #pairtable[0] = length of structure, 0 if no basepair, 
		energy= fc.eval_structure_pt(pt) /100; #/100 for dcal
		
		self.assertEqual("%6.2f" % energy, "%6.2f" % -2.5);
		print  struct_pt, "[%6.2f" % energy,"]\n";
	
	# Testing with filehandler and with stdout	
	def test_eval_structure_verbose(self):
		print "test_eval_structure_verbose";
		fc = RNA.fold_compound(seq1);
		try:
			f = open(filename, "w")
			print filename ," is opened for writing\n";
			energy = fc.eval_structure_verbose(struct1,f);
			energy2 = fc.eval_structure_verbose(struct1,None);
			
			self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60);
			print  struct1, "[%6.2f" % energy,"]\n";
		except IOError:
			print "Could not open ",filename;	
	
	
	#not working because of pairtbale and filehandler, and the energy is also not right(should be changed)
	def test_eval_structure_pt_verbose(self):
		print "test_eval_structure_pt_verbose\n";
		try:
			f = open(filename, "w")
			print filename ," is opened for writing\n";
			fc=RNA.fold_compound(seq_pt);
			pt = [len(struct_pt),10,9,8,0,0,0,0,3,2,1];  #pairtable[0] = length of structure, 0 if no basepair, 
			energy = fc.eval_structure_pt_verbose(pt,f)/100;  # / 100 for dcal
			energy2 = fc.eval_structure_pt_verbose(pt,None)/100;  # / 100 for dcal
			
			self.assertEqual("%6.2f" % energy, "%6.2f" % -2.5);
			print  struct_pt, "[%6.2f" % energy,"]\n";
		except IOError:
			print "Could not open ",filename;
	
	#!!!is not working, results segmentation fault
	def test_eval_covar_structure(self):
		print "test_eval_covar_structure\n";
		s1="CCCCAAAACGGG";
		s2="CCCGAAAAGGGG";
		s3="CCCCAAAAGGGG";
		ali = [s1,s2,s3];
		covarStructure = "((((....))))";
		

		fc = RNA.fold_compound(ali);
		pseudoEScore=fc.eval_covar_structure2(covarStructure);
		print covarStructure, "[ %6.2f" %pseudoEScore ,"]\n";
		self.assertTrue(pseudoEScore);
	
	#not workgin because of pairtable
	def test_eval_loop_pt(self):
		print "test_eval_loop_pt";
	
	def test_eval_move_del(self):
		print "test_eval_move_del";
		fc = RNA.fold_compound(seq1);
		energy = fc.eval_move(struct1,-7,-11);  # remove basepair (6,11) ,  energy change should be 2.10
		self.assertEqual("%6.2f" % energy, "%6.2f" % 2.10);
		print "\n", struct1, " moveset (-7,-11) --> [%6.2f" % energy ,"]\n";
		
	def test_eval_move_ins(self):
		print "test_eval_move_ins";
		fc = RNA.fold_compound(seq1);
		energy = fc.eval_move(struct11,7,11);  # add basepair (6,11) ,  energy change should be -2.10
		self.assertEqual("%6.2f" % energy, "%6.2f" % -2.10);
		print "\n", struct11, " moveset (7,11) --> [%6.2f" % energy,"]\n";	
	
	#not working, because of pairtable
	def test_eval_move_pt_del(self):
		print "test_eval_move_pt_del";
	
	#not workgin due to segmentation fault
	#def test_centroid(self):
		#print "test_centroid\n";
		#fc=RNA.fold_compound(align,None,RNA.VRNA_OPTION_MFE | RNA.VRNA_OPTION_EVAL_ONLY);
		#(cs,dist) = fc.centroid();
		#print cs, "\tDistance of :  %6.2f" %dist ,"\n";
	
	
	
	####contraints.h
	
	
	def test_constraints_add(self):
		seq_con  =  	"CCCAAAAGGGCCCAAAAGGG";
		str_con_def=	"(((....)))(((....)))";
		str_con=	"..........(((....)))";
		
		hc_file=	"hc.txt";
		print "test_constraints_add";
		fc = RNA.fold_compound(seq_con);
		fc.constraints_add(hc_file);
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ss,str_con);
		
		fc.hc_init();
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ss,str_con_def);
		
		
if __name__ == '__main__':
	unittest.main();
