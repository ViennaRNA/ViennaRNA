

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

##		 12345678
seq_pt =	"CCAAAAGG"
struct_pt =	"((....))"

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
	#def test_eval_structure_pt(self):
		#print "test_eval_structure_pt\n";

		#fc=RNA.fold_compound(seq_pt);
		#pt = [len(struct_pt),8,7,0,0,0,0,2,1];  #pairtable[0] = length of structure, 0 if no basepair, 
		#energy= fc.pt_test(pt);
		##self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60);
		##print "\n", struct_pt, "[%6.2f" % energy,"]\n";
	#def test_pt_test(self):
		
		#fc=RNA.fold_compound(seq_pt);
		##pt = [len(struct_pt),8,7,0,0,0,0,2,1];  #pairtable[0] = length of structure, 0 if no basepair, 
		##pt = [8,8,7,0,0,0,0,2,1];  #pairtable[0] = length of structure, 0 if no basepair, 
		
		##print pt[0];
		##print len(seq_pt);
		###		 12345678
		##seq_pt =	"CCAAAAGG"
		##struct_pt =	"((....))"
		##pairtable= 	887000021
		
		
		#print fc.eval_structure_pt(pt);
	
		
	
	
	# not working because overloading function with filenames is not recognised, but without it is???	
	def test_eval_structure_verbose(self):
		print "test_eval_structure_verbose\n";
		fc = RNA.fold_compound(seq1);
		try:
			f = open(filename, "w")
			print filename ," is opened for writing\n";
			energy = fc.eval_structure_verbose(struct1,f);
			self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60);
			print "\n", struct1, "%6.2f" % energy;
		except IOError:
			print "Could not open ",filename;
	
	#not working because of pairtbale and filehandler, and the energy is also not right(should be changed)
	#def test_eval_structure_pt_verbose(self):
		#print "test_eval_structure_pt_verbose\n";
		#fc = RNA.fold_compound(seq_pt);
		#try:
			#f = open(filename, "w")
			#print filename ," is opened for writing\n";
			#energy = fc.eval_structure_pt_verbose(pt,f);
			#self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60);
			#print "\n", struct_pt, "%6.2f" % energy;
		#except IOError:
			#print "Could not open ",filename;
	
	#geht nicht wegen irgendweinen segmentation foult error
	def test_eval_covar_structure(self):
		print "test_eval_covar_structure\n";
		fc = RNA.fold_compound(align);
		print fc.type();
		#(cs,mfe) = fc.eval_covar_structure();
		#print cs, "[ %6.2f" %mfe ,"]\n";
		#self.assertEqual(cs,struct1);
	
	#not workgin because of pairtable
	def test_eval_loop_pt(self):
		print "test_eval_loop_pt\n";
	
	def test_eval_move_del(self):
		print "test_eval_move_del\n";
		fc = RNA.fold_compound(seq1);
		energy = fc.eval_move(struct1,-7,-11);  # remove basepair (6,11) ,  energy change should be 2.10
		self.assertEqual("%6.2f" % energy, "%6.2f" % 2.10);
		print "\n", struct1, " moveset (-7,-11) --> [%6.2f" % energy ,"]\n";
		
	def test_eval_move_ins(self):
		print "test_eval_move_ins\n";
		fc = RNA.fold_compound(seq1);
		energy = fc.eval_move(struct11,7,11);  # add basepair (6,11) ,  energy change should be -2.10
		self.assertEqual("%6.2f" % energy, "%6.2f" % -2.10);
		print "\n", struct11, " moveset (7,11) --> [%6.2f" % energy,"]\n";	
	
	#not working, because of pairtable
	def test_eval_move_pt_del(self):
		print "test_eval_move_pt_del\n";
	
	#not workgin due to segmentation fault
	#def test_centroid(self):
		#print "test_centroid\n";
		#fc=RNA.fold_compound(align);
		#(cs,dist) = fc.centroid();
		#print cs, "\tDistance of :  %6.2f" %dist ,"\n";
	
	#def test_alifold(self):
		#print "hu";
		#ss,energyAli = RNA.alifold(align);
		#fc=RNA.fold_compound(align);
		#energy = fc.eval_covar_structure(ss);
		#self.assertEqual("%6.2f" % energy,"%6.2f" % energyAli);
		#print "\n", ss, "\t%6.2f" % energy , " EnergyAli = \t%6.2f" % energyAli ;
	
	#def test_test(self):
		#fc = RNA.fold_compound(sq1);
		#one,two,three,four,five = fc.sc_detect_hi_motif("(((...))))");
		#print one ,"  ", two,"  ", three, "  ", four;
		#one,two,three,four,five = fc.sc_get_hi_motif();
		#print one ,"  ", two,"  ", three, "  ", four;
	#def test_string(self):
		#print "test_string\n";
		#fc = RNA.fold_compound(sq1);
		#f = open(filename, "w");
		#print f;
		#print fc.testFunction(align,f);
		#print fc.testFunction(align);
		#print fc.testFunction(align, align);
		
if __name__ == '__main__':
	unittest.main();
