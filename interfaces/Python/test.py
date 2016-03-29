

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
struct1_pt =	 [len(struct1),16,15,14,0,13,12,11,0,0,0,7,6,5,3,2,1];

seq_con  =  	"CCCAAAAGGGCCCAAAAGGG";
str_con_def=	"(((....)))(((....)))";
str_con=	"..........(((....)))";


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
	
	## 
	def test_eval_structure_pt(self):
		print "test_eval_structure_pt\n";
		fc=RNA.fold_compound(seq1);
		energy= fc.eval_structure_pt(struct1_pt) /100; #/100 for dcal
		
		self.assertEqual("%6.2f" % energy, "%6.2f" % -5.60);
		print  struct1, "[%6.2f" % energy,"]\n";
	
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
	
	
	
	#def test_eval_structure_pt_verbose(self):
		print "test_eval_structure_pt_verbose\n";
		try:
			f = open(filename, "w")
			print filename ," is opened for writing\n";
			fc=RNA.fold_compound(seq1);
			
			energy = fc.eval_structure_pt_verbose(struct1_pt,f)/100;  # / 100 for dcal
			energy2 = fc.eval_structure_pt_verbose(struct1_pt,None)/100;  # / 100 for dcal
			
			self.assertEqual("%6.2f" % energy, "%6.2f" % -5.6);
			print  struct1, "[%6.2f" % energy,"]\n";
		except IOError:
			print "Could not open ",filename;
	
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
	
	def test_eval_loop_pt(self):
		print "test_eval_loop_pt";		
		fc= RNA.fold_compound(seq1);
		energy= fc.eval_loop_pt(6,struct1_pt) /100; #/100 for dcal
		print "[ %6.2f" %energy ,"]\n";
		self.assertEqual("%6.2f" % energy,"%6.2f" % -3.3);
	
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
	
	def test_eval_move_pt_del(self):
		print "test_eval_move_pt_del";
		fc = RNA.fold_compound(seq1);
		energy = fc.eval_move_pt(struct1_pt,-7,-11) /100;  # remove basepair (6,11) ,  energy change should be 2.10
		self.assertEqual("%6.2f" % energy, "%6.2f" % 2.10);
		print "\n", struct1, " moveset (-7,-11) --> [%6.2f" % energy ,"]\n";
		
	

	#not workgin due to segmentation fault
	#def test_centroid(self):
		#print "test_centroid\n";
		#fc=RNA.fold_compound(align);
		#(cs,dist) = fc.centroid();
		#print dist;
		#print cs, "\tDistance of :  %6.2f" %dist ,"\n";
	
	
	
	####constraints.h
	
	
	def test_constraints_add(self):
		#seq_con  =  	"CCCAAAAGGGCCCAAAAGGG";
		#str_con_def=	"(((....)))(((....)))";
		#hc.txt=	"xxx.................";
		#str_con=	"..........(((....)))";
		
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
	
	###!!!!try toget the COntraint option to python
	def test_hc_add_up(self):
		#seq_con  =  	"CCCAAAAGGGCCCAAAAGGG";
		#str_con_def=	"(((....)))(((....)))";
		#str_con=	"..........(((....)))";
		print "test_hc_add_up\n";
		fc = RNA.fold_compound(seq_con);
		#fc.hc_add_up(1,RNA.VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);
		fc.hc_add_up(1);
		
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ss,".((....)).(((....)))");
		
	def test_hc_add_bp_nonspecific(self):
		print "test_hc_add_bp_nonspecific";
		#GGGCCCCCCCCCCCCCCCCC
		#(((......)))........
		fc=RNA.fold_compound("GGGCCCCCCCCCCCCCCCCC");
		fc.hc_add_bp_nonspecific(20,-1); # force the last base to pair with some bases upstream
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ss,"(((..............)))");
	
	#adding the constraint is not correct, function is wrong????
	def test_hc_add_bp(self):
		print "test_hc_add_bp";
		seq_con  =  	"CCCAAAAGGGCCCAAAAGGG";
		str_con_def=	"(((....)))(((....)))";
		fc=RNA.fold_compound(seq_con);
		fc.hc_add_bp(3,20);
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ss,"...........((....)).");
		
	
	#try to get the VRNA_CONSTRAINT_DB_DEFAULT options to python
	def test_hc_add_from_db(self):
		print "test_hc_add_from_db";
		#seq_con  =  	"CCCAAAAGGGCCCAAAAGGG";
		#str_con_def=	"(((....)))(((....)))";
		#hc.txt=	"xxx.................";
		#str_con=	"..........(((....)))";
		fc = RNA.fold_compound(seq_con);
		fc.hc_add_from_db("xxx.................");
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ss,str_con);
		
		
	#change to realistic data
	def test_sc_add_SHAPE_deigan(self):
		print "test_sc_add_SHAPE_deigan";
		seq_con  =  	"CCCAAAAGGGCCCAAAAGGG";
		fc=RNA.fold_compound(seq_con);
		reactivities =[4.4,3.3,3.3];  # just random stuff, should be changed to realistic data
		ret = fc.sc_add_SHAPE_deigan(reactivities,1.8,-0.6,RNA.VRNA_OPTION_MFE);
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ret,1);
	
	#only for testing the function, none real functionality
	
	
	# we need real shape files
	#def test_sc_add_SHAPE_deigan_ali(self):
		#print "test_sc_add_SHAPE_deigan_ali";
		#fc=RNA.fold_compound(align);   
		#files = ["shape_deigan_file1","shape_deigan_file2","shape_deigan_file3"];
		#assoc = [1,2,3];
		#ret = fc.sc_add_SHAPE_deigan_ali(files, assoc,1.8,-0.6);
		#(ss,mfe) = fc.mfe();
		#print ss, "[ %6.2f" %mfe ,"]\n";
		#self.assertEqual(ret,1);
					
	
	def test_sc_add_SHAPE_zarringhalam(self):
		print "test_sc_add_SHAPE_zarringhalam";
		seq_con  =  	"CCCAAAAGGGCCCAAAAGGG";
		fc=RNA.fold_compound(seq_con);
		reactivities =[4.4,3.3,3.3];  
		ret = fc.sc_add_SHAPE_zarringhalam(reactivities,0.6,0.5,"M"); # just random stuff, should be changed to realistic data
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ret,1);
	
		
if __name__ == '__main__':
	unittest.main();
