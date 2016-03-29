

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


def getShapeDataFromFile(filepath):
	
	retVec = [];
	with open(filepath, 'r') as f:
		lines = f.readlines();
		
		for line in lines:
			retVec.append(float(line.split(' ')[1]));
	return retVec;

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
	
	#adding the constraint is not correct, function is wrong????, with enforce options
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
		
		
	#
	#def test_sc_add_SHAPE_deigan(self):
		#print "test_sc_add_SHAPE_deigan";
		#seq_telomerase  =  	"GGGCUGUUUUUCUCGCUGACUUUCAGCCCAACACAAAAAAAGUCAGC";  # directly /home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/sequences_ct_files/Telomerase_pseudoknot_human.fa
		#fc=RNA.fold_compound(seq_telomerase);
		#reactivities = getShapeDataFromFile("/home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/shape_data/Telomerase_pseudoknot_human.shape_2rows");
		#ret = fc.sc_add_SHAPE_deigan(reactivities,1.8,-0.6,RNA.VRNA_OPTION_MFE);
		#(ss,mfe) = fc.mfe();
		#print ss, "[ %6.2f" %mfe ,"]\n";
		#self.assertEqual(ret,1);
		
	def test_sc_add_SHAPE_deigan(self):
		print "test_sc_add_SHAPE_deigan";
		seq_ribo  =  	"GACUCGGGGUGCCCUUCUGCGUGAAGGCUGAGAAAUACCCGUAUCACCUGAUCUGGAUAAUGCCAGCGUAGGGAAGUUC";  # directly /home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/sequences_ct_files/TPP_riboswitch_E.coli.db
		fc=RNA.fold_compound(seq_ribo);
		reactivities = getShapeDataFromFile("/home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/shape_data/TPP_riboswitch_E.coli.shape_2rows");
		ret = fc.sc_add_SHAPE_deigan(reactivities,1.9,-0.7,RNA.VRNA_OPTION_MFE);
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ret,1);

	
	#only for testing the function, none real functionality
	
	
	# we need real shape files with alignments, now we have unequal sequence length
	#def test_sc_add_SHAPE_deigan_ali(self):
		#print "test_sc_add_SHAPE_deigan_ali";
		##Sequence Ali from /home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/sequences_ct_files/5domain16S_rRNA_E.coli.fa
		## and    	   /home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/sequences_ct_files/5domain16S_rRNA_H.volcanii.fa
		## shape data from 
		#shapeData1 = "/home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/5domain16S_rRNA_E.coli.shape_2rows";
		#shapeData2 = "/home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/5domain16S_rRNA_H.volcanii.shape_2rows";
		#shapeAli = ["GAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGUCUGGGAAACUGCCUGAUGGAGGGGGAUAACUACUGGAAACGGUAGCUAAUACCGCAUAACGUCGCAAGACCAAAGAGGGGGACCUUCGGGCCUCUUGCCAUCGGAUGUGCCCAGAUGGGAUUAGCUAGUAGGUGGGGUAACGGCUCACCUAGGCGACGAUCCCUAGCUGGUCUGAGAGGAUGACCAGCCACACUGGAACUGAGACACGGUCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCACAAUGGGCGCAAGCCUGAUGCAGCCAUGCCGCGUGUAUGAAGAAGGCCUUCGGGUUGUAAAGUACUUUCAGCGGGGAGGAAGGGAGUAAAGUUAAUACCUUUGCUCAUUGACGUUACCCGCAGAAGAAGCACCGGCUAACUCCGUGCCAGCAGCCGCGGUAAUACGGAGGGUGCAAGCGUUAAUC","GGUCAUUGCUAUUGGGGUCCGAUUUAGCCAUGCUAGUUGCACGAGUUCAUACUCGUGGCGAAAAGCUCAGUAACACGUGGCCAAACUACCCUACAGAGAACGAUAACCUCGGGAAACUGAGGCUAAUAGUUCAUACGGGAGUCAUGCUGGAAUGCCGACUCCCCGAAACGCUCAGGCGCUGUAGGAUGUGGCUGCGGCCGAUUAGGUAGACGGUGGGGUAACGGCCCACCGUGCCGAUAAUCGGUACGGGUUGUGAGAGCAAGAGCCCGGAGACGGAAUCUGAGACAAGAUUCCGGGCCCUACGGGGCGCAGCAGGCGCGAAACCUUUACACUGCACGCAAGUGCGAUAAGGGGACCCCAAGUGCGAGGGCAUAUAGUCCUCGCUUUUCUCGACCGUAAGGCGGUCGAGGAAUAAGAGCUGGGCAAGACCGGUGCCAGCCGCCGCGGUAAUACCGGCAGCUCAAGUGAUGACC"];
		#fc=RNA.fold_compound(shapeAli);   
		#assoc = [1,2];
		#ret = fc.sc_add_SHAPE_deigan_ali([shapeData1,shapeData2], assoc,1.8,-0.6);
		#(ss,mfe) = fc.mfe();
		#print ss, "[ %6.2f" %mfe ,"]\n";
		#self.assertEqual(ret,1);
					
	
	def test_sc_add_SHAPE_zarringhalam(self):
		print "test_sc_add_SHAPE_zarringhalam";
		seq_con  =  	"GACUCGGGGUGCCCUUCUGCGUGAAGGCUGAGAAAUACCCGUAUCACCUGAUCUGGAUAAUGCCAGCGUAGGGAAGUUC";# directly /home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/sequences_ct_files/TPP_riboswitch_E.coli.db
		fc=RNA.fold_compound(seq_con);
		reactivities = getShapeDataFromFile("/home/mescalin/mario/projects/interfaces/ShapeKnots_DATA/shape_data/TPP_riboswitch_E.coli.shape_2rows"); 
		ret = fc.sc_add_SHAPE_zarringhalam(reactivities,0.5,0.5,"M"); 
		(ss,mfe) = fc.mfe();
		print ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ret,1);
	
	# start with ligand.h
	# test theophylline ligand binding interface
	
	def test_sc_add_hi_motif(self):
		print "test_sc_add_hi_motif";
		fc= RNA.fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG");
		#struct =	      ".(((((..((((((((((((((((((....(((((((............)))))))........)))))))))))))...))))))))))..............";
		#structWithMotif=      "((((...((((((((......)))))...)))...))))...((.((((((((((.((((((.((.((((....)))).))..)))))).))))))))))))..";
		ret = fc.sc_add_hi_motif("GAUACCAG&CCCUUGGCAGC","(...((((&)...)))...)",-9.22);
		(ss,mfe) = fc.mfe();
		print ret,"\t",ss, "[ %6.2f" %mfe ,"]\n";
		self.assertEqual(ret,1);
		
	
	# test not possible because data from soft constraints from folding comopund has to be set, with set data
	#def test_sc_detect_hi_motif(self):
		#print "test_sc_detect_hi_motif";
		#fc= RNA.fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG");
		##(ret,i,j,k,l)= fc.sc_detect_hi_motif("(...((((&)...)))...)");
		#ret= fc.sc_detect_hi_motif("(...((((&)...)))...)");
		##print "\n",ret,"\t",i,"\t",j,"\t",k,"\t",l;
		#self.assertEqual(ret,1);
	
	# test not possible because data from soft constraints from folding comopund has to be set, with set data
	#def test_sc_get_hi_motif(self):
		#print "test_sc_get_hi_motif";
		#fc= RNA.fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG");
		#(ret,i,j,k,l)= fc.sc_get_hi_motif();
		#print "\n",ret,"\t",i,"\t",j,"\t",k,"\t",l;
		#self.assertEqual(ret,1);
		
if __name__ == '__main__':
	unittest.main();



