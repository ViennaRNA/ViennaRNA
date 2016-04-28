import RNApath

RNApath.addSwigInterfacePath(3)


import RNA # all global variables are in RNA.cvar
import unittest
from struct import *

seq1  =     "CGCAGGGAUACCCGCG"
struct1=    "(((.(((...))))))"
struct11 =  "(((.((.....)))))"
seq2  =     "GCGCCCAUAGGGACGC"
struct2=    "((((((...))).)))"
seq3  =     "GCGCACAUAGUGACGC"
struct3=    "(..(((...)))...)"
align=[seq1,seq2,seq3]
(struct, mfe) = RNA.fold(seq1)






class GeneralTests(unittest.TestCase):
    def test_hammingDistance(self):
        print("test_hammingDistance \t calculate a hamming distance")
        self.assertEqual(RNA.hamming(seq1,seq2),16)
        self.assertEqual(RNA.bp_distance(struct1,struct2),6)
    def test_temperature(self):
        print("test_temperature\n") 
        self.assertEqual(RNA.cvar.temperature,37) #!!!NOT WORKING !!! WHY
    def test_foldASequence(self):
        print("test_foldASequence\n")
        # new better interface
        (struct, mfe) = RNA.fold(seq1)
        self.assertEqual(struct,struct1)
        # check energy
        self.assertEqual(RNA.energy_of_struct(seq1,struct1), mfe)
    def test_constrained_folding(self):
        print("test_constrained_folding\n")
        RNA.cvar.fold_constrained = 1
        (struct,cmfe) = RNA.fold(seq1,"....xx....xx....")
        self.assertEqual(struct,'(((..........)))')
        self.assertEqual(RNA.energy_of_struct(seq1,struct),cmfe)
        RNA.cvar.fold_constrained = 0
    def test_cofold(self):
        print("test_cofold\n")
        RNA.cvar.cut_point = len(seq1)+1
        (costruct,comfe) = RNA.cofold(seq1 + seq2)
        self.assertEqual(costruct,'(((.(((...))))))((((((...))).)))')
        cmfe = RNA.energy_of_struct(seq1 + seq2, costruct)
        self.assertTrue(abs(comfe - cmfe) < 1e-5)
        (x,ac,bc,fcab,cf) = RNA.co_pf_fold(seq1 + seq2, struct)
        self.assertTrue(cf<comfe)
        self.assertTrue(comfe-cf<1.3)

        (x,usel1, usel2, fcaa, usel3)= RNA.co_pf_fold(seq1 + seq1, struct)
        RNA.cvar.cut_point = len(seq2)+1
        (x,usel1, usel2, fcbb, usel3)= RNA.co_pf_fold(seq2 + seq2, struct)
        (AB,AA,BB,A,B) = RNA.get_concentrations(fcab, fcaa, fcbb,ac, bc, 1e-5, 1e-5)
        AB/=2e-5
        AA/=2e-5
        BB/=2e-5
        A/=2e-5
        B/=2e-5
        
        ret = (abs(AB-0.0)+abs(AA-0.00578)+abs(BB-0.01100)+abs(A-0.48843)+abs(B-0.47801))
        self.assertTrue(ret<0.0001)
        RNA.cvar.cut_point=-1
        #pf_fo ld
        s,f = RNA.pf_fold(seq1,struct)
        self.assertTrue(f < mfe)
        self.assertTrue(mfe-f<0.8)
        
    def test_tree_distance(self):
        print("test_tree_distance\n");
        xstruc = RNA.expand_Full(struct1)
        T1 = RNA.make_tree(xstruc)
        xstruc = RNA.expand_Full(struct2)
        T2 = RNA.make_tree(xstruc)
        RNA.edit_backtrack = 1
        tree_dist = RNA.tree_edit_distance(T1, T2)
        # print RNA::get_aligned_line(0), RNA::get_aligned_line(1),"\n"
        self.assertEqual(tree_dist,2)
    
    #def test_check_access_C_array(self): # !!!NOT working, NONE in cvar.iindx
        #print("test_check_access_C_array\n");
        #print("WHAT ",RNA.cvar.iindx)
        #ret = RNA.intP_getitem(RNA.cvar.iindx,3)
        #print(ret);
        #self.assertEqual(ret,108);     
        #RNA.ushortP_setitem(RNA.cvar.xsubi, 0, 171);
        #RNA.ushortP_setitem(RNA.cvar.xsubi, 1, 42);
        #RNA.ushortP_setitem(RNA.cvar.xsubi, 2, 93);
        #self.assertEqual(RNA.cdata(RNA.cvar.xsubi, 6),pack('HHH',171,42,93));
        
     #def test_bp_prop(self): #geht auch nicht da iindex nicht geht 
         #print("test_bp_prop\n");
         

        #p1 = RNA.get_pr(2,15);
        #ii = RNA.intP_getitem($RNA::iindx, 2);
        #my $p2 = (RNA::pf_float_precision() != 0) ? RNA::floatP_getitem($RNA::pr, $ii-15) : RNA::doubleP_getitem($RNA::pr, $ii-15);
        #ok(($p1<0.999) && ($p1>0.99) && (abs($p1-$p2)<1.2e-7));

        #my $bpf = RNA::Make_bp_profile(length($seq1));
        #my @bpf = unpack("f*",RNA::cdata($bpf, length($seq1)*4*3));
        #ok (($bpf[2*3]+$bpf[2*3+1]>.99999)&&$bpf[2*3+2]==0 &&
            #($bpf[2*3+1]>=$p1));
        #my $pack = RNA::pack_structure($struc1);
        #is (RNA::unpack_structure($pack), $struc1);
    def test_parse_structure(self):
        print("test_parse_structure\n")
        RNA.parse_structure(struct1)
        self.assertEqual(RNA.cvar.loops,2)
        self.assertEqual(RNA.cvar.pairs,6)
        self.assertEqual(RNA.cvar.unpaired,4)
        self.assertEqual(RNA.intP_getitem(RNA.cvar.loop_degree,1),2)
    
    def test_rna_plots(self):
        print("test_rna_plots\n")
        RNA.PS_rna_plot(seq1, struct1, "test_ss.ps")
        anote = "2 15 1 gmark\n" + "3 cmark\n"
        RNA.PS_rna_plot_a(seq1, struct1, "test_ss_a.ps", None, anote)
        RNA.PS_dot_plot(seq1, "test_dp.ps")
        RNA.ssv_rna_plot(seq1, struct, "test.coord")
        print("please check the two postscript files test_ss.ps and test_dp.ps\n")
        RNA.write_parameter_file("test.par")
    
    def test_different_symbol_set(self):
        print("test_different_symbol_set\n")
        RNA.cvar.symbolset = "GC"
        start = RNA.random_string(len(struct1), "GC")
        (sinv, cost) = RNA.inverse_fold(start, struct1)
        (ss, en) = RNA.fold(sinv)
        self.assertEqual(ss, struct1)
    
    def test_suboptimal(self):
        print("test_suboptimal\n")
            
            
        RNA.free_pf_arrays()
        RNA.free_arrays()
        RNA.free_co_arrays()

        RNA.cvar.subopt_sorted = 1
        RNA.cvar.noLonelyPairs = 1
        solution = RNA.subopt(seq1, None, 500, None)

        print("%d suboptimals\n" % solution.size());
        for x in range(0,solution.size()):
        # the last access should produce a "value out of range" warning
            if(solution.get(x).structure) :
                print("%s,%6.2f\n" % (solution.get(x).structure,solution.get(x).energy))
            

        ## test native array output of subopt()
        solution = RNA.subopt(seq1, 500)
        print("%d suboptimals \n" % len(solution))
        for s in solution:
            print("%s %6.2f\n" % (s.structure,s.energy))
        

        solution = ""

        RNA.cvar.cut_point = 3
        e =  RNA.energy_of_struct("GCGC", "(())")
        self.assertEqual(int(e*100+0.5), 70)

        duplex = RNA.duplexfold(seq1, seq2)

        self.assertEqual(duplex.structure, ".(((.....(((.&)))))).")
        duplex=None
        
        align = ["GCCAUCCGAGGGAAAGGUU", "GAUCGACAGCGUCU-AUCG", "CCGUCUUUAUGAGUCCGGC"]
        (css, cen) = RNA.alifold(align)
        self.assertEqual(css,"(((.(((...)))..))).")
        #print(align[0]);
        #RNA.consens_mis(["GCCAUCCGAGGGAAAGGUU", "GAUCGACAGCGUCU-AUCG", "CCGUCUUUAUGAGUCCGGC"])  !!!TypeError: list must contain strings, WTF definetly string
        #self.assertEqual(RNA.consens_mis(align), "SMBHBHYDRBGDVWmVKBB")
        #RNA.free_alifold_arrays()
        
    #def test_moveSets(self): # !!!not working, because struct1_move is apperently not passed as reference
        #print("test_moveSets") 

        #RNA.cvar.cut_point=-1
        #struct1_move = "(..............)"
        ## move_standar( sequence, start structure, move_type(GRADIENT, FIRST, ADAPTIVE), verbosity, shifts, noLP)
        #RNA.move_standard(seq1, struct1_move, 0, 1, 0, 0)
        ##self.assertEqual(struct1_move, "................")


        #struct1_move = "(..............)"
        #RNA.move_standard(seq1, struct1_move, 1, 1, 0, 0)
        #self.assertEqual(struct1_move, "(((.((....)).)))")
    def test_simplexycoordinates(self):
        print("test_simplexycoordinates\n")
        
        coords = RNA.simple_xy_coordinates(struct1)

        for c in (coords):
            print(c.X, ",", c.Y, "\n")
        
    def test_model_details_structure(self):
        print("test_model_details_parameter_structure\n")
    
        # check model details structure
        md = RNA.md(); # default values
        self.assertEqual(int(md.dangles), 2)
        self.assertEqual(md.temperature, 37.0)

        RNA.cvar.dangles     = 0
        RNA.cvar.temperature = 40.1
        md = RNA.md("global") # global values
        self.assertEqual(int(md.dangles), 0)
        self.assertEqual(int(RNA.cvar.dangles), 0)
        self.assertEqual(md.temperature, 40.1)

        #reset globals to default
        RNA.cvar.dangles = 2
        RNA.cvar.temperature= 37.0
        
        # check parameter structures
        params = RNA.param()
        self.assertEqual(params.get_temperature(),37.0)
        params = RNA.param(md)
        self.assertEqual(params.get_temperature(),40.1)
        
        pf_params = RNA.exp_param()
        self.assertEqual(pf_params.get_temperature(),37.0)
        pf_params = RNA.exp_param(md)
        self.assertEqual(pf_params.get_temperature(),40.1)
        md = None
                
class FoldCompoundTest(unittest.TestCase):

    def test_create_fold_compound_Single(self):
        print("test_create_fold_compound_Single\n")

        fc = RNA.fold_compound(seq1)
        self.assertEqual(fc.type(),0)


    def test_create_fold_compound_Align(self):
        print("test_create_fold_compound_Align\n")
        fc= RNA.fold_compound(align)
        self.assertEqual(fc.type(),1)


    def test_create_fold_compound_2D(self):
        print("test_create_fold_compound_2D\n")
        fc= RNA.fold_compound(seq1,seq2,seq3)
        self.assertTrue(fc)


    ###centroid.h
    def test_centroid(self):
        print("test_centroid\n")
        fc=RNA.fold_compound(align)
        fc.pf()
        (sc,dist) = fc.centroid()
        print( sc,"\tDistance of :  %6.2f" %dist ,"\n")
        self.assertTrue(sc and dist)


    ## partition function from here
    def test_pf(self):
        print("test_pf")
        fc= RNA.fold_compound(seq1)
        (ss,gfe) = fc.pf()
        print(ss, "[ %6.2f" %gfe ,"]\n")
        self.assertTrue(ss)
        bp_dis = fc.mean_bp_distance()
        print(seq1 ,"\t meanBPDistance : ", bp_dis,"\n")
        self.assertTrue(bp_dis)


    # hairpin_loops.h from here
    #def test_eval_hp_loop(self):
        #print("test_eval_hp_loop")
        #seq1  =      "GCAAAAGG"
        #struct1=    ".(....)."

        #fc=RNA.fold_compound(seq1)
        ##ehair = fc.eval_hp_loop(2,7)
        #ehair = fc.E_hp_loop(2,7)
        #print(seq1, " 2,7  = [ %6.2f" %ehair ,"] \n")
        #self.assertEqual("%6.2f" %ehair,"%6.2f" % +410)

        #Exterior loop evalution not working
        #eExt = fc.E_ext_hp_loop(0,0)
        #print(seq1, " 2,6  = [ %6.2f" %eExt ,"] \n")
        #self.assertEqual("%6.2f" %eExt,"%6.2f" % -140)

        #length = 8
        #External loop                           :  -140
        #Hairpin  loop (  2,  7) CG              :   410
        #GCAAAAGG
        #.(....).
        #energy =   2.70


    #def test_exp_E_hp_loop(self):
        #print("test_exp_E_hp_loop")
        #seq1  =      "GCAAAAGG"
        #struct1=    ".(....)."

        #fc=RNA.fold_compound(seq1)
        #fc.pf()
        #ehair = fc.exp_E_hp_loop(2,7)
        #print(seq1, " 2,7  = [ %6.2f" %ehair ,"] \n")
        #self.assertEqual("%6.2f" %ehair,"%6.2f" % +410)

    #interior_loops.h from here
    #def test_E_int_loop(self):
        #print("test_E_int_loop")
        ##    "123456789012"
        #seq1 =  "AGACAAAAGACA"
        #struct1=".(.(....).)."
        #fc=RNA.fold_compound(seq1, None, RNA.OPTION_MFE)
        #e = fc.E_int_loop(2,11)
        #print(seq1, " 2,7  = [ %6.2f" %e ,"] \n")
        #self.assertEqual("%6.2f" %e,"%6.2f" % +80)


if __name__ == '__main__':
    unittest.main()
