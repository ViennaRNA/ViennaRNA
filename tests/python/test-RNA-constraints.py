import RNApath

RNApath.addSwigInterfacePath()

import RNA
import unittest

seq_con     = "CCCAAAAGGGCCCAAAAGGG"
str_con     = "..........(((....)))"
str_con_def = "(((....)))(((....)))"

datadir = RNApath.getDataDirPath()

def getShapeDataFromFile(filepath):
    retVec = []
    retVec.append(-999.0);  # data list is 1-based, so we add smth. at pos 0
    count=1
    with open(filepath, 'r') as f:
        lines = f.readlines()

        for line in lines:
            pos = int(line.split(' ')[0])
            value = float(line.split(' ')[1])

            if(pos==count):
                retVec.append(value)
            else:
                for i in range(pos-count):
                    retVec.append(-999.0)
                retVec.append(value)
                count=pos
            count+=1
    return retVec


def getShapeSequenceFromFile(filepath):
    retSeq=""
    with open(filepath,'r') as f:
        lines = f.readlines()

    return lines[0]

class constraintsTest(unittest.TestCase):

    def test_constraints_add(self):
        #seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        #str_con_def=    "(((....)))(((....)))"
        #hc.txt=    "P 1 0 2"
        #str_con=    "..........(((....)))"

        hc_file = datadir + "hc.txt"
        print "test_constraints_add"
        fc = RNA.fold_compound(seq_con)
        fc.constraints_add(hc_file)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual(ss,str_con)

        fc.hc_init()
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual(ss,str_con_def)

        #sc.txt = E 3 8 1 -5
        sc_file = datadir + "sc.txt"
        fc.sc_init()
        fc.constraints_add(sc_file)
        (ss,mfeNew) = fc.mfe()
        print ss, "[ %6.2f" %mfeNew ,"]\n"
        self.assertEqual("%6.2f" %mfe, "%6.2f" % (mfeNew +5))


    def test_hc_add_up(self):
        #seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        #str_con_def=    "(((....)))(((....)))"
        #str_con=    "..........(((....)))"
        print "test_hc_add_up\n"
        fc = RNA.fold_compound(seq_con)
        fc.hc_add_up(1,RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual(ss,".((....)).(((....)))")


    def test_hc_add_bp_nonspecific(self):
        print "test_hc_add_bp_nonspecific"
        #GGGCCCCCCCCCCCCCCCCC
        #(((......)))........
        fc=RNA.fold_compound("GGGCCCCCCCCCCCCCCCCC")
        fc.hc_add_bp_nonspecific(20,-1); # force the last base to pair with some bases upstream
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual(ss,"(((..............)))")


    def test_hc_add_bp(self):
        print "test_hc_add_bp"
        seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        str_con_def=    "(((....)))(((....)))"
        fc=RNA.fold_compound(seq_con)
        fc.hc_add_bp(1,20,RNA.CONSTRAINT_CONTEXT_ENFORCE | RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual(ss,"(((..............)))")


    def test_hc_add_from_db(self):
        print "test_hc_add_from_db"
        #seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
        #str_con_def=    "(((....)))(((....)))"
        #hc.txt=    "xxx................."
        #str_con=    "..........(((....)))"
        fc = RNA.fold_compound(seq_con)
        fc.hc_add_from_db("xxx.................")
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual(ss,str_con)


    def test_sc_set_up(self):
        print "test_sc_set_up"

        #        "1234567890
        seq_sc  =      "CCCAAAAGGG"
        fc = RNA.fold_compound(seq_sc)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual(ss,"(((....)))")

        #fc.sc_init()

        m= [0.0,-5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];  #E 1 0 1 -1 ,  position 1 gets -5 if unpaired ,".((....))." structure should be prefered!!Attention vector starts with position 0

        fc.sc_set_up(m)
        (ss,mfeNew) = fc.mfe()
        print ss, "[ %6.2f" %mfeNew ,"]\n"
        self.assertEqual("%6.2f" %mfeNew,"%6.2f" % -5.70)


    def test_sc_set_bp(self):
        print "test_sc_set_bp"

        #add energy of -5 to basepair 1-9 if formed, prefed structure should now be ((.....))., with a energy of -4.90
        m = [[0 for x in range(11)] for y in range(11)]
        m[1][9] = -5.0; # base 1-9 should get -5.0 if basepair
        m[9][1] = -5.0; # base 1-9 should get -5.0 if basepair


        seq_sc = "CCCAAAAGGG"
        fc = RNA.fold_compound(seq_sc)
        fc.sc_set_bp(m)
        (ss,mfeNew) = fc.mfe()
        print ss, "[ %6.2f" %mfeNew ,"]\n"
        self.assertEqual("%6.2f" %mfeNew,"%6.2f" % -4.90)


    def test_sc_add_deigan(self):
        print "test_sc_add_deigan"
        seq  =  getShapeSequenceFromFile(datadir + "Lysine_riboswitch_T._martima.db")
        reactivities = getShapeDataFromFile(datadir + "Lysine_riboswitch_T._martima.shape_2rows")

        fc=RNA.fold_compound(seq)
        print reactivities

        ret = fc.sc_add_SHAPE_deigan(reactivities,1.9,-0.7,RNA.OPTION_DEFAULT)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.4f" %mfe ,"]\n"
        self.assertEqual("%6.4f" %mfe,"%6.4f" % -121.55 )


    def test_sc_add_SHAPE_deigan2(self):
        print "test_sc_add_SHAPE_deigan2"
        seq_ribo  = getShapeSequenceFromFile(datadir + "TPP_riboswitch_E.coli.db")
        reactivities = getShapeDataFromFile(datadir + "TPP_riboswitch_E.coli.shape_2rows")

        fc=RNA.fold_compound(seq_ribo)
        print reactivities

        ret = fc.sc_add_SHAPE_deigan(reactivities,1.9,-0.7,RNA.OPTION_DEFAULT);  # these values were copied from ronnys Talk about constraints
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual("%6.2f" %mfe,"%6.2f" % -52.61 )


    # I just added completely randomly choosen sequences and shape data, this test only checks if the function can be called, not if it results correct results.
    def test_sc_add_SHAPE_deigan_ali(self):
        print "test_sc_add_SHAPE_deigan_ali"
        # shape data from
        shapeSeq1 = "CCCAAAAGGG"
        shapeSeq2 = "AAUAAAAAUU"

        shapeData1 = [-999.0,0.04,1.12,1.1,-999.0,0.05,0.5,0.3,-999.0,1.4]
        shapeData2 = [-999.0,-999.0,-999.0,1.23,1.4,0.05,0.5,0.3,-999.0,1.4]
        shapeAli = [shapeSeq1,shapeSeq2]
        fc=RNA.fold_compound(shapeAli)
        assoc = [-1,1,2]
        ret = fc.sc_add_SHAPE_deigan_ali(shapeAli, assoc,1.8,-0.6)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual(ret,1)


    def test_sc_add_SHAPE_zarringhalam(self):
        print "test_sc_add_SHAPE_zarringhalam"
        seq_ribo  = getShapeSequenceFromFile(datadir + "TPP_riboswitch_E.coli.db")
        fc=RNA.fold_compound(seq_ribo)
        reactivities = getShapeDataFromFile(datadir + "TPP_riboswitch_E.coli.shape_2rows")
        print reactivities
        ret = fc.sc_add_SHAPE_zarringhalam(reactivities,0.5,0.5,"O"); # these values were copied from ronnys Talk about constraints, O = default value
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual("%6.2f" %mfe,"%6.2f" % -5.28 )


    def test_sc_add_hi_motif(self):
        print "test_sc_add_hi_motif"
        fc= RNA.fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG")
        #struct =          ".(((((..((((((((((((((((((....(((((((............)))))))........)))))))))))))...)))))))))).............."
        #structWithMotif=      "((((...((((((((......)))))...)))...))))...((.((((((((((.((((((.((.((((....)))).))..)))))).)))))))))))).."
        ret = fc.sc_add_hi_motif("GAUACCAG&CCCUUGGCAGC","(...((((&)...)))...)",-9.22)
        (ss,mfe) = fc.mfe()
        print ret,"\t",ss, "[ %6.2f" %mfe ,"]\n"
        self.assertEqual(ret,1)


    def test_theophylline_ligand_binding_interface(self):
        print("test_theophylline_ligand_binding_interface\n")
        RNA.noLonelyPairs = 0
        fc = RNA.fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG")
        (ss, mfe) = fc.mfe()
        print("%s [ %6.2f ]\n" % (ss, mfe))

        fc.sc_add_hi_motif("GAUACCAG&CCCUUGGCAGC", "(...((((&)...)))...)", -9.22)
        (ss, mfe) = fc.mfe()
        print("%s [ %6.2f ]\n" % (ss, mfe))

        fc.sc_remove()

        fc.sc_add_hi_motif("GAAAAAU", "(.....)", -19)
        (ss, mfe) = fc.mfe()
        print("%s [ %6.2f ]\n" %(ss, mfe))


if __name__ == '__main__':
    unittest.main()
