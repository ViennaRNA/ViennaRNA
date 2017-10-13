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
        lines = f.read().splitlines()

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
        lines = f.read().splitlines()

    return lines[0]

class constraintsTest(unittest.TestCase):

    def test_sc_add_deigan(self):
        print "test_sc_add_deigan"
        seq          = getShapeSequenceFromFile(datadir + "Lysine_riboswitch_T._martima.db")
        reactivities = getShapeDataFromFile(datadir + "Lysine_riboswitch_T._martima.shape_2rows")

        fc=RNA.fold_compound(seq)
        print reactivities

        ret = fc.sc_add_SHAPE_deigan(reactivities, 1.9, -0.7, RNA.OPTION_DEFAULT)
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.4f ]" % mfe
        self.assertEqual("%6.4f" %mfe,"%6.4f" % -121.55 )


    def test_sc_add_SHAPE_deigan2(self):
        print "test_sc_add_SHAPE_deigan2"
        seq_ribo     = getShapeSequenceFromFile(datadir + "TPP_riboswitch_E.coli.db")
        reactivities = getShapeDataFromFile(datadir + "TPP_riboswitch_E.coli.shape_2rows")

        fc=RNA.fold_compound(seq_ribo)
        print reactivities

        ret = fc.sc_add_SHAPE_deigan(reactivities,1.9,-0.7,RNA.OPTION_DEFAULT);  # these values were copied from ronnys Talk about constraints
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f ]" % mfe
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
        print ss, "[ %6.2f ]" % mfe
        self.assertEqual(ret,1)


    def test_sc_add_SHAPE_zarringhalam(self):
        print "test_sc_add_SHAPE_zarringhalam"
        seq_ribo     = getShapeSequenceFromFile(datadir + "TPP_riboswitch_E.coli.db")
        fc           = RNA.fold_compound(seq_ribo)
        reactivities = getShapeDataFromFile(datadir + "TPP_riboswitch_E.coli.shape_2rows")

        ret = fc.sc_add_SHAPE_zarringhalam(reactivities, 0.5, 0.5, "O"); # these values were copied from ronnys Talk about constraints, O = default value
        (ss,mfe) = fc.mfe()
        print ss, "[ %6.2f ]" % mfe
        self.assertEqual("%6.2f" %mfe,"%6.2f" % -5.28 )



if __name__ == '__main__':
    unittest.main()
