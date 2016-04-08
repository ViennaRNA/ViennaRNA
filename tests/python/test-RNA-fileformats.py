import RNApath

RNApath.addSwigInterfacePath()

import RNA
import unittest

def getShapeDataFromFile(filepath):
    retVec = []
    retVec.append(-999.0); # data list is 1-based, so we add smth. at pos 0
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
    for i,v in enumerate(retVec):
        print "%d\t%6.6f" % (i, v)

    return retVec


class constraintsTest(unittest.TestCase):

    def test_file_SHAPE_read(self):
        print "test_file_SHAPE_read"
        reactivities = getShapeDataFromFile("data/TPP_riboswitch_E.coli.shape_2rows")

        (a,b,c) = RNA.file_SHAPE_read("data/TPP_riboswitch_E.coli.shape_2rows", 79, -1)
        print "read file:"
        print a
        print b
        print c

        print reactivities
        print a


if __name__ == '__main__':
    unittest.main()
