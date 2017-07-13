import RNApath

RNApath.addSwigInterfacePath(3)


import RNA
import unittest

shortseq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUA"
longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"

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


def pf_window_callback(v, v_size, i, maxsize, what, data=None):
    if what & RNA.PROBS_WINDOW_UP:
        data['up'].append({ 'i': i, 'up': v})
    else:
        data['bpp'].extend([{'i': i, 'j': j, 'p': p} for j, p in enumerate(v) if (p is not None) and (p >= 0.01)])


def bpp_callback(v, v_size, i, maxsize, what, data=None):
    if what & RNA.PROBS_WINDOW_BPP:
        data.extend([{'i': i, 'j': j, 'p': p} for j, p in enumerate(v) if (p is not None) and (p >= 0.01)])


class pf_window_functionTest(unittest.TestCase):

    def test_pfl_fold(self):
        print("test_pfl_fold")
        bpp = RNA.pfl_fold(longseq, 200, 150, 0.01)

        # sanity check for base pair probabilities
        self.assertTrue(len(bpp) == 640);
        self.assertTrue(min([v.p for v in bpp]) >= 0.01)
        self.assertTrue(max([v.p for v in bpp]) <= 1.0)
        self.assertTrue(max([v.j for v in bpp]) == 300)


    def test_pfl_fold_up(self):
        print("test_pfl_fold_up")
        up = RNA.pfl_fold_up(longseq, 10, 200, 150)
        # sanity check for unpaired probabilities
        self.assertTrue(len(up) == 301)
        self.assertTrue(min([len(v) for v in up if v is not None]) == 11)
        self.assertTrue(max([len(v) for v in up if v is not None]) == 11)
        self.assertTrue(min([min([u for u in v]) for v in up]) >= 0.)
        self.assertTrue(max([max([u for u in v]) for v in up]) <= 1.)


    def test_probs_window(self):
        print("test_probs_window")
        data = { 'bpp': [], 'up': []}
        md = RNA.md()
        md.max_bp_span = 150
        md.window_size = 200
        fc = RNA.fold_compound(longseq, md, RNA.OPTION_WINDOW)
        fc.probs_window(10, RNA.PROBS_WINDOW_BPP | RNA.PROBS_WINDOW_UP, pf_window_callback, data)

        # sanity check for base pair probabilities
        self.assertTrue(len(data['bpp']) == 640)
        self.assertTrue(min([v['p'] for v in data['bpp']]) >= 0.01)
        self.assertTrue(max([v['p'] for v in data['bpp']]) <= 1.0)
        self.assertTrue(max([v['j'] for v in data['bpp']]) == 300)
        # sanity check for unpaired probabilities
        self.assertTrue(len(data['up']) == 300)
        self.assertTrue(min([len(v['up']) for v in data['up']]) == 11)
        self.assertTrue(max([len(v['up']) for v in data['up']]) == 11)
        self.assertTrue(min([min([vv for vv in v['up'] if vv is not None]) for v in data['up']]) >= 0.)
        self.assertTrue(max([max([vv for vv in v['up'] if vv is not None]) for v in data['up']]) <= 1.)


    def test_pfl_SHAPE(self):
        print("test_pfl_SHAPE")
        benchmark_set = ["Lysine_riboswitch_T._martima", "TPP_riboswitch_E.coli" ]

        for b in benchmark_set:
            seq  =  getShapeSequenceFromFile(datadir + b + ".db")
            reactivities = getShapeDataFromFile(datadir + b + ".shape_2rows")
            data = { 'bpp': [], 'up': []}

            md = RNA.md()
            md.max_bp_span = len(seq)
            md.window_size = len(seq)

            # compute pairing probabilities using local partition function implementation
            # with L = W = len(seq)
            fc = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)
            fc.sc_add_SHAPE_deigan(reactivities, 1.9, -0.7, RNA.OPTION_WINDOW)
            fc.probs_window(0, RNA.PROBS_WINDOW_BPP, pf_window_callback, data)

            # compute pairing probabilities using global parition function
            fc2 = RNA.fold_compound(seq)
            fc2.sc_add_SHAPE_deigan(reactivities, 1.9, -0.7, RNA.OPTION_DEFAULT)
            fc2.pf()
            bpp = fc2.bpp()

            # check for differences in pairing probabilities between local and global implementation
            # Hint: There must not be any!
            for prob in data['bpp']:
                p = prob['p']
                i = prob['i']
                j = prob['j']
                p2 = bpp[i][j]
                self.assertEqual("%2.8f" % abs(p - p2), "%2.8f" % 0.)


    """
    Compute partition function and base pair probabilities both, using the implementation
    for local structure prediction and global structure prediction.
    When comparing both results, equilibrium probabilities must not have changed!
    """
    def test_pfl_window_full(self):
        print("test_probs_window_full")

        data = []

        md = RNA.md()
        md.max_bp_span = len(shortseq)
        md.window_size = len(shortseq)

        # compute pairing probabilities using local partition function implementation
        # with L = W = len(seq)
        fc = RNA.fold_compound(shortseq, md, RNA.OPTION_WINDOW)
        fc.probs_window(0, RNA.PROBS_WINDOW_BPP, bpp_callback, data)

        # compute pairing probabilities using global parition function
        fc2 = RNA.fold_compound(shortseq)
        fc2.pf()
        bpp = fc2.bpp()

        # check for differences in pairing probabilities between local and global implementation
        # Hint: There must not be any!
        for prob in data:
            p = prob['p']
            i = prob['i']
            j = prob['j']
            p2 = bpp[i][j]
            self.assertEqual("%g" % p, "%g" % p2)


    """
    Compute partition function and base pair probabilities both, constrained
    and unconstrained, where the constraint simply shifts the free energy base
    line by -1 kcal/mol per nucleotide.
    When comparing both results, equilibrium probabilities must not have changed,
    except for free energy of the ensemble!
    """
    def test_pfl_sc(self):
        print("test_pfl_sc")

        fc = RNA.fold_compound(longseq, None, RNA.OPTION_WINDOW)

        data = []
        # unconstrained partition function
        fc.probs_window(0, RNA.PROBS_WINDOW_BPP, bpp_callback, data)

        # add constraints
        for i in range(1, len(longseq) + 1):
          fc.sc_add_up(i, -1.0, RNA.OPTION_WINDOW)

        for i in range(1, len(longseq)):
            for j in range(i + 1, len(longseq) + 1):
                fc.sc_add_bp(i, j, -2, RNA.OPTION_WINDOW)

        data2 = []
        # constrained partition function
        fc.probs_window(0, RNA.PROBS_WINDOW_BPP, bpp_callback, data2)

        # check pairing probabilities
        for du in data:
            dc = next(i for i in data2 if (i['i'] == du['i'] and i['j'] == du['j']))
            self.assertTrue(dc != None)
            self.assertEqual("%g" % du['p'], "%g" % dc['p'])



if __name__ == '__main__':
    unittest.main()
