import os
from operator import add
import math
import unittest
import sys

if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA


kT = 0.61632077549999997
# maximum allowed difference beteen compared probabilties
allowed_diff = 1e-7

# some sequences to work with
shortseq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUA"
longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"

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


def up_callback(v, v_size, i, maxsize, what, data):
    if what & RNA.PROBS_WINDOW_UP:
        data.append({ 'i': i, 'up': v})


def up_split_callback(v, v_size, i, maxsize, what, data):
    if what & RNA.PROBS_WINDOW_UP:
        what = what & ~RNA.PROBS_WINDOW_UP
        dat = []
        # Non-split case:
        if what == RNA.ANY_LOOP:
                dat = data
        # all the cases where probability is split into different loop contexts
        elif what == RNA.EXT_LOOP:
                dat = data['ext']
        elif what == RNA.HP_LOOP:
                dat = data['hp']
        elif what == RNA.INT_LOOP:
                dat = data['int']
        elif what == RNA.MB_LOOP:
                dat = data['mb']
        dat.append({'i': i, 'up': v})


class pf_window_functionTest(unittest.TestCase):
    DATADIR = os.environ.get('VRNA_TEST_DATA', os.sep.join(["tests", "data"]))

    def test_pfl_fold(self):
        """RNA.pfl_fold() - sanity check for base pair probabilities"""
        bpp = RNA.pfl_fold(longseq, 200, 150, 0.01)

        # sanity check for base pair probabilities
        self.assertTrue(len(bpp) == 640);
        self.assertTrue(min([v.p for v in bpp]) >= 0.01)
        self.assertTrue(max([v.p for v in bpp]) <= 1.0)
        self.assertTrue(max([v.j for v in bpp]) == 300)


    def test_pfl_fold_up(self):
        """RNA.pfl_fold_up() - sanity check for unpaired probabilities"""
        up = RNA.pfl_fold_up(longseq, 10, 200, 150)
        # sanity check for unpaired probabilities
        self.assertTrue(len(up) == 301)
        self.assertTrue(min([len(v) for v in up if v is not None]) == 11)
        self.assertTrue(max([len(v) for v in up if v is not None]) == 11)
        self.assertTrue(min([min([u for u in v]) for v in up]) >= 0.)
        self.assertTrue(max([max([u for u in v]) for v in up]) <= 1.)


    def test_probs_window(self):
        """fold_compound.probs_window() - sanity checks for probabiities"""
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
        """SHAPE data"""
        benchmark_set = ["Lysine_riboswitch_T._martima", "TPP_riboswitch_E.coli" ]

        for b in benchmark_set:
            seq  =  getShapeSequenceFromFile(os.sep.join([self.DATADIR, b + ".db"]))
            reactivities = getShapeDataFromFile(os.sep.join([self.DATADIR, b + ".shape_2rows"]))
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


    def test_pfl_window_full(self):
        """Compare results from local and global predictions
        
        Compute partition function and base pair probabilities both, using the implementation
        for local structure prediction and global structure prediction.
        When comparing both results, equilibrium probabilities must not have changed!
        """

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


    def test_pfl_sc(self):
        """Soft constraints - energy landscape shift - base pair probabilities
        
        Compute partition function and base pair probabilities both, constrained
        and unconstrained, where the constraint simply shifts the free energy base
        line by -1 kcal/mol per nucleotide.
        When comparing both results, equilibrium probabilities must not have changed,
        except for free energy of the ensemble!
        """

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


    @unittest.skipIf(sys.platform.startswith("win"), "May fail under Windows due to numeric imprecision")
    def test_probs_window(self):
        """Soft constraints - energy landscape shift - unpaired probabilities
        
        Compute unpaired probabilities both, constrained and unconstrained, where the
        constraint simply shifts the free energy base line by -1 kcal/mol per nucleotide.
        When comparing both results, equilibrium probabilities must not have changed!
        """

        ulength = 45
        data = []
        data2 = []

        md = RNA.md()
        md.max_bp_span = 150
        md.window_size = 200

        fc = RNA.fold_compound(longseq, md, RNA.OPTION_WINDOW)

        # unconstrained unpaired probabilities
        fc.probs_window(ulength, RNA.PROBS_WINDOW_UP, up_callback, data)

        # add constraints
        for i in range(1, len(longseq) + 1):
          fc.sc_add_up(i, -1.0, RNA.OPTION_WINDOW)

        for i in range(1, len(longseq)):
            for j in range(i + 1, len(longseq) + 1):
                fc.sc_add_bp(i, j, -2, RNA.OPTION_WINDOW)

        # constrained unpaired probabilities
        fc.probs_window(ulength, RNA.PROBS_WINDOW_UP, up_callback, data2)

        for i in range(1, len(longseq) + 1):
            for u in range (1, ulength + 1):
                if i - u + 1 <= 0:
                    break
                self.assertEqual("%g" % data[i - 1]['up'][u], "%g" % data2[i - 1]['up'][u])

    def test_up_split_sum(self):
        """Reporting check for probs_window() callback

        Compute unpaired probabilities both, split into different loop contexts and full probability
        for all loop contexts. This check verifies that the sum of the individual loop types actually
        adds up to that obtained for any loop
        """
        ulength = 10
        data_split = {'ext': [], 'hp': [], 'int': [], 'mb': [] }
        data_all = []

        md = RNA.md()
        md.max_bp_span = 150
        md.window_size = 200
        fc = RNA.fold_compound(longseq, md, RNA.OPTION_WINDOW)

        fc.probs_window(ulength, RNA.PROBS_WINDOW_UP | RNA.PROBS_WINDOW_UP_SPLIT, up_split_callback, data_split)
        fc.probs_window(ulength, RNA.PROBS_WINDOW_UP, up_split_callback, data_all)

        # convert data to something we can more easily compare
        pfl_pu_ext = [ d['up'] for d in data_split['ext'] for i in range(0, len(shortseq) + 1) if d['i'] == i ]
        pfl_pu_hp = [ d['up'] for d in data_split['hp'] for i in range(0, len(shortseq) + 1) if d['i'] == i ]
        pfl_pu_int = [ d['up'] for d in data_split['int'] for i in range(0, len(shortseq) + 1) if d['i'] == i ]
        pfl_pu_mb = [ d['up'] for d in data_split['mb'] for i in range(0, len(shortseq) + 1) if d['i'] == i ]
        pfl_pu_tot = [ d['up'] for d in data_all for i in range(0, len(shortseq) + 1) if d['i'] == i ]

        # replace None values with 0
        pfl_pu_ext = [ [0 if v is None else v for v in d] for d in pfl_pu_ext ]
        pfl_pu_hp = [ [0 if v is None else v for v in d] for d in pfl_pu_hp ]
        pfl_pu_int = [ [0 if v is None else v for v in d] for d in pfl_pu_int ]
        pfl_pu_mb = [ [0 if v is None else v for v in d] for d in pfl_pu_mb ]
        pfl_pu_tot = [ [0 if v is None else v for v in d] for d in pfl_pu_tot ]

        # sum up unpaired probabilities of individual loop contexts
        # Note: (pfl_pu_* entries start at i = 1)
        pfl_pu_tot_sum = [ map(lambda x: x[0] + x[1] + x[2] + x[3], zip(pfl_pu_ext[i],pfl_pu_hp[i],pfl_pu_int[i],pfl_pu_mb[i])) for i in range(0, len(shortseq)) ]

        # compute difference between sum of individual loop contexts
        # and full probability
        pfl_pu_tot_diff = [ map(lambda x: abs(x[0] - x[1]), zip(pfl_pu_tot_sum[i], pfl_pu_tot[i])) for i in range(0, len(shortseq)) ]
        max_diff = max(map(max, zip(*pfl_pu_tot_diff)))

        self.assertTrue(max_diff < allowed_diff)


    def test_up_1_diff_pf(self):
        """Unpaired probabilities for u=1 are the same as obtained from global fold"""
        ulength = 1
        data = []

        md = RNA.md()
        md.max_bp_span = len(shortseq)
        md.window_size = len(shortseq)
        fc = RNA.fold_compound(shortseq, md, RNA.OPTION_WINDOW)

        fc.probs_window(ulength, RNA.PROBS_WINDOW_UP, up_split_callback, data)

        # compute pairing probability and from that
        # unpaired probability for individual nucleotides
        fc2 = RNA.fold_compound(shortseq)
        fc2.pf()
        bpp = fc2.bpp()

        up_single = [ 0 for i in range(0, len(shortseq) + 1) ]
        for i in range(1, len(shortseq) + 1):
            for j in range(i, len(shortseq) + 1):
                up_single[i] = up_single[i] + bpp[i][j]
                up_single[j] = up_single[j] + bpp[i][j]

        up_single = [ d for d in map(lambda x: 1 - x, up_single) ]
        # remove 0th element from list
        up_single.pop(0)

        # check diff between unpaired probs u = 1 and that retrieved from regular partition function
        pfl_pu_1 = [ d['up'][1] for d in data for i in range(0, len(shortseq) + 1) if d['i'] == i ]
        max_diff = max(map(lambda x: abs(x[0] - x[1]), zip(pfl_pu_1, up_single)))

        self.assertTrue(max_diff < allowed_diff)


    def test_up_exhaustive(self):
        """Unpaired probabilties for short RNAs are identical to what we get from exhaustive enumeration"""
        RNA.init_rand()
        randseq_length = 20
        ulength = 5

        md = RNA.md()
        md.max_bp_span = randseq_length
        md.window_size = randseq_length
        # turn off dangles due to scaling issues that prevent us
        # from comparing against exhaustive results
        md.dangles = 0

        # repreat the test 3 times
        for trial in range(0, 3):
            data_split = {'ext': [], 'hp': [], 'int': [], 'mb': [] }
            randseq = RNA.random_string(randseq_length, "ACGU")
            fc = RNA.fold_compound(randseq, md, RNA.OPTION_WINDOW)
            fc.probs_window(ulength, RNA.PROBS_WINDOW_UP | RNA.PROBS_WINDOW_UP_SPLIT, up_split_callback, data_split)

            # convert data to something we can more easily compare
            pfl_pu_ext = [ d['up'] for d in data_split['ext'] for i in range(0, randseq_length + 1) if d['i'] == i ]
            pfl_pu_hp = [ d['up'] for d in data_split['hp'] for i in range(0, randseq_length + 1) if d['i'] == i ]
            pfl_pu_int = [ d['up'] for d in data_split['int'] for i in range(0, randseq_length + 1) if d['i'] == i ]
            pfl_pu_mb = [ d['up'] for d in data_split['mb'] for i in range(0, randseq_length + 1) if d['i'] == i ]

            # replace None values with 0
            pfl_pu_ext = [ [0 if v is None else v for v in d] for d in pfl_pu_ext ]
            pfl_pu_hp = [ [0 if v is None else v for v in d] for d in pfl_pu_hp ]
            pfl_pu_int = [ [0 if v is None else v for v in d] for d in pfl_pu_int ]
            pfl_pu_mb = [ [0 if v is None else v for v in d] for d in pfl_pu_mb ]

            # compute unpaired probabilities from exhaustive enumeration
            pf = 0.0
            pu_ext = [ [ 0 for u in range(0, ulength + 1) ] for i in range(0, randseq_length + 1) ]
            pu_hp = [ [ 0 for u in range(0, ulength + 1) ] for i in range(0, randseq_length + 1) ]
            pu_int = [ [ 0 for u in range(0, ulength + 1) ] for i in range(0, randseq_length + 1) ]
            pu_mb = [ [ 0 for u in range(0, ulength + 1) ] for i in range(0, randseq_length + 1) ]

            # 1st, compute partition functions for loop contexts using subopt()
            fc = RNA.fold_compound(randseq, md)
            for s in fc.subopt(5000):
                pf = pf + math.exp(-s.energy/kT)
                ss = RNA.db_to_element_string(s.structure)
                for u in range(1, ulength + 1):
                    for j in range(1, randseq_length + 1):
                        if j - u + 1 <= 0:
                            continue
                        
                        is_ext = 1
                        for i in range(j - u + 1, j + 1):
                            if ss[i - 1] != 'e':
                                is_ext = 0
                                break
                        if is_ext:
                            pu_ext[j][u] = pu_ext[j][u] + math.exp(-s.energy/kT)
                        
                        is_hp = 1
                        for i in range(j - u + 1, j + 1):
                            if ss[i - 1] != 'h':
                                is_hp = 0
                                break
                        if is_hp:
                            pu_hp[j][u] = pu_hp[j][u] + math.exp(-s.energy/kT)
                        
                        is_int = 1
                        for i in range(j - u + 1, j + 1):
                            if ss[i - 1] != 'i':
                                is_int = 0
                                break
                        if is_int:
                            pu_int[j][u] = pu_int[j][u] + math.exp(-s.energy/kT)
                        
                        is_mb = 1
                        for i in range(j - u + 1, j + 1):
                            if ss[i - 1] != 'm':
                                is_mb = 0
                                break
                        if is_mb:
                            pu_mb[j][u] = pu_mb[j][u] + math.exp(-s.energy/kT)

            # 2nd, get unpaired probabilities
            for i in range(1, randseq_length + 1):
                for u in range(1, ulength + 1):
                    pu_ext[i][u] = pu_ext[i][u] / pf
                    pu_hp[i][u] = pu_hp[i][u] / pf
                    pu_int[i][u] = pu_int[i][u] / pf
                    pu_mb[i][u] = pu_mb[i][u] / pf

            pu_ext.pop(0)
            pu_hp.pop(0)
            pu_int.pop(0)
            pu_mb.pop(0)

            diff_ext = [ list(map(lambda x: abs(x[0] - x[1]), zip(pu_ext[i], pfl_pu_ext[i]))) for i in range(0, randseq_length) ]
            diff_hp = [ list(map(lambda x: abs(x[0] - x[1]), zip(pu_hp[i], pfl_pu_hp[i]))) for i in range(0, randseq_length) ]
            diff_int = [ list(map(lambda x: abs(x[0] - x[1]), zip(pu_int[i], pfl_pu_int[i]))) for i in range(0, randseq_length) ]
            diff_mb = [ list(map(lambda x: abs(x[0] - x[1]), zip(pu_mb[i], pfl_pu_mb[i]))) for i in range(0, randseq_length) ]

            # flatten lists
            diff_ext = [ d for dd in diff_ext for d in dd ]
            diff_hp = [ d for dd in diff_hp for d in dd ]
            diff_int = [ d for dd in diff_int for d in dd ]
            diff_mb = [ d for dd in diff_mb for d in dd ]

            # compute average difference
            diff_ext_avg = sum(diff_ext) / float(len(diff_ext))
            diff_hp_avg = sum(diff_hp) / float(len(diff_hp))
            diff_int_avg = sum(diff_int) / float(len(diff_int))
            diff_mb_avg = sum(diff_mb) / float(len(diff_mb))

            # compute variance
            diff_ext_var = sum(map(lambda x: (x - diff_ext_avg)**2, diff_ext)) / float(len(diff_ext))
            diff_hp_var = sum(map(lambda x: (x - diff_hp_avg)**2, diff_hp)) / float(len(diff_hp))
            diff_int_var = sum(map(lambda x: (x - diff_int_avg)**2, diff_int)) / float(len(diff_int))
            diff_mb_var = sum(map(lambda x: (x - diff_mb_avg)**2, diff_mb)) / float(len(diff_mb))

            print("\tTrial %d - Difference between pflfold unpaired probs and exhaustive enumeration" % trial)
            print("\tExterior loop    (avg, var)\t%g\t%g" % (diff_ext_avg, diff_ext_var))
            print("\tHairpin loop     (avg, var)\t%g\t%g" % (diff_hp_avg, diff_hp_var))
            print("\tInterior loop    (avg, var)\t%g\t%g" % (diff_int_avg, diff_int_var))
            print("\tMultibranch loop (avg, var)\t%g\t%g" % (diff_mb_avg, diff_mb_var))

            self.assertTrue(diff_ext_avg < allowed_diff)
            self.assertTrue(diff_hp_avg < allowed_diff)
            self.assertTrue(diff_int_avg < allowed_diff)
            self.assertTrue(diff_mb_avg < allowed_diff)



if __name__ == '__main__':
    pf_window_functionTest.DATADIR = RNApath.getDataDirPath()
    unittest.main(testRunner=taprunner.TAPTestRunner())
