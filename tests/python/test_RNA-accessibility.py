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


# maximum allowed difference beteen compared probabilties
allowed_diff  = 1e-8
cutoff        = 1e-10

# some sequences to work with
shortseq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUA"
longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"

aln = [ "CAAUCUUAACGAUGGAUGUCUUGGUUCCUAUAGCGAUGAAGGCCGCAGCAAAGUGCGAUAUGCAAUGA-----AAAAUG--CA-AUUA--------CUGUGAAUCAUCAGAAUGCUGAAUGUAAACUAUACCAUAUUUACCC----UUAUGGGCAAAUUAACGUGGUAUUCC-UACAGAA-----A",
        "AACUUUCAGCGACGGAUGUCUCGGCUCGAACAACGAUGAAGGGCGCAGCGAAAUGUGAUAAGCAAUGU-----GAAUUG--CAGAAUU--------CCGUGAACCAAUAGGGACUUGAACGUAUACUGCGCUUUCGGGA-------UAUC--------CCUGAAGGCAUGCC-UGCUUCAGUGUCU",
        "AACCCU-AGUGAUGGAUGUCUAGGCUCCCGUAUCGAUGAAGAACGUAGCAGAAAGCGACAUGCGAUGAUACGCGCAUAGACCAGGAUUCCGCGGACUCGUGAAUCAUUGAAUCUUCGAACGCAUCUGGCGCUCAGCAGUGGAGCCGUAGC--------AGGCUGAGCAUGCCCUGGAGUAGUGACA"
      ]

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


class bpp_functionTest(unittest.TestCase):
    def test_bpp_consistency(self):
        """Base pair probabilities from callback vs. dynamic programming table entries"""
        s = shortseq

        fc = RNA.fold_compound(s)
        ss, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        bpp_dp  = sorted(sorted(fc.plist_from_probs(cutoff = cutoff), key = lambda ep: ep.j ), key=lambda ep: ep.i)
        bpp_cb  = sorted(sorted(fc.eqm_probs_bpp(cutoff = cutoff), key = lambda ep: ep.j ), key=lambda ep: ep.i)

        self.assertTrue(len(bpp_dp) == len(bpp_cb))

        for pair in zip(bpp_dp, bpp_cb):
            self.assertEqual(pair[0].i, pair[1].i)
            self.assertEqual(pair[0].j, pair[1].j)
            self.assertTrue(abs(pair[0].p - pair[1].p) < allowed_diff)


    def test_bpp_consistency(self):
        """Base pair probabilities (circular RNA) from callback vs. dynamic programming table entries"""
        s = shortseq

        fc = RNA.fold_compound(s, RNA.md(circ = 1))
        ss, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        bpp_dp  = sorted(sorted(fc.plist_from_probs(cutoff = cutoff), key = lambda ep: ep.j ), key=lambda ep: ep.i)
        bpp_cb  = sorted(sorted(fc.eqm_probs_bpp(cutoff = cutoff), key = lambda ep: ep.j ), key=lambda ep: ep.i)

        self.assertTrue(len(bpp_dp) == len(bpp_cb))

        for pair in zip(bpp_dp, bpp_cb):
            self.assertEqual(pair[0].i, pair[1].i)
            self.assertEqual(pair[0].j, pair[1].j)
            self.assertTrue(abs(pair[0].p - pair[1].p) < allowed_diff)


    def test_bpp_consistency_comparative(self):
        """Consensus base pair probabilities from callback vs. dynamic programming table entries"""
        fc = RNA.fold_compound(aln)
        ss, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        bpp_dp  = sorted(sorted(fc.plist_from_probs(cutoff = cutoff), key = lambda ep: ep.j ), key=lambda ep: ep.i)
        bpp_cb  = sorted(sorted(fc.eqm_probs_bpp(cutoff = cutoff), key = lambda ep: ep.j ), key=lambda ep: ep.i)

        self.assertTrue(len(bpp_dp) == len(bpp_cb))

        for pair in zip(bpp_dp, bpp_cb):
            self.assertEqual(pair[0].i, pair[1].i)
            self.assertEqual(pair[0].j, pair[1].j)
            self.assertTrue(abs(pair[0].p - pair[1].p) < allowed_diff)


class accessibility_functionTest(unittest.TestCase):
    def test_up_consistency(self):
        """Unpaired probabilities from local vs. global prediction"""
        s = longseq
        n = len(s)
        ulength = n
        data = []

        # compute accessibility (local implementation)
        md = RNA.md()
        md.max_bp_span = n
        md.window_size = n
        fc = RNA.fold_compound(s, md, RNA.OPTION_WINDOW)
        fc.probs_window(ulength, RNA.PROBS_WINDOW_UP, up_callback, data)

        # transform accessibilities obtained from local folding int
        # a list of vrna_ep_t entries
        alocal = []

        for d in data:
            j = d['i']
            a = d['up']
            for k, v in enumerate(a[1:], start = 1):
                if v != None and v >= cutoff:
                    alocal.append(RNA.ep(j + 1 - k, j, v, RNA.PLIST_TYPE_UNPAIRED))

        # compute accessibility (global implementation)
        fc = RNA.fold_compound(s)
        ss, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        a = fc.eqm_probs_up(max_u = ulength, cutoff = cutoff)

        # sort values
        alocal = sorted(sorted(alocal, key = lambda ep: ep.j ), key=lambda ep: ep.i)
        a = sorted(sorted(a, key = lambda ep: ep.j ), key = lambda ep: ep.i)

        # compare results
        self.assertTrue(len(alocal) == len(a))

        for pair in zip(alocal, a):
            self.assertEqual(pair[0].i, pair[1].i)
            self.assertEqual(pair[0].j, pair[1].j)
            self.assertTrue(abs(pair[0].p - pair[1].p) < allowed_diff)


if __name__ == '__main__':
    accessibility_functionTest.DATADIR = RNApath.getDataDirPath()
    unittest.main(testRunner=taprunner.TAPTestRunner())
