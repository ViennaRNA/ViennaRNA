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
allowed_diff = 1e-8

# some sequences to work with
shortseq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUA"
longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"



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
                if v != None:
                    alocal.append(RNA.ep(j + 1 - k, j, v, RNA.PLIST_TYPE_UNPAIRED))

        # compute accessibility (global implementation)
        fc = RNA.fold_compound(s)
        ss, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        a = fc.eqm_probs_up(max_u = ulength, cutoff = 0)

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
