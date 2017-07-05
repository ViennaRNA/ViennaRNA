import RNApath

RNApath.addSwigInterfacePath(3)


import RNA
import unittest

longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU"

def pf_window_callback(v, v_size, i, maxsize, what, data=None):
    if what & RNA.PROBS_WINDOW_UP:
        data['up'].append({ 'i': i, 'up': v})
    else:
        data['bpp'].extend([{'i': i, 'j': j, 'p': p} for j, p in enumerate(v) if (p is not None) and (p >= 0.01)])


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


if __name__ == '__main__':
    unittest.main()
