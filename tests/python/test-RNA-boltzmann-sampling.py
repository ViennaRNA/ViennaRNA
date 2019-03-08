import RNApath

RNApath.addSwigInterfacePath()

import RNA
import unittest

sequence = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUACUUAUACCGCCUGUGCGGACUACUAUCCUGACCACAUAGU"

def store_structure(s, data):
    if s:
        data.append(s)


def prepare_fc():
    md          = RNA.md()
    md.uniq_ML  = 1
    fc          = RNA.fold_compound(sequence, md)
    (ss, mfe)   = fc.mfe()

    fc.exp_params_rescale(mfe)
    fc.pf()

    return fc

class GeneralTests(unittest.TestCase):
    def test_pbacktrack5(self):
        print "test_pbacktrack    (single sub-structure, a.k.a. pbacktrack5)"
        fc = prepare_fc()

        s  = fc.pbacktrack(10)
        self.assertEqual(len(s), 10)

        s  = fc.pbacktrack(50)
        self.assertEqual(len(s), 50)

    def test_pbacktrack5_multi(self):
        print "test_pbacktrack    (multiple sub-structures, a.k.a. pbacktrack5_num)"
        fc = prepare_fc()

        ss  = fc.pbacktrack(10, 20)
        self.assertEqual(len(ss), 20)
        for s in ss:
            self.assertEqual(len(s), 10)

        ss  = fc.pbacktrack(50, 100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), 50)

    def test_pbacktrack(self):
        print "test_pbacktrack    (single structure, a.k.a. pbacktrack)"
        fc = prepare_fc()

        s  = fc.pbacktrack()
        self.assertEqual(len(s), len(sequence))

    def test_pbacktrack_multi(self):
        print "test_pbacktrack    (multiple structures, a.k.a. pbacktrack_num)"
        fc = prepare_fc()

        ss  = fc.pbacktrack(0, 100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), len(sequence))

    def test_pbacktrack_nr(self):
        print "test_pbacktrack_nr"
        fc = prepare_fc()

        ss  = fc.pbacktrack_nr(100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), len(sequence))
        # check for uniqueness, i.e. no duplicates
        sss = list(set(ss))
        self.assertEqual(len(ss), len(sss))

    def test_pbacktrack_nr_resume(self):
        print "test_pbacktrack_nr_resume"
        num_samples = 500
        iterations  = 5
        fc      = prepare_fc()
        d       = None
        s_list  = list()
        for i in range(0, iterations):
            d, ss   = fc.pbacktrack_nr_resume(num_samples, d)
            s_list  = s_list + list(ss)

        self.assertEqual(len(s_list), iterations * num_samples)

        for s in s_list:
            self.assertEqual(len(s), len(sequence))

        # check for uniqueness, i.e. no duplicates
        sss = list(set(s_list))
        self.assertEqual(len(s_list), len(sss))

    def test_pbacktrack5_cb(self):
        print "test_pbacktrack_cb (multiple sub-structures, a.k.a. pbacktrack5_cb)"
        fc = prepare_fc()
        ss = list()
        i = fc.pbacktrack_cb(10, 20, store_structure, ss)
        self.assertEqual(i, 20)
        self.assertEqual(len(ss), 20)
        for s in ss:
            self.assertEqual(len(s), 10)

        ss = list()
        i = fc.pbacktrack_cb(50, 100, store_structure, ss)
        self.assertEqual(i, 100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), 50)

    def test_pbacktrack_cb(self):
        print "test_pbacktrack_cb (multiple structures, a.k.a. pbacktrack_cb)"
        fc = prepare_fc()
        ss = list()
        i = fc.pbacktrack_cb(0, 100, store_structure, ss)
        self.assertEqual(i, 100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), len(sequence))

    def test_pbacktrack_nr_cb(self):
        print "test_pbacktrack_nr_cb"
        fc = prepare_fc()
        ss = list()
        i = fc.pbacktrack_nr_cb(100, store_structure, ss)
        self.assertEqual(i, 100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), len(sequence))

        # check for uniqueness, i.e. no duplicates
        sss = list(set(ss))
        self.assertEqual(len(ss), len(sss))

    def test_pbacktrack_nr_resume_cb(self):
        print "test_pbacktrack_nr_resume_cb"
        num_samples = 500
        iterations  = 5
        fc = prepare_fc()
        ss = list()
        d  = None
        for i in range(0, iterations):
            d, i = fc.pbacktrack_nr_resume_cb(num_samples, store_structure, ss, d)
            self.assertEqual(i, num_samples)

        self.assertEqual(len(ss), iterations * num_samples)
        for s in ss:
            self.assertEqual(len(s), len(sequence))

        # check for uniqueness, i.e. no duplicates
        sss = list(set(ss))
        self.assertEqual(len(ss), len(sss))


if __name__ == '__main__':
    unittest.main();
