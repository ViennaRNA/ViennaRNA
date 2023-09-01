import unittest
import sys

if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA


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
        """Boltzmann sampling (single sub-structure, a.k.a. pbacktrack5)"""
        fc = prepare_fc()

        s  = fc.pbacktrack5(10)
        self.assertEqual(len(s), 10)

        s  = fc.pbacktrack5(50)
        self.assertEqual(len(s), 50)

    def test_pbacktrack5_multi(self):
        """Boltzmann sampling (multiple sub-structures, a.k.a. pbacktrack5_num)"""
        fc = prepare_fc()

        ss  = fc.pbacktrack5(20, 10)
        self.assertEqual(len(ss), 20)
        for s in ss:
            self.assertEqual(len(s), 10)

        ss  = fc.pbacktrack5(100, 50)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), 50)

    def test_pbacktrack(self):
        """Boltzmann sampling (single structure, a.k.a. pbacktrack)"""
        fc = prepare_fc()

        s  = fc.pbacktrack()
        self.assertEqual(len(s), len(sequence))

    def test_pbacktrack_multi(self):
        """Boltzmann sampling (multiple structures, a.k.a. pbacktrack_num)"""
        fc = prepare_fc()

        ss  = fc.pbacktrack(100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), len(sequence))

    @unittest.skipIf(sys.platform.startswith("win"), "May fail under Windows")
    def test_pbacktrack_nr(self):
        """Boltzmann sampling - non redundant"""
        try:
            import pytest
            if sys.platform.startswith("win"):
                pytest.skip("May fail under Windows")
        except:
            pass

        fc = prepare_fc()

        ss  = fc.pbacktrack(100, RNA.PBACKTRACK_NON_REDUNDANT)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), len(sequence))
        # check for uniqueness, i.e. no duplicates
        sss = list(set(ss))
        self.assertEqual(len(ss), len(sss))

    @unittest.skipIf(sys.platform.startswith("win"), "May fail under Windows")
    def test_pbacktrack_nr_resume(self):
        """Boltzmann sampling - non redundant resume"""
        try:
            import pytest
            if sys.platform.startswith("win"):
                pytest.skip("May fail under Windows")
        except:
            pass

        num_samples = 500
        iterations  = 15
        fc      = prepare_fc()
        d       = None
        s_list  = list()
        for _ in range(0, iterations):
            d, ss   = fc.pbacktrack(num_samples, d, RNA.PBACKTRACK_NON_REDUNDANT)
            s_list  = s_list + list(ss)

        self.assertEqual(len(s_list), iterations * num_samples)

        for s in s_list:
            self.assertEqual(len(s), len(sequence))

        # check for uniqueness, i.e. no duplicates
        sss = list(set(s_list))
        self.assertEqual(len(s_list), len(sss))

    def test_pbacktrack5_cb(self):
        """Boltzmann sampling (multiple sub-structures, a.k.a. pbacktrack5_cb)"""
        fc = prepare_fc()
        ss = list()
        i = fc.pbacktrack5(20, 10, store_structure, ss)
        self.assertEqual(i, 20)
        self.assertEqual(len(ss), 20)
        for s in ss:
            self.assertEqual(len(s), 10)

        ss = list()
        i = fc.pbacktrack5(100, 50, store_structure, ss)
        self.assertEqual(i, 100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), 50)

    def test_pbacktrack_cb(self):
        """Boltzmann sampling (multiple structures, a.k.a. pbacktrack_cb)"""
        fc = prepare_fc()
        ss = list()
        i = fc.pbacktrack(100, store_structure, ss)
        self.assertEqual(i, 100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), len(sequence))

    @unittest.skipIf(sys.platform.startswith("win"), "May fail under Windows")
    def test_pbacktrack_nr_cb(self):
        """Boltzmann sampling - non redundant - callback"""
        try:
            import pytest
            if sys.platform.startswith("win"):
                pytest.skip("May fail under Windows")
        except:
            pass

        fc = prepare_fc()
        ss = list()
        i = fc.pbacktrack(100, store_structure, ss, RNA.PBACKTRACK_NON_REDUNDANT)
        self.assertEqual(i, 100)
        self.assertEqual(len(ss), 100)
        for s in ss:
            self.assertEqual(len(s), len(sequence))

        # check for uniqueness, i.e. no duplicates
        sss = list(set(ss))
        self.assertEqual(len(ss), len(sss))

    @unittest.skipIf(sys.platform.startswith("win"), "May fail under Windows")
    def test_pbacktrack_nr_resume_cb(self):
        """Boltzmann sampling - non redundant resume - callback"""
        try:
            import pytest
            if sys.platform.startswith("win"):
                pytest.skip("May fail under Windows")
        except:
            pass

        num_samples = 500
        iterations  = 15
        fc = prepare_fc()
        ss = list()
        d  = None
        for _ in range(0, iterations):
            d, i = fc.pbacktrack(num_samples, store_structure, ss, d, RNA.PBACKTRACK_NON_REDUNDANT)
            self.assertEqual(i, num_samples)

        self.assertEqual(len(ss), iterations * num_samples)
        for s in ss:
            self.assertEqual(len(s), len(sequence))

        # check for uniqueness, i.e. no duplicates
        sss = list(set(ss))
        self.assertEqual(len(ss), len(sss))


if __name__ == '__main__':
    unittest.main(testRunner=taprunner.TAPTestRunner())
