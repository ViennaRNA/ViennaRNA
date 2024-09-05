import unittest

if __name__ == '__main__':
    from py_include import taprunner
    import RNApath
    RNApath.addSwigInterfacePath()

import RNA

seq = "CCCGCCCCACCACUCCAGCAGACCUUGCCCCUUGUGAGCUGGAUAGACUUGGGAUGGGGAGGGAGGGAGUUUUGUCUGUCUCCCUCCCCUCUCAGAACAUACUGAUUGGGAGGUGCG"

def t(x, i, n):
    """
    Circular rotation of sequence position x by i positions for sequence of length n
    """
    return ((x + i) - 1) % n + 1

def c(x, y, i, n):
    """
    Circular rotation of base pair positions (x,y) by i positions for sequence of length m
    """
    xx, yy = t(x, i, n), t(y, i, n)
    return (xx, yy) if yy > xx else (yy, xx)


def rotate_structure(db, i):
    """
    Circular rotation of a dot-bracket structure by i positions
    
    Simply rotating the db string like:
    
    db = db[n-i:] + db[:n - i]

    doesn't work sinze some base pairs may end up as ) ( pair
    instead of ( ). So, we transform the structure at pair
    table level and convert it back to dot-bracket for later
    comparison.
    
    In addition, any non-dot-bracket characters will be re-inserted
    as they are not present at the pair table level
    """
    n = len(db)

    if (i % n) == 0:
        return str(db)

    pt_db = [ e for e in RNA.ptable(db) ]
    pt    = [ 0 for _ in range(n + 1) ]
    pt[0] = n

    for j, p in enumerate(pt_db[1:], start = 1):
        if p > j:
            ii, jj = c(j, p, i, n)
            pt[jj], pt[ii] = ii, jj

    ss = list(RNA.db_from_ptable(pt))

    # re-insert any special characters, such as gquads
    for j, cc in enumerate(list(db), start = 1):
        if cc not in ["(", ")", "."]:
            ss[t(j, i, n) - 1] = cc

    return "".join(ss)


def rotate_pair_probs(plist, i, n, ptype = RNA.PLIST_TYPE_BASEPAIR):
    """
    Circular rotation of a pair probability list by i positions
    """
    return sorted([ (*c(p.i, p.j, i, n), f"{p.p:.10f}") for p in plist if p.type == ptype ])


class CircMFETest(unittest.TestCase):
    def test_mfe(self):
        """Circular RNA MFE and structure"""

        md        = RNA.md()
        md.circ   = 1
        n         = len(seq)
        fc        = RNA.fold_compound(seq, md)
        ss0, mfe0 = fc.mfe()

        for i in range(1, n):
            # rotate sequence by i
            s = seq[i:] + seq[:i]

            fc = RNA.fold_compound(s, md)
            ss, mfe = fc.mfe()

            # rotate structure by i
            ss = rotate_structure(ss, i)

            self.assertEqual(mfe, mfe0)
            self.assertEqual(ss, ss0)

    def test_mfe_gquad(self):
        """Circular RNA MFE and structure (G-Quadruplexes)"""

        md        = RNA.md()
        md.circ   = 1
        md.gquad  = 1
        n         = len(seq)
        fc        = RNA.fold_compound(seq, md)
        ss0, mfe0 = fc.mfe()

        for i in range(1, n):
            # rotate sequence by i
            s = seq[i:] + seq[:i]

            fc = RNA.fold_compound(s, md)
            ss, mfe = fc.mfe()

            # rotate structure by i
            ss = rotate_structure(ss, i)

            self.assertEqual(mfe, mfe0)
            self.assertEqual(ss, ss0)


class CircPartFunTest(unittest.TestCase):
    def test_pf(self):
        """Circular RNA partition function and pair probabilities"""
        cutoff    = 1e-5

        md        = RNA.md()
        md.circ   = 1
        n         = len(seq)
        fc        = RNA.fold_compound(seq, md)
        ss0, mfe0 = fc.mfe()
        fc.exp_params_rescale(mfe0)
        pp0, dG0  = fc.pf()
        pairs0    = rotate_pair_probs(fc.plist_from_probs(cutoff), 0, n)

        for i in range(1, n):
            # rotate sequence by i
            s = seq[i:] + seq[:i]

            fc      = RNA.fold_compound(s, md)
            ss, mfe = fc.mfe()
            fc.exp_params_rescale(mfe)
            pp, dG  = fc.pf()
            pairs   = rotate_pair_probs(fc.plist_from_probs(cutoff), i, n)

            # rotate structure by i
            self.assertEqual(dG, dG0)

            # rotate pairing propensity string
            # pp = pp[n-i:] + pp[:n - i]
            # this doesn't always work, and for now, we don't have
            # a sane way to transform similar to MFE structure above
            #self.assertEqual(pp, pp0)

            # compare base pair probabilities
            self.assertEqual(pairs, pairs0)

    def test_pf_gquad(self):
        """Circular RNA partition function and pair probabilities (G-Quadruplexes)"""
        cutoff    = 1e-5

        md        = RNA.md()
        md.circ   = 1
        md.gquad  = 1
        n         = len(seq)
        fc        = RNA.fold_compound(seq, md)
        ss0, mfe0 = fc.mfe()
        fc.exp_params_rescale(mfe0)
        pp0, dG0  = fc.pf()
        pairs0    = rotate_pair_probs(fc.plist_from_probs(cutoff), 0, n)
        qgs0      = rotate_pair_probs(fc.plist_from_probs(cutoff), 0, n, RNA.PLIST_TYPE_GQUAD)

        for i in range(1, n):
            # rotate sequence by i
            s = seq[i:] + seq[:i]

            fc      = RNA.fold_compound(s, md)
            ss, mfe = fc.mfe()
            fc.exp_params_rescale(mfe)
            pp, dG  = fc.pf()
            pairs   = rotate_pair_probs(fc.plist_from_probs(cutoff), i, n)
            qgs     = rotate_pair_probs(fc.plist_from_probs(cutoff), i, n, RNA.PLIST_TYPE_GQUAD)

            # rotate structure by i
            self.assertEqual(dG, dG0)

            # rotate pairing propensity string
            # pp = pp[n-i:] + pp[:n - i]
            # this doesn't always work, and for now, we don't have
            # a sane way to transform similar to MFE structure above
            #self.assertEqual(pp, pp0)

            # compare base pair probabilities
            self.assertEqual(pairs, pairs0)

            # compare gquad probabilities
            self.assertEqual(qgs, qgs0)


if __name__ == '__main__':
    unittest.main(testRunner=taprunner.TAPTestRunner())

