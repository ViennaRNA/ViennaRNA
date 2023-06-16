import RNApath

RNApath.addSwigInterfacePath()


import RNA
import unittest
from py_include import taprunner
from struct import *
import locale
import json
import re

locale.setlocale(locale.LC_ALL, 'C')

def extract_duplex_data(duplex, data, one_letter_code, fallback):
    e  = data['dG37']
    s1 = duplex[:data['length1']]
    s2 = duplex[data['length1']:]
    s2 = s2[::-1]
    mpos = [m.start() + 1 for m in re.finditer(one_letter_code, s1 + s2)]
    seq = s1 + "&" + s2
    seq = seq.replace(one_letter_code, fallback)

    return e, seq, mpos


class ModifiedBaseTests(unittest.TestCase):
    """Some tests for modified bases"""

    def test_pseudouridine_duplexes(self):
        """Pseudouridine duplex energies"""
        params = json.loads(RNA.parameter_set_rna_mod_pseudouridine_parameters)
        duplexes = params['duplexes']
        code     = params['modified_base']['one_letter_code']
        fallback = params['modified_base']['fallback']
        for k,v in duplexes.items():
            e, seq, mpos = extract_duplex_data(k, v, code, fallback)
            # prepare computation
            fc = RNA.fold_compound(seq)
            fc.sc_mod_pseudouridine(mpos)
            e_diff = abs(fc.mfe()[1] - e)
            self.assertTrue(e_diff < 0.7)


if __name__ == '__main__':
    unittest.main(testRunner=taprunner.TAPTestRunner())

