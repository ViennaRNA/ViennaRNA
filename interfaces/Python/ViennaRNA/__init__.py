try:
    from RNA import RNA
    from RNA.RNA import *
except Exception as E:
    raise ImportError('ViennaRNA / RNAlib not found')

ViennaRNA = RNA
viennarna = RNA
RNA       = RNA
rna       = RNA

def main():
    return RNA

if __name__ == '__main__':
    main()
