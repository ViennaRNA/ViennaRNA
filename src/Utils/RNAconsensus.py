#!/usr/bin/env  python3
#
# 2024, Ronny Lorenz <ronny@tbi.univie.ac.at>
#
#


# dictionary of canonical base pairs
BPs = {
  "AU" : 5,
  "GC" : 1,
  "CG" : 2,
  "UA" : 6,
  "GU" : 3,
  "GT" : 3,
  "TG" : 4,
  "UG" : 4,
  "AT" : 5,
  "TA" : 6
}


def read_alignment(filename):
    alignment = None

    try:
        import RNA
        ret, names, seqs, id, structure = RNA.file_msa_read(filename, RNA.FILE_FORMAT_MSA_DEFAULT | RNA.FILE_FORMAT_MSA_QUIET)
        if ret > 0:
            alignment = [ (name, seq) for (name, seq) in zip(names, seqs) ]
    except:
        print("Failed to parse alignment file")

    return alignment


def structure_from_alifold(filename):
    import re

    structure     = None
    structure_pat = re.compile(r"^([\.\(\)\+]+)")

    with open(filename) as f:
        for line in f.readlines():
            line = line.strip()
            m = structure_pat.match(line)
            if m:
                structure = m.group(1)
                break

    return structure


def structure_from_dotplot(filename, threshold = 0.9):
    structure = None
    
    return structure


def map_structure_to_seqs(alignment, structure, turn = 3):
    pt = None
    has_vrna = False
    try:
        import RNA
        pt = RNA.ptable(structure)
        has_vrna = True
    except:
        pt = make_pair_table(structure)

    for name, seq in alignment:
        seql = list("-" + seq)
        cons = list("x" + structure)

        for i in range(1, pt[0] + 1):
            if seql[i] == '-':
                # mark position for removal
                cons[i] = 'x'
                # if this is an opening base pair, make the pairing partner unpaired
                if pt[i] > i:
                    cons[pt[i]] = '.'
            elif pt[i] > i:
                # check if this is an allowed base pair in the current sequence
                if seql[i] + seql[pt[i]] not in BPs:
                    cons[i] = "."
                    cons[pt[i]] = "."

        seql = "".join(seql).replace("-", "")
        cons = "".join(cons).replace("x", "")

        # enforce minimum hairpin length constraint
        pts = None
        if has_vrna:
            pts = RNA.ptable(cons)
        else:
            pts = make_pair_table(cons)

        for i in range(1, pts[0] + 1):
            if pts[i] > i and (pts[i] - i - 1) < turn:
                pts[pts[i]] = 0
                pts[i] = 0

        # convert back
        cons = [ "." for _ in range(pts[0] + 1)]
        for i in range(1, pts[0] + 1): 
            if pts[i] > i:
                cons[i - 1] = '('
                cons[pts[i] - 1] = ')'

        cons = "".join(cons)

        print(f">{name}\n{seql}\n{cons}")


def make_pair_table(structure):
    """
    create a pair table for an input structure in dot-bracket format
    """
    pt    = [0 for _ in range(len(structure) + 1)]
    pt[0] = len(structure)
    stack = []

    for j, c in enumerate(list(structure), start = 1):
        if c == '(':
            stack.append(j)
        elif c == ')':
            if stack:
                i = stack.pop()
            else:
                print("unbalanced brackets in consensus structure")
                return None

            pt[j] = i
            pt[i] = j

    if stack:
        print("unbalanced brackets in consensus structure")
        return None

    return pt


def refold(args):
    """
    Legacy refold.pl functionality
    """
    # try parsing the input alignment to work on
    aln = None
    
    if args.alignment:
        aln = read_alignment(args.alignment)
    elif args.alnfile:
        aln = read_alignment(args.alnfile)

    # only continue if we successfully read an input alignment
    if aln != None:
        # read consensus structure information, either from
        # RNAalifold output or dot-plot. 
        # Try dotplot first...
        structure = None
        if args.dotplot:
            structure = structure_from_dotplot(args.dotplot, args.threshold)
        elif args.structfile:
            structure = structure_from_dotplot(args.structfile, args.threshold)

        # Try RNAalifold ouput if dot-plot parsing was unsuccessful
        if structure == None:
            if args.alifold_output:
                structure = structure_from_alifold(args.alifold_output)
            elif args.structfile:
                structure = structure_from_alifold(args.structfile)

        if structure:
            map_structure_to_seqs(aln, structure, args.turn)

            return 0 # Success

    return 1  # Failure


def main():
    import sys
    import argparse

    def set_default_subparser(parser,
                              default_subparser):
        """default subparser selection. Call after setup, just before parse_args()

        parser: the name of the parser you're making changes to
        default_subparser: the name of the subparser to call by default"""
        subparser_found = False
        for arg in sys.argv[1:]:
            if arg in ['-h', '--help']:  # global help if no subparser
                break
        else:
            for x in parser._subparsers._actions:
                if not isinstance(x, argparse._SubParsersAction):
                    continue
                for sp_name in x._name_parser_map.keys():
                    if sp_name in sys.argv[1:]:
                        subparser_found = True
            if not subparser_found:
                # insert default in first position before all other arguments
                sys.argv.insert(1, default_subparser) 


    parser = argparse.ArgumentParser(
        prog = "RNAconsens",
        description = "An program to predict RNA secondary structures for single sequences based on the information gained from a multiple sequence alignment of homologous sequences",
        epilog = "If in doubt our program is right, nature is at fault.  Comments should be sent to rna@tbi.univie.ac.at.")

    # global options for all modes
    parser.add_argument("-a", "--alignment", metavar = "filename", type=str, help="A multiple sequence alignment file name")
    parser.add_argument("--turn", type=int, help="", default = 3)

    # create subparsers for the different operating modes
    subparsers = parser.add_subparsers(title='Modes',
                                       description=None,
                                       help='available modes')

    # 1. Legacy refold.pl mode where a structure prediction from RNAalifold is taken as
    # input to create hard constraints for predictions of structures for the single
    # sequences within the alignment. For that purpose, the program prints each
    # sequence and the consensus structure mapped to the sequence such that the
    # output can be used as input for RNAfold.
    parser_a = subparsers.add_parser('refold', help='The legacy refold.pl mode')

    group1 = parser_a.add_argument_group('Consensus MFE', 'Map the consensus MFE structure as predicted by RNAalifold to each sequence in the input alignment')
    group1.add_argument("-f", "--alifold-output", type=str, help="The output of RNAalifold for the input alignment", default = None)

    group2 = parser_a.add_argument_group('Consensus Probabilities', 'Create structure constraints for each sequence in the input alignment from consensus base pair probabilities as predicted by RNAalifold')
    group2.add_argument("-d", "--dotplot", metavar = "filename", type=str, help="The dot-plot created by RNAalifold for the input alignment", default = None)
    group2.add_argument("-t", "--threshold", type=float, help="", default = 0.9)

    parser_a.add_argument("alnfile", nargs='?', metavar = "alignment_file")
    parser_a.add_argument("structfile", nargs='?', metavar = "dotplot_or_alifold_file")
    parser_a.set_defaults(func=refold)

    set_default_subparser(parser, "refold")

    args = parser.parse_args()

    # call the function corresponding to the chosen mode
    # and use its return value as exit code
    sys.exit(args.func(args))



if __name__ == '__main__':
    main()
