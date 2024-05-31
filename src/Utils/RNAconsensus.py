#!/usr/bin/env  python3
#
# 2024, Ronny Lorenz <ronny@tbi.univie.ac.at>
#
#

import sys


# dictionary of canonical base pairs
_BPs = {
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


def _parse_clustal(filename):
    """
    Helper function that parses a ClustalW or Stockholm formatted
    MSA file and returns the alignment as list of pairs of name
    and aligned sequence
    """
    import re
    alignment = dict()
    order     = 0

    with open(filename) as f:
        for line in f.readlines():
            line = line.strip()
            # stop processing if stockholm file entry ends
            if line.startswith("//"):
                break

            split = line.split()

            # only consider lines where we have two separate columns of data
            # and the second column looks more or less like a sequence
            if len(split) == 2:
                if re.match(r'^[a-zA-Z\-\_]+$', split[1]):
                    if split[0] in alignment:
                        alignment[split[0]][0] += split[1]
                    else:
                        alignment[split[0]] = (split[1], order)
                        order += 1

    if len(alignment) > 0:
        # reformat alignment and preserver order
        return [ (k, v[0]) for k, v in sorted(alignment.items(), key=lambda item: item[1][1]) ]
    else:
        return None


def _read_alignment(args):
    """
    Wrapper that tries reading an alignment file via the ViennaRNA Python interface
    and falls-back to the _parse_clustal() function in this script if the Python interface
    is not accessible
    """
    alignment = None

    filename = None

    if args.alignment:
        filename = args.alignment
    elif args.alnfile:
        filename = args.alnfile

    try:
        import RNA
        ret, names, seqs, id, structure = RNA.file_msa_read(filename, RNA.FILE_FORMAT_MSA_DEFAULT | RNA.FILE_FORMAT_MSA_QUIET)
        if ret > 0:
            alignment = [ (name, seq) for (name, seq) in zip(names, seqs) ]
    except:
        alignment = _parse_clustal(filename)

    return alignment


def _structure_from_alifold(args):
    """
    A helper function to parse RNAalifold output and extract the predicted
    MFE consensus structrue
    """
    import re

    filename      = None
    structure     = None
    structure_pat = re.compile(r"^([\.\(\)\+]+)")

    if args.alifold_output:
        filename = args.alifold_output
    elif args.structfile:
        filename = args.structfile

    try:
        with open(filename) as f:
            for line in f.readlines():
                line = line.strip()
                m = structure_pat.match(line)
                if m:
                    structure = m.group(1)
                    break
    except:
        pass

    return structure


def _hard_constraints_output(name, sequence, constraint, args):
    """
    Default callback for the hard constraints prediction that
    either runs the predictions for the single sequences or
    simply prints the constraint strings.
    """
    if args.run:
        ss  = None
        mfe = None

        # Try using the ViennaRNA Python API first
        try:
            import RNA
            sequence  = sequence.replace("T", "U")
            sequence  = sequence.replace("t", "u")
            fc        = RNA.fold_compound(sequence)
            options   = RNA.CONSTRAINT_DB_DEFAULT

            if args.enforce:
                options |= RNA.CONSTRAINT_DB_ENFORCE_BP

            fc.hc_add_from_db(constraint, options)

            ss, mfe = fc.mfe()

            args.output_stream.write(f">{name}\n{sequence}\n{ss} ({mfe:6.2f})\n")
        except:
            import subprocess
            # default to using the command line program RNAfold if Python API is not accessible
            run = ['RNAfold', '-C', '--noPS']
            if args.enforce:
                run.append('--enforceConstraint')

            data = f">{name}\n{sequence}\n{constraint}\n".encode('ascii')
            p = subprocess.Popen(run, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
            ret = p.communicate(input=data)
            args.output_stream.write(ret[0].decode())
    else:
        args.output_stream.write(f">{name}\n{sequence}\n{constraint}\n")


def _make_pair_table(structure):
    """
    Helper function to create a pair table from a dot-bracket structure.
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


def structure_from_dotplot(fname, n, threshold = 0.9):
    """
    Parse an RNAalifold dot-plot and retrieve all base pairs of the predicted MFE 
    structure that have probability above a certain threshold.
  
    Parameters
    ----------

    fname: str
        The dot-plot file name
    n: int
        The expected length of the output structure
    threshold: double
        A threshold between 0 and 1

    Returns
    -------
    str
        The dot-bracket string as obtained from the dot-plot or None on any error
    """
    import re
    structure = None

    prob_pat = re.compile(r"^\s*(\d+\.?\d*)\s+(\d+\.?\d*)\s+hsb\s+(\d+)\s+(\d+)\s+(\d+\.?\d*)\s+lbox")

    try:
        with open(fname) as f:
            if f.readline().startswith("%!PS"):
                structure = [ '.' for _ in range(n) ]

                for line in f.readlines():
                    line = line.strip()

                    # omit comment lines
                    if line.startswith("%"):
                        continue

                    m = prob_pat.match(line)
                    if m:
                        if m.lastindex == 5 and \
                           float(m.group(5)) > threshold and \
                           int(m.group(4)) <= n:
                            structure[int(m.group(3)) - 1] = "("
                            structure[int(m.group(4)) - 1] = ")"

                structure = "".join(structure)
    except:
        pass

    return structure


def map_structure_to_seqs(alignment, structure, cb = _hard_constraints_output, cb_data = None, min_hp_size = 3, canonical_pairs = _BPs):
    """
    Given a multiple sequence alignment (MSA) and a corresponding consensus structure,
    map the consensus structure to each sequence in the alignment.
    
    The mapping is such that for each individual sequence, retained base pairs must
    be canonical, i.e. Watson-Crick or Wobble pairs, and they must not involve gap
    positions. Otherwise, they are removed Finally, gaps are removed and the
    three results `name`, `sequence`, and `structure` are returned through a callback
    mechanism. Additionally, the callback receives as fourth argument the `cb_data`
    object for further processing if required.

    Parameters
    ----------
    
    alignment:  list(tuple())
        A list of (name, aligned_seq) pairs representing the MSA
    structure:  str
        The consensus structure for the MSA
    cb: func
        A callback function of the form cb(name, sequence, structure, cb_data)
    cb_data: object
        Any data passed through to the callback functionAn argument object specifying additional options
    min_hp_size: int
        The minimum number of unpaired nucleotides in a hairpin loop
    canonical_pairs:  iterable
        A list, tuple or dictionary containing upper-case entries "XY" where X and Y are nucleotides forming a canonical base pair (X,Y)
    """
    pt        = None
    has_vrna  = False

    try:
        import RNA
        pt        = RNA.ptable(structure)
        has_vrna  = True
    except:
        pt = _make_pair_table(structure)

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
                if seql[i].upper() + seql[pt[i]].upper() not in canonical_pairs:
                    cons[i]     = "."
                    cons[pt[i]] = "."

        seql = "".join(seql).replace("-", "")
        cons = "".join(cons).replace("x", "")

        # enforce minimum hairpin length constraint
        pts = None
        if has_vrna:
            pts = RNA.ptable(cons)
        else:
            pts = _make_pair_table(cons)

        for i in range(1, pts[0] + 1):
            if pts[i] > i and \
               (pts[i] - i - 1) < min_hp_size:
                pts[pts[i]] = 0
                pts[i] = 0

        # convert back to dot-bracket
        cons = [ "." for _ in range(pts[0]) ]
        for i in range(1, pts[0] + 1): 
            if pts[i] > i:
                cons[i - 1] = '('
                cons[pts[i] - 1] = ')'

        cons = "".join(cons)

        # execture callback
        cb(name, seql, cons, cb_data)


def hard_constraints(args):
    """
    Legacy refold.pl strategie using hard constraints
    """
    # try parsing the input alignment to work on
    aln = _read_alignment(args)

    # only continue if we successfully read an input alignment
    if aln != None:
        n = len(aln[0][1])
        # read consensus structure information, either from
        # RNAalifold output or dot-plot. 

        # Try dotplot first...
        fname = None

        if args.dotplot:
            fname = args.dotplot
        elif args.structfile:
            fname = args.structfile

        structure = structure_from_dotplot(fname, n, args.threshold)

        # Try RNAalifold ouput if dot-plot parsing was unsuccessful
        if structure == None:
            structure = _structure_from_alifold(args)

        if structure:
            map_structure_to_seqs(aln, structure, min_hp_size = args.turn, cb_data = args)

            return 0 # Success

    return 1  # Failure


def main():
    import sys
    import argparse

    def set_default_subparser(parser,
                              default_subparser):
        """Set a default subparser by subparser name.
        
        This function goes through the argument list sys.argv and checks
        whether a subparser is specified. If this is not the case, it
        checks whether any non-global arguments have been provided. Among
        this set, it tests if they belong to the default parser as specified.
        If this is the case or if the set is empty it adds the default
        parser as subparser, i.e. first argument to sys.argv. Otherwise
        it returns without any changes.
        
        Note
        ----
        
            Call after setup, just before parse_args()

        Parameters
        ----------
        
        parser: class(ArgumentParser)
            The name of the parser you're making changes to
        default_subparser: str
            The name of the subparser to call by default
        """
        subparser_found   = False
        subparser_options = False
        global_args       = parser._option_string_actions.keys()
        default_sp_args   = []
        non_global_args   = []

        # determine default subparser arguments
        for x in parser._subparsers._actions:
                if not isinstance(x, argparse._SubParsersAction):
                    continue
                if default_subparser in x._name_parser_map.keys():
                    sp = x._name_parser_map[default_subparser]
                    default_sp_args = sp._option_string_actions.keys()

        # create list of arguments we need for checking
        non_global_args   = [ arg for arg in sys.argv[1:] if arg not in global_args ]
        non_default_args  = [ arg for arg in non_global_args if arg not in default_sp_args ]

        # now check whether a left-over non_default_arg actually is
        # a subparser command
        for x in parser._subparsers._actions:
            if not isinstance(x, argparse._SubParsersAction):
                continue
            for sp_name in x._name_parser_map.keys():
                if sp_name in non_default_args:
                    subparser_found = True

        # If no subparser has been specified and we have left-over non_global_args
        # insert default in first position before all other arguments
        if not subparser_found and non_global_args:
            sys.argv.insert(1, default_subparser) 


    core = argparse.ArgumentParser(add_help=False)

    # global options for all modes
    core.add_argument("-a",
                      "--alignment",
                      metavar = "filename",
                      type = str,
                      help="A multiple sequence alignment file name")
    core.add_argument("-o",
                      "--output",
                      type = str,
                      help="A file or directory where to store the output",
                      default = None)
    core.add_argument("--turn",
                      type = int,
                      help = "Minimum hairpin length",
                      default = 3)

    main_parser = argparse.ArgumentParser(
        prog = "RNAconsensus",
        description = """
        A program to predict RNA secondary structures for single sequences based on the information
        gained from a multiple sequence alignment of homologous sequences.
        """,
        epilog = "If in doubt our program is right, nature is at fault. Comments should be sent to rna@tbi.univie.ac.at.",
        parents=[core])


    # create subparsers for the different operating modes
    sub_parsers = main_parser.add_subparsers(title = 'Prediction Stratgies',
                                             description = """
                                             This programm allows for different strategies for the way the consensus structure
                                             information is incorporated for the single sequence predictions
                                             """,
                                             help='Available strategies')

    # 1. Legacy refold.pl mode where a structure prediction from RNAalifold is taken as
    # input to create hard constraints for predictions of structures for the single
    # sequences within the alignment. For that purpose, the program prints each
    # sequence and the consensus structure mapped to the sequence such that the
    # output can be used as input for RNAfold.
    hard_parser = sub_parsers.add_parser('hardcons',
                                         description = """
                                         In this mode, the program reads a multiple sequence alignment and a consensus
                                         secondary structure, either in the form of the standard output of RNAalifold,
                                         or in the form of a dot plot as produced by RNAalifold -p. For each sequence
                                         in the alignment it writes the name, sequence, and constraint structure to
                                         stdout in a format suitable for piping into RNAfold -C.

                                         The constraint string is produced by removing from the consensus structure
                                         all gaps, as well as all pairs not compatible with the particular sequence.
                                         If the structure is read from a dot plot file, only pairs with probability
                                         larger than some threshold (default 0.9) are used.
                                         """,
                                         help='The legacy refold.pl mode',
                                         parents=[core])

    hard_parser.add_argument("-r",
                             "--run",
                             help="Run the predictions instead of producing output that can be used with RNAFold",
                             action="store_true")
    hard_parser.add_argument("--enforce",
                             help = "Enforce the structure constraint for single sequence predictions",
                             action="store_true")

    hard_g1 = hard_parser.add_argument_group("Consensus MFE",
                                             "Map the consensus MFE structure as predicted by RNAalifold to each sequence in the input alignment")
    hard_g1.add_argument("-f",
                         "--alifold-output",
                         type = str,
                         help = "The output of RNAalifold for the input alignment",
                         default = None)

    hard_g2 = hard_parser.add_argument_group("Consensus Probabilities",
                                             """
                                             Create structure constraints for each sequence in the input alignment from
                                             consensus base pair probabilities as predicted by RNAalifold""")
    hard_g2.add_argument("-d",
                         "--dotplot",
                         metavar = "filename",
                         type = str,
                         help = "The dot-plot created by RNAalifold for the input alignment",
                         default = None)
    hard_g2.add_argument("-t",
                         "--threshold",
                         type = float,
                         help = "Base pair probability threshold. Take only base pairs with probability larger than threshold into account",
                         default = 0.9)

    hard_parser.add_argument("alnfile",
                             nargs='?',
                             metavar = "alignment_file",
                             help = "The multiple sequence alignment file name if not specified by -a")
    hard_parser.add_argument("structfile",
                             nargs='?',
                             metavar = "structure_file",
                             help = "The RNAalifold output or dot plot file name")
    hard_parser.set_defaults(func = hard_constraints)

    soft_parser = sub_parsers.add_parser('softcons',
                                         description = """
                                         """,
                                         help='RNAsoftcons mode',
                                         parents=[core])


    set_default_subparser(main_parser, "hardcons")

    if len(sys.argv)==1:
        main_parser.print_help(sys.stderr)
        sys.exit(1)

    args = main_parser.parse_args()

    # set default output stream and closing method
    if args.output:
        args.output_stream = open(args.output, "w")
    else:
        args.output_stream = sys.stdout

    # call the function corresponding to the chosen mode
    # and use its return value as exit code
    ret = args.func(args)

    # close output stream if necessary
    if args.output_stream != sys.stdout:
        args.output_stream.close()

    sys.exit(ret)



if __name__ == '__main__':
    main()
