[AUTHOR]

Dominik Luntzer, Ronny Lorenz

[EXAMPLES]

RNApvmin acceptes a SHAPE file and a corresponding nucleotide sequence, which is read form stdin.

.nf
.ft CW
  RNApvmin sequence.shape < sequence.fasta > sequence.pv
.ft
.fi

The normalized SHAPE reactivity data has to be stored in a text file, where each line contains the position
and the reactivity for a certain nucleotide ([position] [nucleotide] [SHAPE reactivity]).

.nf
.ft CW
  1 A 1.286
  2 U 0.383
  3 C 0.033
  4 C 0.017
  ...
  ...
  98 U 0.234
  99 G 0.885
.ft
.fi

The nucleotide information in the SHAPE file is optional and will be used to cross check the given input sequence if present.
If SHAPE reactivities could not be determined for every nucleotide, missing values can simply be omited.

The progress of the minimization will be printed to stderr. Once a solution was found, the calculated perturbation vector
will be print to stdout and can then further be used to constrain RNAfold's MFE/partition function calculation by applying
the perturbation energies as soft constraints.

.nf
.ft CW
  RNAfold --shape=sequence.pv --shapeMethod=W < sequence.fasta
.ft
.fi

[>REFERENCES]

S. Washietl, I.L. Hofacker, P.F. Stadler, M. Kellis (2012)
"RNA folding with soft constraints: reconciliation of probing data and thermodynamics secondary structure prediction"
Nucl Acids Res: 40(10), pp 4261-4272

