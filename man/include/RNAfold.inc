[AUTHOR]

Ivo L Hofacker, Walter Fontana, Sebastian Bonhoeffer, Peter F Stadler, Ronny Lorenz

[EXAMPLES]
Single line sequence input and calculation of partition function and MEA structure

.nf
.ft CW
  $ RNAfold --MEA -d2 -p
.ft
.fi
  
The program will then prompt for sequence input. Using the example sequence
"CGACGTAGATGCTAGCTGACTCGATGC" and pressing ENTER the output of the program will be
similar to

.nf
.ft CW
  CGACGUAGAUGCUAGCUGACUCGAUGC
  (((.((((.......)).)))))....
   minimum free energy =  -1.90 kcal/mol
  (((.((((.......))},})))....
   free energy of ensemble =  -2.86 kcal/mol
  (((.(.((.......))..)))).... {  0.80 d=2.81}
  (((.((((.......))).)))).... { -1.90 MEA=22.32}
   frequency of mfe structure in ensemble 0.20997; ensemble diversity 4.19
.ft
.fi


Here, the first line just repeats the sequence input. The second line contains a
MFE structure in dot bracket notation followed by the minimum free energy. After
this, the pairing probabilities for each nucleotide are shown in a pseudo dot-bracket
notation followed by the free energy of ensemble. The next two lines show the centroid
structure with its free energy and its distance to the ensemble as well as the MEA structure,
its free energy and the maximum expected accuracy, respectively. The last line finally
contains the frequency of the MFE representative in the complete ensemble of secondary
structures and the ensemble diversity. For further details about the calculation and
interpretation of the given output refer to the reference manual of RNAlib.

Since version 2.0 it is also possible to provide FASTA file sequence input. Assume
you have a file containing two sequences in FASTA format, e.g

.nf
.ft CW
  $ cat sequences.fa
  >seq1
  CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUG
  GAACGAUCUAUAACACGACUUCACUCUU
  >seq2
  GAAUGACCCGAUAACCCCGUAAUAUUUGGAACGAUCUA
  UAACACGACUUCACUCUU
.ft
.fi

In order to compute the MFE for the two sequences the user can use the following
command

.nf
.ft CW
  $ RNAfold < sequences.fa
.ft
.fi

which would result in an output like this

.nf
.ft CW
  >seq1
  CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUGGAACGAUCUAUAACACGACUUCACUCUU
  .((.(((...((((..(((((........)))))))))...))).))................... ( -5.40)
  >seq2
  GAAUGACCCGAUAACCCCGUAAUAUUUGGAACGAUCUAUAACACGACUUCACUCUU
  .......((((..............))))........................... ( -2.00)
.ft
.fi

[CONSTRAINT EXAMPLES]

Secondary structure constraints may be given in addition to the sequence information, too.
Using the first sequence of the previous example and restricting the nucleotides of the
outermost helix to be unpaired, i.e. base pairs (2,47) and (3,46) the input file should
have the following form

.nf
.ft CW
  $ cat sequence_unpaired.fa
  >seq1
  CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUG
  GAACGAUCUAUAACACGACUUCACUCUU
  .xx...................................
  .......xx...................
.ft
.fi

Calling RNAfold with the structure constraint option -C it shows the following result

.nf
.ft CW
  $ RNAfold -C < sequence_unpaired.fa
  >seq1
  CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUGGAACGAUCUAUAACACGACUUCACUCUU
  ....(((...((((..(((((........)))))))))...)))...................... ( -4.20)
.ft
.fi

This represents the minimum free energy and a structure representative of the RNA
sequence given that nucleotides 2,3,46 and 47 must not be involved in any base pair.
For further information about constrained folding refer to the details of the -C option
and the reference manual of RNAlib.

Since version 2.2 the ViennaRNA Package distinguishes hard and soft constraints.
As a consequence, structure predictions are easily amenable to a versatile set of constraints,
such as maximal base pair span, incorporation of SHAPE reactivity data, and RNA-ligand binding
to hairpin, or interior loop motifs.

Limit the maximum distance of basepairs, for instance to 600nt:

.nf
.ft CW
  $ RNAfold --maxBPspan 600
.ft
.fi

Use SHAPE reactivity data to guide secondary structure prediction:

.nf
.ft CW
  $ RNAfold --shape=reactivities.dat < sequence.fa
.ft
.fi

where the file reactivities.dat is a two column text file with sequence positions (1-based)
and normalized reactivity values. Missing values may be left out, or assigned a negative
score:

.nf
.ft CW
  $ cat reactivities.dat
  9    -999       # No reactivity information
  10   -999
  11   0.042816   # normalized SHAPE reactivity
  12   0          # also a valid SHAPE reactivity
  15   0.15027    # Missing data for pos. 13-14
  ...
  42   0.16201
.ft
.fi

Note, that RNAfold will only process the first sequence in the input file, when provided
with SHAPE reactivity data!

Structure constraints beyond those that can be expressed with a pseudo-dot bracket notation
may be provided in a constraints definition file:

.nf
.ft CW
  $ RNAfold --constraint=constraints.txt < sequence.fa
.ft
.fi

where the constraints definition file syntax is a generalization of those used in
UNAfold/mfold. Each line starts with a single character command, followed by three
or four numbers. In addition, optional auxiliary modifier characters may be used
to limit the constraint to specific loop types. For base pair specific constraints,
we currently distinguish pairs in exterior loops (E), closing pairs of hairpin
loops (H), closing (I) and enclosed (i) pairs of interior loops, and closing (M)
and enclosed (m) pairs of multibranch loops. Nucleotide-wise constraints may be
limited to their loop context using the corresponding uppercase characters. The
default is to apply a constraint to all (A) loop types. Furthermore, pairing
constraints for single nucleotides may be limited to upstream (U), or downstream (D)
orientation. The constraint definition file specification is as follows:

.nf
.ft CW
  F i 0 k   [TYPE] [ORIENTATION] # Force nucleotides i...i+k-1 to be paired
  F i j k   [TYPE] # Force helix of size k starting with (i,j) to be formed
  P i 0 k   [TYPE] # Prohibit nucleotides i...i+k-1 to be paired
  P i j k   [TYPE] # Prohibit pairs (i,j),...,(i+k-1,j-k+1)
  P i-j k-l [TYPE] # Prohibit pairing between two ranges
  C i 0 k   [TYPE] # Nucleotides i,...,i+k-1 must appear in context TYPE
  C i j k          # Remove pairs conflicting with (i,j),...,(i+k-1,j-k+1)
  E i 0 k e        # Add pseudo-energy e to nucleotides i...i+k-1
  E i j k e        # Add pseudo-energy e to pairs (i,j),...,(i+k-1,j-k+1)

                   # [TYPE]        = { E, H, I, i, M, m, A }
                   # [ORIENTATION] = { U, D }
.ft
.fi

Again, RNAfold by default only processes the first sequence in the input sequence
when provided with a constraints definition file. To apply the exact same constraints
to each of the input sequences in a multi FASTA file, use the batch mode commandline
option:

.nf
.ft CW
  $ RNAfold --constraint=constraints.txt --batch < sequences.fa
.ft
.fi


[REFERENCES]

M. Zuker, P. Stiegler (1981),
"Optimal computer folding of large RNA sequences using thermodynamic and auxiliary information",
Nucl Acid Res: 9, pp 133-148

J.S. McCaskill (1990),
"The equilibrium partition function and base pair binding probabilities for RNA secondary structures",
Biopolymers: 29, pp 1105-1119

I.L. Hofacker & P.F. Stadler (2006),
"Memory Efficient Folding Algorithms for Circular RNA Secondary Structures",
Bioinformatics

A.F. Bompfuenewerer, R. Backofen, S.H. Bernhart, J. Hertel, I.L. Hofacker, P.F. Stadler, S. Will (2007),
"Variations on {RNA} Folding and Alignment: Lessons from Benasque",
J. Math. Biol.

D. Adams (1979),
"The hitchhiker's guide to the galaxy",
Pan Books, London

The calculation of mfe structures is based on dynamic programming algorithm originally developed by M. Zuker and P. Stiegler. The partition function algorithm is based on work by J.S. McCaskill.
