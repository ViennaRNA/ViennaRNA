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

[>EXAMPLES]
.SH "CONSTRAINT EXAMPLES"

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

.I Restricting the maximal span of a base pair

A convenience commandline option allows you to easily limit the distance (j - i + 1) between
two nucleotides i and j that form a basepair. For instance a limit of 600nt can be accomplished
using:

.nf
.ft CW
  $ RNAfold --maxBPspan 600
.ft
.fi

.I Guide structure prediction with SHAPE reactivity data

Use SHAPE reactivity data to guide secondary structure prediction:

.nf
.ft CW
  $ RNAfold --shape=reactivities.dat < sequence.fa
.ft
.fi

where the file reactivities.dat is a two column text file with sequence positions (1-based)
and normalized reactivity values (usually between 0 and 2. Missing values may be left out,
or assigned a negative score:

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

.I Complex structure constraints and grammar extensions

Structure constraints beyond those that can be expressed with a pseudo-dot bracket notation
may be provided in a so-called command file:

.nf
.ft CW
  $ RNAfold --commands=constraints.txt < sequence.fa
.ft
.fi

The command file syntax is a generalization of constraints as used in
UNAfold/mfold. Each line starts with a one or two letter command followed
by command parameters. For structure constraints, this amounts to a single
command character followed by three or four numbers. In addition, optional
auxiliary modifier characters may be used to limit the constraint to specific
loop types. For base pair specific constraints,
we currently distinguish pairs in exterior loops (E), closing pairs of hairpin
loops (H), closing (I) and enclosed (i) pairs of interior loops, and closing (M)
and enclosed (m) pairs of multibranch loops. Nucleotide-wise constraints may be
limited to their loop context using the corresponding uppercase characters. The
default is to apply a constraint to all (A) loop types. Furthermore, pairing
constraints for single nucleotides may be limited to upstream (U), or downstream (D)
orientation. The command file specification is as follows:

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
  UD m e    [LOOP] # Add ligand binding to unstructured domains with motif
                   # m and binding free energy e

                   # [LOOP]        = { E, H, I, M, A }
                   # [TYPE]        = [LOOP] + { i, m }
                   # [ORIENTATION] = { U, D }
.ft
.fi

Again, RNAfold by default only processes the first sequence in the input sequence
when provided with constraints in a command file. To apply the exact same constraints
to each of the input sequences in a multi FASTA file, use the batch mode commandline
option:

.nf
.ft CW
  $ RNAfold --constraint=constraints.txt --batch < sequences.fa
.ft
.fi

.I Ligand binding contributions to specific hairpin/interior loop motifs

A convenience function allows one to specify a hairping/interior loop motif where a ligand
is binding with a particular binding free energy dG.
Here is an example that adds a theophylline binding motif. Free energy contribution of
this motif of dG=-9.22kcal/mol is derived from k_d=0.32umol/l, taken from Jenison et al.
1994. Although the structure motif consists of a symmetric interior loop of size 6,
followed by a small helix of 3 basepairs, and a bulge of 3 nucleotides, the entire
structure can still be represented by one interior loop.
See the below mofif description where the '&' character splits the motif into a 5' and
a 3' part. The first line gives the sequences motif, the second line shows the actual
structure motif of the aptamer pocket, and the third line is the interior loop motif
that fully encapsulates the theophylline aptamer:

.nf
.ft CW
  GAUACCAG&CCCUUGGCAGC
  (...((((&)...)))...)
  (......(&).........)
.ft
.fi

To use the above information in the folding recursions of RNAfold, one only needs to
provide the motif itself, and binding free energy:

.nf
.ft CW
  $ RNAfold --motif="GAUACCAG&CCCUUGGCAGC,(...((((&)...)))...),-9.22" < sequences.fa
.ft
.fi

.I Ligand binding contributions to unpaired segments of the RNA structure

The extension of the RNA folding grammar with unstructured domains allows for an easy
incorporation of ligands that bind to unpaired stretches of an RNA structure. To
model such interactions only two parameters are required: (i) a sequence motif in
IUPAC notation that specifies where the ligand binds to, and (ii) a binding free
energy that can be derived from the association/dissociation constant of the ligand.
With these two parameters in hand, the modification of RNAfold to include the competition
of regular intramolecular base pairing and ligand interaction is as easy as writing
a simple command file of the form:

.nf
.ft CW
  UD m e    [LOOP]
.ft
.fi

where m is the motif string in upper-case IUPAC notation, and e the binding free energy
in kcal/mol and optional loop type restriction [LOOP]. See also the command file specification as defined above.

For instance, having a protein with a 4-nucleotide footprint binding 'AAAA', a
binding free energy e = -5.0 kcal/mol, and a binding restriction to exterior- and
multibranch loops results in a command file:

.nf
.ft CW
  $ cat commands.txt
  UD AAAA -5.0  ME
.ft
.fi

and the corresponding call to RNAfold to compute MFE and equilibrium probabilities becomes:

.nf
.ft CW
  $ RNAfold --commands=commands.txt -p < sequence.fa
.ft
.fi

The resulting MFE plot will be annotated to display the binding site(s) of the ligand,
and the base pair probability dot-plot is extended to include the probability that
a particular nucleotide is bound by the ligand.


[>REFERENCES]

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
