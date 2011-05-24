[AUTHOR]

Ivo L Hofacker, Walter Fontana, Sebastian Bonhoeffer, Peter F Stadler, Ronny Lorenz

[EXAMPLES]
Single line sequence input and calculation of partition function and MEA structure

  $ RNAfold --MEA -d2 -p
  
The program will then prompt for sequence input. Using the example sequence
"CGACGTAGATGCTAGCTGACTCGATGC" and pressing ENTER the output of the program will be
similar to

  CGACGUAGAUGCUAGCUGACUCGAUGC
  (((.((((.......)).)))))....
   minimum free energy =  -1.90 kcal/mol
  (((.((((.......))},})))....
   free energy of ensemble =  -2.86 kcal/mol
  (((.(.((.......))..)))).... {  0.80 d=2.81}
  (((.((((.......))).)))).... { -1.90 MEA=22.32}
   frequency of mfe structure in ensemble 0.20997; ensemble diversity 4.19

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

  $ cat sequences.fa
  >seq1
  CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUG
  GAACGAUCUAUAACACGACUUCACUCUU
  >seq2
  GAAUGACCCGAUAACCCCGUAAUAUUUGGAACGAUCUA
  UAACACGACUUCACUCUU

In order to compute the MFE for the two sequences the user can use the following
command

  $ RNAfold < sequences.fa

which would result in an output like this

  >seq1
  CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUGGAACGAUCUAUAACACGACUUCACUCUU
  .((.(((...((((..(((((........)))))))))...))).))................... ( -5.40)
  >seq2
  GAAUGACCCGAUAACCCCGUAAUAUUUGGAACGAUCUAUAACACGACUUCACUCUU
  .......((((..............))))........................... ( -2.00)

Secondary structure constraints may be given in addition to the sequence information, too.
Using the first sequence of the previous example and restricting the nucleotides of the
outermost helix to be unpaired, i.e. base pairs (2,47) and (3,46) the input file should
have the following form

  $ cat sequence_unpaired.fa
  >seq1
  CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUG
  GAACGAUCUAUAACACGACUUCACUCUU
  .xx...................................
  .......xx...................

Calling RNAfold with the structure constraint option -C it shows the following result

  $ RNAfold -C < sequence_unpaired.fa
  >seq1
  CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUGGAACGAUCUAUAACACGACUUCACUCUU
  ....(((...((((..(((((........)))))))))...)))...................... ( -4.20)

This represents the minimum free energy and a structure representative of the RNA
sequence given that nucleotides 2,3,46 and 47 must not be involved in any base pair.
For further information about constrained folding refer to the details of the -C option
and the reference manual of RNAlib.

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
