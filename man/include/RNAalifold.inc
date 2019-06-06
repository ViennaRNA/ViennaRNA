[AUTHOR]

Ivo L Hofacker, Stephan Bernhart, Ronny Lorenz

[EXAMPLES]

A simple call to compute consensus MFE structure, ensemble free energy,
base pair probabilities, centroid structure, and MEA structure for a
multiple sequence alignment (MSA) provided as Stockholm formatted file
alignment.stk might look like:

.nf
.ft CW
  $ RNAalifold -p --MEA alignment.stk
.ft
.fi

Consider the following MSA file for three sequences

.nf
.ft CW
  # STOCKHOLM 1.0

  #=GF AC   RF01293
  #=GF ID   ACA59
  #=GF DE   Small nucleolar RNA ACA59
  #=GF AU   Wilkinson A
  #=GF SE   Predicted; WAR; Wilkinson A
  #=GF SS   Predicted; WAR; Wilkinson A
  #=GF GA   43.00
  #=GF TC   44.90
  #=GF NC   40.30
  #=GF TP   Gene; snRNA; snoRNA; HACA-box;
  #=GF BM   cmbuild -F CM SEED
  #=GF CB   cmcalibrate --mpi CM
  #=GF SM   cmsearch --cpu 4 --verbose --nohmmonly -E 1000 -Z 549862.597050 CM SEQDB
  #=GF DR   snoRNABase; ACA59;
  #=GF DR   SO; 0001263; ncRNA_gene;
  #=GF DR   GO; 0006396; RNA processing;
  #=GF DR   GO; 0005730; nucleolus;
  #=GF RN   [1]
  #=GF RM   15199136
  #=GF RT   Human box H/ACA pseudouridylation guide RNA machinery.
  #=GF RA   Kiss AM, Jady BE, Bertrand E, Kiss T
  #=GF RL   Mol Cell Biol. 2004;24:5797-5807.
  #=GF WK   Small_nucleolar_RNA
  #=GF SQ   3


  AL031296.1/85969-86120     CUGCCUCACAACGUUUGUGCCUCAGUUACCCGUAGAUGUAGUGAGGGUAACAAUACUUACUCUCGUUGGUGAUAAGGAACAGCU
  AANU01225121.1/438-603     CUGCCUCACAACAUUUGUGCCUCAGUUACUCAUAGAUGUAGUGAGGGUGACAAUACUUACUCUCGUUGGUGAUAAGGAACAGCU
  AAWR02037329.1/29294-29150 ---CUCGACACCACU---GCCUCGGUUACCCAUCGGUGCAGUGCGGGUAGUAGUACCAAUGCUAAUUAGUUGUGAGGACCAACU
  #=GC SS_cons               -----((((,<<<<<<<<<___________>>>>>>>>>,,,,<<<<<<<______>>>>>>>,,,,,))))::::::::::::
  #=GC RF                    CUGCcccaCAaCacuuguGCCUCaGUUACcCauagguGuAGUGaGgGuggcAaUACccaCcCucgUUgGuggUaAGGAaCAgCU
  //
.ft
.fi


Then, the above program call will produce this output:

.nf
.ft CW
  3 sequences; length of alignment 84.
  >ACA59
  CUGCCUCACAACAUUUGUGCCUCAGUUACCCAUAGAUGUAGUGAGGGUAACAAUACUUACUCUCGUUGGUGAUAAGGAACAGCU
  ...((((((.(((((((((...........))))))))).))))))..........(((((......)))))............ (-12.54 = -12.77 +   0.23)
  ...((((((.(((((((((...........))))))))).)))))){{,.......{{{{,......}))))............ [-14.38]
  ...((((((.(((((((((...........))))))))).))))))..........((((........))))............ {-12.44 = -12.33 +  -0.10 d=10.94}
  ...((((((.(((((((((...........))))))))).))))))..........((((........))))............ {-12.44 = -12.33 +  -0.10 MEA=66.65}
   frequency of mfe structure in ensemble 0.368739; ensemble diversity 17.77 
.ft
.fi

Here, the first line is written to \fIstderr\fR and simply states the number of sequences and
the length of the alignment. This line can be suppressed using the \fB\-\-quiet\fR option.
The main output then consists of 7 lines, where the first two resemble the FASTA header
with the ID as read from the input data set, followed by the consensus sequence in the
second line. The third line consists of the consensus secondary structure in dot-bracket
notation followed by the averaged minimum free energy in parenthesis. This energy is
composed of two major contributions, the actual free energies derived from the Nearest
Neighbor model, and the covariance pseudo-energy term, which are both displayed after
the equal sign. The fourth line shows the base pair propensity in pseudo dot-bracket
notation followed by the ensemble free energy dG = -kT ln(Z) in square brackets.
Similarly, the next two lines state the controid- and the MEA structure in dot-bracket
notation, followed by their corresponding free energy contributions, the mean distance
(d) to the ensemble as well as the maximum expected accuracy (MEA). Again, the free
energies are split into Nearest Neighbor contribution and the covariance pseudo-energy
term.

Furthermore, RNAalifold will produce three output files: ACA59_ss.ps, ACA59_dp.ps, and
ACA59_ali.out that contain the secondary structure drawing, the base pair probability
dot-plot, and a detailed table of base pair probabilities, respectively.


.SH "THE ALIOUT FILE"

When computing base pair probabilities (\fB\-\-partfunc\fR option), RNAalifold will produce
a file with the suffix `ali.out`. This file contains the base pairing probabilities between
different alignment columns together with some detailed statistics for the individual
sequences within the alignment. The file is a simple text file with a two line header that
states the number of sequences and length of the alignment. The first couple of lines
of this file may look like:

.nf
.ft CW
  3 sequence; length of alignment 84
  alifold output
     14    36  0  92.7%   0.212 CG:1    UA:2   
     13    37  0  92.7%   0.218 GU:1    AU:2   
     12    38  0  92.7%   0.231 CG:3   
     15    35  0  91.9%   0.239 UG:3   
     16    34  0  85.2%   0.434 UA:2    --:1   
      8    42  0  80.7%   0.526 AU:3   +
      9    41  0  80.4%   0.542 CG:3   +
      7    43  1  80.1%   0.541 CG:2   +
.ft
.fi

Starting with the third row, there are at least six and at most 13 columns separated by
whitespaces stating: (1) the i-position and (2) the j-position of a potential base pair
(i, j), followed by (3) the number of counter examples, i.e. the number of sequences in
the alignment that can't form a canonical base pair with their respective sequence positions.
Next is (4) the base pair probabilitiy in percent, (5) a pseudo entropy measure
S_ij = S_i + S_j - p_ij ln(p_ij), where S_i and S_j are the positional entropies for the
two alignment columns i and j, and p_ij is the base pair probability. Finally, the last
columns (6-12) state the number of particular base pairs for the individual sequences in
the alignment. Here, we distinguish the base pairs "GC","CG","AU","UA","GU","UG", and
the special case "\-\-" that represents gaps at both positions i and j.
Finally, base pairs that are not part of the MFE structure are marked by an additional
"+" sign in the last column.


[>REFERENCES]

The algorithm is a variant of the dynamic programming algorithms of M. Zuker and P. Stiegler (mfe)
and J.S. McCaskill (pf) adapted for sets of aligned sequences with covariance information.

Ivo L. Hofacker, Martin Fekete, and Peter F. Stadler (2002),
"Secondary Structure Prediction for Aligned RNA Sequences",
J.Mol.Biol.: 319, pp 1059-1066.

Stephan H. Bernhart, Ivo L. Hofacker, Sebastian Will, Andreas R. Gruber, and Peter F. Stadler (2008),
"RNAalifold: Improved consensus structure prediction for RNA alignments",
BMC Bioinformatics: 9, pp 474

[SEE ALSO]

The ALIDOT package http://www.tbi.univie.ac.at/RNA/Alidot/
