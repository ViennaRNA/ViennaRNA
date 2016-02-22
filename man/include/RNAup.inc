/over all possible conformations./

RNAup provides two different modes: By default RNAup computes accessibilities,
in terms of the free energies needed to open a region (default length 4). It
prints the region of highest accessibility and its opening energy to stdout,
opening energies for all other regions are written to a file.

.br
In interaction mode the interaction between two RNAs is calculated. It is
invoked if the input consists of two sequences concatenated with an "&",
or if the options -X[pf] or -b are given. Unless the -b option is specified
RNAup assumes that the longer RNA is a structured target sequence
while the shorter one is an unstructured small RNA.
.br
Additionally, for every position along the target sequence we write the best
free energy of binding for an interaction that includes this position to the
the output file.
Output to stdout consists of the location and free energy, dG,
for the optimal region of interaction. The binding energy dG is also split into
its components the interaction energy dGint and the opening energy dGu_l (and
possibly dGu_s for the shorter sequence).
.br
In addition we print the optimal interaction structure as computed by RNAduplex
for this region. Note that it can happen that the RNAduplex computed optimal
interaction does not coincide with the optimal RNAup region. If the two
predictions don't match the structure string is replaced by a run of "."
and a message is written to stderr.
.br

Each sequence should be in 5' to 3' direction. If the sequence is preceded
by a line of the form
.nf
.ft CW
> name
.ft
.fi

the output file "name_ux_up.out" is produced, where the "x" in "ux" is the
value set by the -u option. Otherwise the file name defaults to
RNA_ux_up.out. The output is concatenated if a file with the same name exists.
.br

RNA sequences are read from stdin as strings of characters. White space and
newline within a sequence cause an error! Newline is used to separate
sequences. The program will continue to read new sequences until a line
consisting of the single character @ or an end of file condition is
encountered.

[AUTHOR]

Ivo L Hofacker, Peter F Stadler, Ulrike Mueckstein, Ronny Lorenz

[EXAMPLES]

.B Output to stdout:

In Interaction mode RNAup prints the most favorable interaction energy
between the two sequences to stdout. The most favorable interaction energy
(dG) depends on the position in the longer sequence (region [i,j]) and the
position in the shorter sequence (region[k,l]): dG[i,j;k,l].  dG[i,j;k,l] is the
largest contribution to dG[i,j] = sum_kl dG[i,j;k,l] which is given in the
output file: therefore dG[i,j;k,l] <= dG[i,j].

.nf
.ft CW
  '....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8'
  > franz
  GGAGUAGGUUAUCCUCUGUU
  > sissi
  AGGACAACCU
  dG = dGint + dGu_l
  (((((.((((&)))).)))))   6,15  :   1,10  (-6.66 = -9.89 + 3.23)
  AGGUUAUCCU&AGGACAACCU
  RNAup output in file: franz_sissi_w25_u3_4_up.out
.ft
.fi

where the result line contains following information

.nf
.ft CW
  RNAduplex results       [i,j]     [k,l]    dG = dGint + dGu_l
  (((((.((((&)))).)))))   6,15   :  1,10     (-6.66=-9.89+3.23)
.ft
.fi


.RD
.B Output to file:

Output to file contains a header including date, the command line of the
call to RNAup, length and names of the input sequence(s) followed
by the sequence(s). The first sequence is the target sequence.
Printing of the header can be turned off using the -nh option.

The line directly after the header gives the column names for the output:

.nf
.ft CW
  position     dGu_l for -u 3      dGu_l for -u 4       dG
#     pos      u3S       u3H       u4S       u4H        dG
.ft
.fi

where all information refers to the target sequence. The dGu_l column contains
information about the -u value (u=3 or u=4) and the contribution to the free
energy to open all structures "S" or only hairpin loops "H", see option -c.
NA means that no results is possible (e.g. column u3S row 2: no region of
length 3 ending at position 2 exists).

.nf
.ft CW
#  Thu Apr 10 09:15:11 2008
#  RNAup -u 3,4 -c SH -b
#  20 franz
#  GGAGUAGGUUAUCCUCUGUU
#  10 sissi
#  AGGACAACCU
#     pos      u3S       u3H       u4S       u4H        dG
       1        NA        NA        NA        NA    -1.540
       2        NA        NA        NA        NA    -1.540
       3     1.371        NA        NA        NA    -1.217
       4     1.754     5.777     1.761        NA    -1.393
       5     1.664     3.140     1.811     5.800    -1.393
.ft
.fi


If the -b option is selected position and dGu_s values for the shorter sequence
are written after the information for the target sequence.

[>REFERENCES]

U. Mueckstein, H. Tafer, J. Hackermueller, S.H. Bernhart, P.F. Stadler, and I.L. Hofacker (2006),
"Thermodynamics of RNA-RNA Binding",
Bioinformatics: 22(10), pp 1177-1182
