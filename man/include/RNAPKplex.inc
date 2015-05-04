/The algorithm furthermore always considers dangle=2 model./

It  also  produces a PostScript file with a plot of the pseudoknot-free
secondary structure graph, in which the bases  forming  the  pseuodknot
are marked red.

Sequences are read in a simple text format where each sequence occupies
a single line. Each sequence may be preceded by a line of the form
.nf
.ft CW
> name
.ft
.fi
to assign a name to the sequence. If a name is given in the input, the
 PostScript file "name.ps" is produced for the structure graph.  Other-
wise  the  file  name defaults to PKplex.ps. Existing files of the same
name will be overwritten.
The input format is similar to fasta except that  even  long  sequences
may  not  be  interrupted  by  line  breaks,  and  the header lines are
optional.  The program will continue to read new sequences until a line
consisting  of  the  single  character @ or an end of file condition is
encountered.

[AUTHOR]

Wolfgang Beyer
