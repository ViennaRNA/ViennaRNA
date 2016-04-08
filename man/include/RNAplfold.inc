/for short stable RNA structures./

Output consists of a dot plot in postscript file, where the averaged pair probabilities
can easily be parsed and visually inspected.

The -u option makes i possible to compute the probability that a stretch of x consequtive
nucleotides is unpaired, which is useful for predicting possible binding sites. Again
this probability is averaged over all windows containing the region.

.B WARNING! Output format changed!!

The output is a plain text matrix containing on each line a position i followed by the
probability that i is unpaired, [i-1..i] is unpaired [i-2..i] is unpaired and so on to
the probability that [i-x+1..i] is unpaired.

[AUTHOR]

Stephan H Bernhart, Ivo L Hofacker, Peter F Stadler, Ronny Lorenz

[>REFERENCES]

S. H. Bernhart, I.L. Hofacker, and P.F. Stadler (2006),
"Local Base Pairing Probabilities in Large RNAs",
Bioinformatics: 22, pp 614-615

A.F. Bompfuenewerer, R. Backofen, S.H. Bernhart, J. Hertel, I.L. Hofacker, P.F. Stadler, S. Will (2007),
"Variations on RNA Folding and Alignment: Lessons from Benasque",
J. Math. Biol.

[SEE ALSO]

RNALfold(1)
