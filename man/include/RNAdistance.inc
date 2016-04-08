/is not recommended./

RNAdistance accepts structures in bracket format, where matching brackets
symbolize base pairs and unpaired bases are represented by a dot ".",
or coarse grained representations where hairpins, interior loops,
bulges, multiloops, stacks and external bases are represented by
(H), (I), (B), (M), (S), and (E), respectively. These can be optionally
weighted. Full structures can be represented in the same fashion using
the identifiers (U) and (P) for unpaired and paired bases, respectively.
We call this the HIT representation (you don't want to know what this means).
For example the following structure consists of 2 hairpins joined by
a multiloop:

.nf
.ft CW
  .((..(((...)))..((..)))).       full structure (usual format);
  (U)((U2)((U3)P3)(U2)((U2)P2)P2) HIT structure;
  ((H)(H)M)  or
  ((((H)S)((H)S)M)S)              coarse grained structure;
  (((((H3)S3)((H2)S2)M4)S2)E2)    weighted coarse grained.
.ft
.fi

The program will continue to read new structures until a line consisting
of the single character @ or an end of file condition is encountered. Input
lines neither containing a valid structure nor starting with > are ignored.


[AUTHOR]

Walter Fontana, Ivo L Hofacker, Peter F Stadler

[>REFERENCES]

B.A. Shapiro (1988),
"An algorithm for comparing multiple RNA secondary structures"
CABIOS: 4, pp 381-393

B.A. Shapiro, K. Zhang (1990),
"Comparing multiple RNA secondary structures using tree comparison",
CABIOS: 6, pp 309-318

W. Fontana, D.A.M. Konings, P.F. Stadler and P. Schuster P (1993),
"Statistics of RNA secondary structures",
Biopolymers: 33, pp 1389-1404
