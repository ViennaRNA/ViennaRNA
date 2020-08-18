# Kinfold - Kinetic Folding Program for Nucleic Acids

  Copyright (C) 2001 Christoph Flamm, Ivo L. Hofacker

Kinfold uses routines from the ViennaRNA Package for energy evaluation.  You
can download the source code of the ViennaRNA Package from
https://www.tbi.univie.ac.at/RNA/

Comments are welcome.
 - Christoph Flamm, <xtof@tbi.univie.ac.at>

## Examples:

### Default mode: first passage time
The start structure is the open chain, the stop structure is the minimum free
energy (MFE) structure. Kinfold reports the first time point at wich the MFE
structure was observed (or the structure found at timepoint `--time`).  The
example output below is a possible trajectory for the sequence ACUGAUCGUAGUCAC.

```
Kinfold --time 100000 < seq.in
...............    0.00      2.202
.(.....).......    3.20      2.276
...............    0.00      3.681
.......(....)..    2.50      3.706
...(........)..    2.50      3.834
..((........)).    1.70      4.115
..(((......))).    0.60      4.332
...((......))..    1.40      4.399
...(((....)))..    0.10      4.625
..((((....)))).   -0.70      4.945 X1
```

The trajectory lists stucture, energy, and the time *until* the given
structure has been observed. (That is, the time at which the structure
shown in the current line transitions to the structure shown in the 
next line.) The X1 signifies that the trajectory terminated in the first 
specified stop structure. In addition the logfile kinout.log contains
information needed to reproduce the simulation results such as start and
stop structures, options, and random seeds for every simulation. Beware 
that Kinfold appends output, if the logfile already exists.

```
#<
#Date: Tue Jul 21 16:21:27 2020
#EnergyModel: dangle=2 Temp=37.0 logML=logarithmic Par=None
#MoveSet: noShift=off noLP=off
#Simulation: num=1 time=100000.00 seed=clock fpt=on rect=off mc=Kawasaki
#Simulation: phi=1 pbounds=0.1 0.1 2
#Output: log=kinout silent=off lmin=off cut=20.00
#ACUGAUCGUAGUCAC
#............... (  0.00)
#..((((....)))). ( -0.70) X01
(60400 33695 55886) X01        4.945

```

All times are given in internal units that can be translated into real time
only by comparison with experiment. Very roughly one time step corresponds to
about 10 micro seconds.

### Simulate transcription process
To run a folding during transcription simulation use the `--grow` option.
Assuming a transcript length of 50 nt, a transcription rate of 100 nt/sec and 1
sec about 10^7 time steps we could use:

```
$ Kinfold --grow 1e5 --glen 10 --time 4e6 < seq.in
```

