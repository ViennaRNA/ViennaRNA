========================
The Program ``RNApvmin``
========================

.. contents:: Table of Contents
    :maxdepth: 1
    :local:


Introduction
============

The program ``RNApvmin`` reads a RNA sequence from *stdin* and uses an iterative minimization
process to calculate a perturbation vector that minimizes the discripancies between predicted
pairing probabilites and observed pairing probabilities (deduced from given shape
reactivities).
The experimental SHAPE data has to be present in the file format described above.
The application will write the calculated vector of perturbation energies to *stdout*,
while the progress of the minimization process is written to *stderr*.
The resulting perturbation vector can be interpreted directly and gives usefull insights into the
discrepancies between thermodynamic prediction and experimentally determined pairing status.
In addition the perturbation energies can be used to constrain folding with ``RNAfold``:

.. code::

  $ RNApvmin rna.shape < rna.seq >vector.csv
  $ RNAfold --shape=vector.csv --shapeMethod=W < rna.seq

The perturbation vector file uses the same file format as the SHAPE data file.
Instead of SHAPE reactivities the raw perturbation energies will be storred in the last column.
Since the energy model is only adjusted when necessary, the calculated perturbation energies may be used
for the interpretation of the secondary structure prediction, since they indicate
which positions require major energy model adjustments in order to yield a prediction
result close to the experimental data. High perturbation energies for just
a few nucleotides may indicate the occurrence of features, which are not explicitly
handled by the energy model, such as posttranscriptional modifications and
intermolecular interactions.
