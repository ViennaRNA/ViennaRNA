=========================
The Program ``RNAsubopt``
=========================

.. contents:: Table of Contents
    :maxdepth: 1
    :local:


Introduction
============

``RNAsubopt`` calculates all suboptimal secondary structures within a
given energy range above the ``MFE`` structure. Be careful, the number
of structures returned grows exponentially with both sequence length and
energy range.

Suboptimal folding
==================

- Generate all suboptimal structures within a certain energy range
  from the ``MFE`` specified by the ``-e`` option::

    $ RNAsubopt -e 1 -s < test.seq 
    CUACGGCGCGGCGCCCUUGGCGA   -500    100
    ...........((((...)))).  -5.00
    ....((((...))))........  -4.80
    (((.((((...))))..)))...  -4.20
    ...((.((.((...)).)).)).  -4.10


The text output shows an energy sorted list (option ``-s``) of all secondary
structures within 1~kcal/mol of the ``MFE`` structure. Our sequence actually
has a ground state structure (-5.70) and three structures within 1~kcal/mol
range. 

``MFE`` folding alone gives no indication that there are actually a number
of plausible structures. Remember that ``RNAsubopt`` cannot automatically
plot structures, therefore you can use the tool ``RNAplot``. Note that you
**can't** simply pipe the  output of ``RNAsubopt`` to ``RNAplot`` using::

  $ RNAsubopt < test.seq | RNAplot

You need to manually create a file for each structure you want to plot.
Here, for example we created a new file named ``suboptstructure.txt``::

  > suboptstructure-4.20
  CUACGGCGCGGCGCCCUUGGCGA
  (((.((((...))))..)))...

The fasta header is optional, but useful (without it the outputfile will
be named ``rna.ps``).

The next two lines contain the sequence and the suboptimal structure you
want to plot; in this case we plotted the structure with the folding
energy of -4.20.

Then plot it with

.. code::

  $ RNAplot < suboptstructure.txt


Note that the number of suboptimal structures grows exponentially with
sequence length and therefore this approach is only tractable for
sequences with less than 100 nt. To keep the number of suboptimal
structures manageable the option ``--noLP`` can be used, forcing
``RNAsubopt`` to produce only structures without isolated base
pairs. While ``RNAsubopt`` produces *all* structures within an
energy range, ``mfold`` produces only a few, hopefully representative,
structures. Try folding the sequence on the mfold
server at http://mfold.rna.albany.edu/?q=mfold.

Sometimes you want to get information about unusual properties of the
Boltzmann ensemble (the sum of all RNA structures possible) for which no
specialized program exists. For example you want to know all fractions 
of a bacterial mRNA in the Boltzmann ensemble where the Shine-Dalgarno (SD)
sequence is unpaired. If the SD sequence is concealed by secondary
structure the translation efficiency is reduced.

In such cases you can resort to drawing a representative sample of
structures from the Boltzmann ensemble by using the option
``-p``. Now you can simply count how many structures in the sample
possess the feature you are looking for. This number divided by the
size of your sample gives you the desired fraction.

The following example calculates the fraction of structures in the
ensemble that have bases 6 to 8 unpaired.


Sampling the Boltzmann Ensemble
===============================

- Draw a sample of size 10,000 from the Boltzmann ensemble
- Calculate the desired property, e.g. by using a ``perl`` script::

    $ RNAsubopt -p 10000 < test.seq > tt
    $ perl -nle '$h++ if substr($_,5,3) eq "...";
      END {print $h/$.}' tt
      0.391960803919608


A far better way to calculate this property is to use ``RNAfold -p``
to get the ensemble free energy, which is related to the partition
function via :math:`F = -RT\ln(Q)`, for the unconstrained (:math:`F_u`)
and the constrained case (:math:`F_c`), where the three bases are not
allowed to form base pairs (use option ``-C``), and evaluate
:math:`p_c = \exp((F_u - F_c)/RT)` to get the desired probability.

So let's do the calculation using ``RNAfold``::

  $ RNAfold -p
  Input string (upper or lower case); @ to quit
  ....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8
  CUACGGCGCGGCGCCCUUGGCGA
  length = 23
  CUACGGCGCGGCGCCCUUGGCGA
  ...........((((...)))).
   minimum free energy =  -5.00 kcal/mol
  ....{,{{...||||...)}}}.
   free energy of ensemble =  -5.72 kcal/mol
  ....................... {  0.00 d=4.66}
   frequency of mfe structure in ensemble 0.311796; ensemble diversity 6.36  


Now we have calculated the free ensemble energy of the ensemble over all structures
:math:`F_u`, in the next step we have to calculate it for the structures using a
constraint (:math:`F_c`).

Following notation has to be used for defining the constraint:

- ``|`` : paired with another base
- ``.`` : no constraint at all
- ``x`` : base must not pair
- ``<`` : base i is paired with a base j<i
- ``>`` : base i is paired with a base j>i
- matching brackets ``( )``: base i pairs base j


So our constraint should look like this::

  .....xxx...............

Next call the application with following command and provide the sequence
and constraint we just created::

  $ RNAfold -p -C

The output should look like this::

  length = 23
  CUACGGCGCGGCGCCCUUGGCGA
  ...........((((...)))).
   minimum free energy =  -5.00 kcal/mol
  ...........((((...)))).
   free energy of ensemble =  -5.14 kcal/mol
  ...........((((...)))). { -5.00 d=0.42}
   frequency of mfe structure in ensemble 0.792925; ensemble diversity 0.79  

Afterwards evaluate the desired probability according to the formula given before
e.g. with a simple ``perl`` script::

  $ perl -e 'print exp(-(5.72-5.14)/(0.00198*310.15))."\n"'


You can see that there is a slight difference between the ``RNAsubopt`` run with 10,000 
samples and the ``RNAfold`` run including all structures. 
