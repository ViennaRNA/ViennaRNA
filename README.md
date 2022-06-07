[![GitHub release](https://img.shields.io/github/release/ViennaRNA/ViennaRNA.svg)](https://www.tbi.univie.ac.at/RNA/#download)
[![Commits since last version](https://img.shields.io/github/commits-since/ViennaRNA/ViennaRNA/v2.5.0)](https://github.com/ViennaRNA/ViennaRNA/compare/v2.5.0...v2.5.1)
[![Build Status](https://travis-ci.com/ViennaRNA/ViennaRNA.svg?branch=master)](https://travis-ci.com/github/ViennaRNA/ViennaRNA)
[![Github All Releases](https://img.shields.io/github/downloads/ViennaRNA/ViennaRNA/total.svg)](https://github.com/ViennaRNA/ViennaRNA/releases)
[![Conda](https://img.shields.io/conda/v/bioconda/viennarna.svg)](https://anaconda.org/bioconda/viennarna)
[![Conda Downloads](https://img.shields.io/conda/dn/bioconda/viennarna.svg)](https://anaconda.org/bioconda/viennarna)
[![AUR](https://img.shields.io/aur/version/viennarna.svg)](https://aur.archlinux.org/packages/viennarna/)

# ViennaRNA Package

A C code library and several stand-alone programs for the prediction and comparison of RNA secondary structures.

Amongst other things, our implementations allow you to:

- predict minimum free energy secondary structures
- calculate the partition function for the ensemble of structures
- compute various equilibrium probabilities
- calculate suboptimal structures in a given energy range
- compute local structures in long sequences
- predict consensus secondary structures from a multiple sequence alignment
- predict melting curves
- search for sequences folding into a given structure
- compare two secondary structures 
- predict interactions between multiple RNA molecules

The package includes `Perl 5`, `Python 2`, and `Python 3` modules that give
access to almost all functions of the C library from within the respective
scripting languages.

There is also a set of programs for analyzing sequence and distance
data using split decomposition, statistical geometry, and cluster methods.
They are not maintained any more and not built by default.

The code very rarely uses static arrays, and all programs should work for 
sequences up to a length of 32,700 (if you have huge amounts of memory that
is).

See the [NEWS](NEWS) and [CHANGELOG.md](CHANGELOG.md) files for changes between versions.

----

## Table of Contents
1. [Availability](#availability)
2. [Installation](#installation)
3. [Configuration](#configuration)
4. [Executable Programs](#executable-programs)
5. [Energy Parameters](#energy-parameters)
6. [References](#references)
7. [License](#license)
8. [Contact](#contact)

----

## Availability

The most recent source code and documentation should always be available through
the [official ViennaRNA website](https://www.tbi.univie.ac.at/RNA) and through
[github](https://github.com/ViennaRNA/ViennaRNA).

----

## Installation

For best portability the ViennaRNA package uses the GNU autoconf and automake
tools. The instructions below are for installing the ViennaRNA package from
source.

*See the file [INSTALL](INSTALL) for a more detailed description of the build
and installation process.*

#### Quick Start

Usually you'll simply unpack the distribution tarball, configure and make:
```
tar -zxvf ViennaRNA-2.5.1.tar.gz
cd ViennaRNA-2.5.1
./configure
make
sudo make install
```

#### User-dir Installation
If you do not have root privileges on your computer, you might want to install
the ViennaRNA Package to a location where you actually have write access to.
Use the `--prefix` option to set the installation prefix like so:
```
./configure --prefix=/home/username/ViennaRNA
make install
```

This will install everything into a new directory `ViennaRNA` directly into
the home directory of user `username`.

*Note, that the actual install destination paths are listed at the end
of the `./configure` output.*

#### Install from git repository

If you attempt to build and install from our git repository, you need to
perform some additional steps __before__ actually running the `./configure`
script:

1. Unpack the `libsvm` archive to allow for SVM Z-score regression with the
   program `RNALfold`:
```
cd src
tar -xzf libsvm-3.25.tar.gz
cd ..
```

2. Unpack the `dlib` archive to allow for concentration dependency computations with the
   program `RNAmultifold`:
```
cd src
tar -xjf dlib-19.23.tar.bz2
cd ..
```

3. Install the additional maintainer tools `gengetopt`, `help2man`,`flex`,`xxd`,
   and `swig` if necessary. For instance, in RedHat based distributions, the
   following packages need to be installed:
    - `gengetopt` (to generate command line parameter parsers)
    - `help2man` (to generate the man pages) 
    - `flex` and `flex-devel` (to generate sources for RNAforester)
    - `vim-common` (for the `xxd` program)
    - `swig` (to generate the scripting language interfaces)

4. Finally, run the autoconf/automake toolchain:
```
autoreconf -i
```

After that, you can compile and install the ViennaRNA Package as if obtained
from the distribution tarball.

----

## Configuration

This release includes the `RNAforester`, `Kinfold`, `Kinwalker`, and `RNAlocmin`
programs, which can also be obtained as independent packages. Running
`./configure` in the ViennaRNA directory will configure these packages as well.
However, for detailed information and compile time options, see the README and
INSTALL files in the respective subdirectories.

#### Cluster Analysis

The programs `AnalyseSeqs` and `AnalyseDists` offer some cluster analysis tools
(split decomposition, statistical geometry, neighbor joining, Ward's method) for
sequences and distance data. To also build these programs add `--with-cluster`
to your configure options.

#### Kinfold
The `kinfold` program can be used to simulate the folding dynamics of an RNA
molecule, and is compiled by default. Use the `--without-kinfold` option to
skip compilation and installation of Kinfold.

#### RNAforester
The `RNAforester` program is used for comparing secondary structures using tree
alignment. Similar to Kinfold, use the `--without-forester` option to skip
compilation and installation of RNAforester.

#### Kinwalker
The `kinwalker` algorithm performs co-transcriptional folding of RNAs, starting
at a user specified structure (default: open chain) and ending at the minimum
free energy structure. Compilation and installation of this program is
deactivated by default. Use the `--with-kinwalker` option to enable building and
installation of Kinwalker.

#### RNAlocmin
The `RNAlocmin` program is part of the `Basin Hopping Graph` Framework and reads
secondary structures and searches for local minima by performing a gradient walk
from each of those structures. It then outputs an energetically sorted list of
local minima with their energies and number of hits to particular minimum, which
corresponds to a size of a gradient basin. Additional output consists of barrier
trees and Arhenius rates to compute various kinetic aspects. Compilation and
installation of this program is activated by default. Use the
`--without-rnalocmin` option to disable building and installation of RNAlocmin.

#### Scripting Interfaces
The ViennaRNA Package comes with scripting language interfaces for `Perl 5`,
`Python 2`, and `Python 3` (provided by `swig`), that allow one to use the
implemented algorithms directly without the need of calling an executable
program. The necessary requirements are determined at configure time and
particular languages may be deactivated automatically, if the requirements are
not met. You may also switch-off particular languages by passing the
`--without-perl`, `--without-python`, and/or `--without-python2` configure
options, e.g.
```
./configure --without-perl --without-python
```
will turn-off the `Perl 5` and `Python 3` interfaces.

Disabling the entire scripting language support alltogether can be accomplished
with the `--without-swig` switch.

#### Streaming SIMD Extension support
Our latest version contains code that implements a faster multibranch loop
decomposition in global MFE predictions, as used e.g. in `RNAfold`. This
implementation makes use of modern processors capability to execute particular
instructions on multiple data simultaneously (SIMD - single instruction multiple
data, thanks to W. B. Langdon for providing the modified code). Consequently,
the time required to assess the minimum of all multibranch loop decompositions
is reduced up to about one half compared to the runtime of the original
implementation. This feature is enabled by default since version 2.4.11 and a
dispatcher ensures that the correct implementation will be selected at runtime.
If for any reason you want to disable this feature at compile-time use the
following configure flag

```
./configure --disable-simd
```

#### Link Time Optimization (LTO)
To increase the performance of our implementations, the ViennaRNA Package tries
to make use of the Link Time Optimization (LTO) feature of modern C-compilers.
If you are experiencing any troubles at make-time or run-time, or the configure
script for some reason detects that your compiler supports this feature although
it doesn't, you can deactivate it using the flag
```
./configure --disable-lto
```

#### OpenMP support
To enable concurrent computation of our implementations and in some cases
parallelization of the algorithms we make use of the OpenMP API. This interface
is well understood by most modern compilers. However, in some cases it might be
necessary to deactivate OpenMP support and therefore transform RNAlib into a
C-library that is not entirely thread-safe. To do so, add the following
configure option
```
./configure --disable-openmp
```

#### POSIX threads (pthread) support
To enable concurrent computation of multiple input data in `RNAfold`,
`RNAcofold`, `RNAalifold`, and for our implementation of the concurrent
unordered insert, ordered output flush data structure `vrna_ostream_t` we make
use of POSIX threads. This should be supported on all modern platforms and
usually does not pose any problems. In case you want to compile without POSIX
threads support for any reason, add the following configure option
```
./configure --disable-pthreads
```

#### SVM Z-score filter in RNALfold
By default, `RNALfold` that comes with the ViennaRNA Package allows for z-score
filtering of its predicted results using a support vector machine (SVM).
However, the library we use to implement this feature (`libsvm`) is statically
linked to our own RNAlib. If this introduces any problems for your own
third-party programs that link against RNAlib, you can safely switch off the
z-scoring implementation using
```
./configure --without-svm
```

#### GNU Scientific Library
The new program `RNApvmin` computes a pseudo-energy pertubation vector that aims
to minimize the discrepancy of predicted, and observed pairing probabilities.
For that purpose it implements several methods to solve the optimization
problem. Many of them are provided by the GNU Scientific Library, which is why
the `RNApvmin` program, and the RNAlib C-library are required to be linked
against `libgsl`. If this introduces any problems in your own third-party
programs that link against RNAlib, you can turn off a larger protion of
available minimizers in `RNApvmin` and linking against `libgsl` alltogether,
using the switch
```
./configure --without-gsl
```

#### Help
For a complete list of all `./configure` options and important environment
variables, type
```
./configure --help
```

----

## Executable Programs

The ViennaRNA Package includes the following executable programs:

| Program         | Description                                                                                                               |
| --------------- | :-------------------------------------------------------------------------------------------------------------------------|
| `RNA2Dfold`     | Compute MFE structure, partition function and representative sample structures of k,l neighborhoods                       |
| `RNAaliduplex`  | Predict conserved RNA-RNA interactions between two alignments                                                             |
| `RNAalifold`    | Calculate secondary structures for a set of aligned RNA sequences                                                         |
| `RNAcofold`     | Calculate secondary structures of two RNAs with dimerization                                                              |
| `RNAdistance`   | Calculate distances between RNA secondary structures                                                                      |
| `RNAdos`        | Compute the density of states for the conformation space of a given RNA sequence                                          |
| `RNAduplex`     | Compute the structure upon hybridization of two RNA strands                                                               |
| `RNAeval`       | Evaluate free energy of RNA sequences with given secondary structure                                                      |
| `RNAfold`       | Calculate minimum free energy secondary structures and partition function of RNAs                                         |
| `RNAheat`       | Calculate the specific heat (melting curve) of an RNA sequence                                                            |
| `RNAinverse`    | Find RNA sequences with given secondary structure (sequence design)                                                       |
| `RNALalifold`   | Calculate locally stable secondary structures for a set of aligned RNAs                                                   |
| `RNALfold`      | Calculate locally stable secondary structures of long RNAs                                                                |
| `RNAmultifold`  | Compute secondary structures and probabilities for multiple interacting RNAs                                              |
| `RNApaln`       | RNA alignment based on sequence base pairing propensities                                                                 |
| `RNApdist`      | Calculate distances between thermodynamic RNA secondary structures ensembles                                              |
| `RNAparconv`    | Convert energy parameter files from ViennaRNA 1.8 to 2.0 format                                                           |
| `RNAPKplex`     | Predict RNA secondary structures including pseudoknots                                                                    |
| `RNAplex`       | Find targets of a query RNA                                                                                               |
| `RNAplfold`     | Calculate average pair probabilities for locally stable secondary structures                                              |
| `RNAplot`       | Draw RNA Secondary Structures in PostScript, SVG, or GML                                                                  |
| `RNApvmin`      | Calculate a perturbation vector that minimizes discripancies between predicted and observed pairing probabilities         |
| `RNAsnoop`      | Find targets of a query H/ACA snoRNA                                                                                      |
| `RNAsubopt`     | Calculate suboptimal secondary structures of RNAs                                                                         |
| `RNAup`         | Calculate the thermodynamics of RNA-RNA interactions                                                                      |
| `AnalyseSeqs`   | Analyse sequence data                                                                                                     |
| `AnalyseDists`  | Analyse distance matrices                                                                                                 |


A couple of useful utilities can be found in the src/Utils directory.

All executables read from stdin and write to stdout and use command line
switches rather than menus to be easily usable in pipes. For more detailed
information see the man pages. Perl utilities contain POD documentation
that can be read by typing e.g.
```
perldoc relplot.pl
```

Together with this version we also distribute the programs

- `kinfold`,
- `RNAforester`,
- `RNAlocmin`, and
- `kinwalker`

See the README files in the respective sub-directories.

----

## References

If you use our software package, you may want to cite the follwing publications:

- [R. Lorenz et al. (2011)](https://almob.biomedcentral.com/articles/10.1186/1748-7188-6-26),
"ViennaRNA Package 2.0", Algorithms for Molecular Biology, 6:26

- [I.L. Hofacker (1994)](https://link.springer.com/article/10.1007/BF00818163),
"Fast folding and comparison of RNA secondary structures",
Monatshefte fuer Chemie, Volume 125, Issue 2, pp 167-188


*Note, that the individual executable programs state their own list of references
in the corresponding man-pages.*

----

## Energy Parameters

Since version 2.0.0 the build-in energy parameters, also available as parameter
file [rna_turner2004.par](misc/rna_turner2004.par), are taken from:

- [D.H. Mathews et al. (2004)](https://doi.org/10.1073/pnas.0401799101),
"Incorporating chemical modification constraints into a dynamic programming
algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

- [D.H. Turner et al. (2009)](https://dx.doi.org/10.1093/nar/gkp892), 
"NNDB: The nearest neighbor parameter database
for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282.

For backward compatibility we also provide energy parameters from Turner et al.
1999 in the file [rna_turner1999.par](misc/rna_turner1999.par).

Additionally, a set of trained RNA energy parameters from Andronescou et al.
2007, [rna_andronescou2007.par](misc/rna_andronescou2007.par), a set of RNA
energy parameters obtained by graft and grow genetic programming from Langdon
et al. 2018, [rna_langdon2018.par](misc/rna_langdon2018.par), and two sets of
DNA parameters, [dna_mathews1999.par](misc/dna_mathews1999.par) and
[dna_mathews2004.par](misc/dna_mathews2004.par), are included in the package as
well.

----

## License

Please read the copyright notice in the file [COPYING](COPYING)!

If you're a commercial user and find these programs useful, please consider
supporting further developments with a donation.

----

## Contact

We need your feedback! Send your comments, suggestions, and questions to
rna@tbi.univie.ac.at

Ivo Hofacker, Spring 2006
