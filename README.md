
# RNAxplorer
RNAxplorer is a multitool, that offers different methods to explore RNA energy landscapes.
The main use case is sampling of representative structures of the RNA conformation space,
in order to compute RNA folding kinetics. The workflow consists of several steps which 
is depicted in the following figure.

![depiction of the RNA folding kinetics workflow](doc/figures/kinetics_workflow.svg)

The sampling method is the most crucial step, because the number of representative structures determines the runtime of the subsequent tasks.
The RNAxplorer employs efficient dynamic programming based Boltzmann sampling, but improves on previous approaches by adding guiding potentials. These guiding potentials are used to avoid already well-sampled regions of the structure space and steer the sampling towards unexplored regions. 
RNAxplorer offers 3 different sampling approaches, which use either attractive potentials that steer sampling towards regions of interest or repellant potentials that help avoid regions that are already well represented in the sample.

Boltzmann sampling based methods for example produce lots of structures in the vicinity of the minimum free energy structure (MFE). The redundancy of similar structures can be avoided by using a repellant potential in order to increase the energy for structures in a certain region or with certain properties. A repellant potential that penalizes the MFE helps to quickly detect energetically higher minima in the energy  landscape. This is depicted in the following figure.

![depiction of guiding potentials](doc/figures/guidingpotential.svg)

## Installation
### From Linux Package
<table><thead><tr>
<th> Arch </th>
<th> Debian </th>
<th> Ubuntu </th>
<th> openSUSE </th>
</tr></thead><tbody><tr>
<td style="vertical-align:top">
<details><summary>Arch_Extra</summary><p><a href="https://download.opensuse.org/repositories/home:/entzian/Arch_Extra/x86_64/RNAxplorer-0.9.0-1-x86_64.pkg.tar.xz"> RNAxplorer - 0.9.0 - x86_64</a></p>
</details></td>
<td style="vertical-align:top">
<details><summary>Debian_10.0</summary>
<p><a href="https://download.opensuse.org/repositories/home:/entzian/Debian_10/amd64/rnaxplorer_0.9.0_amd64.deb"> rnaxplorer - 0.9.0 - 64 bit</a><br>
<a href="https://download.opensuse.org/repositories/home:/entzian/Debian_10/amd64/python3-rnaxplorer_0.9.0_amd64.deb"> python3-rnaxplorer - 0.9.0 - 64 bit</a>
</p>
<details><summary>Debian_9.0</summary><p>
<a href="https://download.opensuse.org/repositories/home:/entzian/Debian_9.0/i386/rnaxplorer_0.9.0_i386.deb"> rnaxplorer - 0.9.0 - 32 bit</a><br>
<a href="https://download.opensuse.org/repositories/home:/entzian/Debian_9.0/i386/python3-rnaxplorer_0.9.0_i386.deb"> python3-rnaxplorer - 0.9.0 - 32 bit</a>
</p>
<p><a href="https://download.opensuse.org/repositories/home:/entzian/Debian_9.0/amd64/rnaxplorer_0.9.0_amd64.deb"> rnaxplorer - 0.9.0 - 64 bit</a><br>
<a href="https://download.opensuse.org/repositories/home:/entzian/Debian_9.0/amd64/python3-rnaxplorer_0.9.0_amd64.deb"> python3-rnaxplorer - 0.9.0 - 64 bit</a>
</p>
</details></td>
<td style="vertical-align:top">
<details><summary>xUbuntu_19.10</summary><p>
<a href="https://download.opensuse.org/repositories/home:/entzian/xUbuntu_19.10/amd64/rnaxplorer_0.9.0_amd64.deb"> rnaxplorer - 0.9.0 - 64 bit</a><br>
<a href="https://download.opensuse.org/repositories/home:/entzian/xUbuntu_19.10/amd64/python3-rnaxplorer_0.9.0_amd64.deb"> python3-rnaxplorer - 0.9.0 - 64 bit</a>
</p></details>
<details><summary>xUbuntu_19.04</summary><p>
<a href="https://download.opensuse.org/repositories/home:/entzian/xUbuntu_19.04/amd64/rnaxplorer_0.9.0_amd64.deb"> rnaxplorer - 0.9.0 - 64 bit</a><br>
<a href="https://download.opensuse.org/repositories/home:/entzian/xUbuntu_19.04/amd64/python3-rnaxplorer_0.9.0_amd64.deb"> python3-rnaxplorer - 0.9.0 - 64 bit</a>
</p></details>
<details><summary>xUbuntu_18.10</summary><p>
<a href="https://download.opensuse.org/repositories/home:/entzian/xUbuntu_18.10/amd64/rnaxplorer_0.9.0_amd64.deb"> rnaxplorer - 0.9.0 - 64 bit</a><br>
<a href="https://download.opensuse.org/repositories/home:/entzian/xUbuntu_18.10/amd64/python3-rnaxplorer_0.9.0_amd64.deb"> python3-rnaxplorer - 0.9.0 - 64 bit</a>
</p>
</details><details><summary>xUbuntu_18.04</summary><p>
<a href="https://download.opensuse.org/repositories/home:/entzian/xUbuntu_18.04/amd64/rnaxplorer_0.9.0_amd64.deb"> rnaxplorer - 0.9.0 - 64 bit</a><br>
<a href="https://download.opensuse.org/repositories/home:/entzian/xUbuntu_18.04/amd64/python3-rnaxplorer_0.9.0_amd64.deb"> python3-rnaxplorer - 0.9.0 - 64 bit</a>
</p>
</details></td>
<td style="vertical-align:top">
<details><summary>openSUSE_Tumbleweed</summary><p>
<a href="https://download.opensuse.org/repositories/home:/entzian/openSUSE_Tumbleweed/x86_64/RNAxplorer-0.9.0-187.1.x86_64.rpm"> RNAxplorer - 0.9.0 - x86_64</a><br>
<a href="https://download.opensuse.org/repositories/home:/entzian/openSUSE_Tumbleweed/x86_64/python3-rnaxplorer-0.9.0-187.1.x86_64.rpm"> python3-RNAxplorer - 0.9.0 - x86_64</a>
</p>
</details></td>
</tr></tbody></table>



### From [release](https://github.com/ViennaRNA/RNAxplorer/releases) Source

You can download the source tar balls for the individual releases from the [release page](https://github.com/ViennaRNA/RNAxplorer/releases).

To configure, compile and install execute the following commands on your command line:
```
./configure [--help for additional configuration options]
make
make install
```
Dependencies:
  - [ViennaRNA library (>= 2.4.14)](https://www.tbi.univie.ac.at/RNA/#download)
  - [lapacke](http://www.netlib.org/lapack/lapacke.html)


### From github Source

To configure, compile and install execute the following commands on your command line:
```
autoreconf -i
./configure [--help for additional configuration options]
make
make install
```
Dependencies:
  - [ViennaRNA library (>= 2.4.14)](https://www.tbi.univie.ac.at/RNA/#download)
  - [gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html)
  - [lapacke](http://www.netlib.org/lapack/lapacke.html)
  
### With Python Interface
In order to install the python interface for the sampling methods, the following
dependencies are required:
  - [swig](http://www.swig.org/)
  - Python3 and Python3-dev

## Use cases
- repellant sampling with guiding potentials on base pair level
- repellant sampling with guiding potentials on structure base pair distance level
- repellant or attractive sampling with reference structures
- retrieve local minima of secondary structures (via gradient walks in parallel)

## Description
 In default mode (or with -M RSH option) it takes an RNA sequence as input (either stdin or --sequence parameter) 
 and outputs sampled secondary RNA structures. The repellant sampling method iteratively penalizes base pairs of 
 local minima of structures that have been seen too often. This results in a diverse sample set with the most
 important low free energy structures.
 A second mode can be used with --penalize-structures, which uses a different kind of repellant guiding potentials.
 Here the guiding potentials of depend on the inherent base pair distances between loop decompositions of overrepresented structures 
 (More details can be found in the supplementary material of the corresponding publication).
 
 Another important sampling method (-M RS option) is based on reference structures (--struc1
 and --struc2). This method produces structures in the vicinity of these two reference
 structures. Arbitrary many references can be added if a fasta file is used as input
 (via stdin).
 
 Often the output of sampling methods has to be coarse grained by local minima that are defined
 by a gradient walk. A parallelized gradient descent procedure can be used to retrieve
 local minima (-M RL option) of sampled structures (input via stdin).
 
## Examples
### Repellant Sampling Method (penalize base pairs)
In order to use this method, we need a sequence or fastafile (text file) with the following content:
``` 
>sv11
GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
```

In order to produce 10 sampled structures call either
```
RNAxplorer -M RSH -n 10 --sequence GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
```
or
```
cat sv11.fasta | RNAxplorer -M RSH -n 10
```

The output is the the samples and the unique local minima. The columns contain the index, structure, free energy and how often this minimum was reached from a sampled structure:
``` 
GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
(null)
(null)
samples so far:     10 /     10 ... done
Samples: 
(((..((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..))))))))...))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..))))))))).)))..
(((..((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..))))))...))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..)))))))).)))..
(((.((((((((.((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...))))))))))))))))))...))))))))..))).
(((.(((((((((((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))..))))))))..))).
(((..((((((.(((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))...))))))...))).
Local minima: 
     GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
   0 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))). -97.70      6
   1 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))). -97.00      1
   2 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..)))))))))..))). -96.70      1
   3 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..)))))))).))).. -96.10      1
   4 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..))))))))).))).. -95.10      1

```

The output can also be stored in separate files:
```
cat sv11.fasta | RNAxplorer -M RSH -n 10 --lmin-file=repellant_sampling.txt
```
This creates two files with the following content:
```
repellant_sampling.txt:
     GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
   0 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))). -97.70      6
   1 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))). -97.00      2
   2 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..)))))))))..))). -96.70      1
   3 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..))))))))).))).. -95.10      1
```
```
repellant_sampling.samples: 
     GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..)))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((....(((((((((((((..((((...))))..)))))))))))))....)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))).
((..(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))..)))))))))..))..
(((.(((((((((((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))..))))))))..))).
```

### Repellant Sampling Method (penalize structures)
The second sampling methods is based on guiding potentials that penalize structres, which have been seen too often, in terms of inherent base pair distances.
This method can be used with the following command:
```
cat sv11.fasta | RNAxplorer -M RSH -n 10 --penalize-structures --lmin-file=repellant_sampling.txt
```
Again the output is stored in the following files
```
repellant_sampling.txt:
     GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
   0 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))). -97.70      7
   1 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))). -97.00      2
   2 .((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))...)). -94.10      1
```

```
repellant_sampling.samples:
     GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
.((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))...)).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..))))))))))))))..))))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).
.((.((((((((((((((((((((((((((..(..(((((((((((((..((((...))))..)))))))))))))..)..)))))))))))))))))))..)))))))..))..
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))))..)))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((...(((...)))...)))))))))))))).).)))))))))))))))))))..))))))))..))).
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).

```

### Sampling with References
In order to produce structures in the vicinity of certain reference structures, you should first add the references to a fasta file:

```
>sv11
GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....)))))))).....
(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).

```
Then you call the following method:
```
cat sv11.fasta | RNAxplorer -M SM -e N -i 10
distortions: d_x0 = 0.3686046157, d_x1 = 0.0000000000 
1	85	-64.80	(1) ((((.((((((...))))))...((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....)))))))).....
2	88	-65.60	(2) (((((((((((...)))))))..((((((((((....))))))))))((....)))))).....((((((((........)))))))).((((((((.....)))))))).....
2	88	-65.60	(2) (((((((((((...)))))))..((((((((((....))))))))))((....)))))).....((((((((........)))))))).((((((((.....)))))))).....
11	77	-64.00	(1) (((.(((((((...)))))))(.((((((((((....)))))))))))...(((...)))....((((((((........)))))))).((((((((.....)))))))).))).
84	2	-95.40	(1) (((.(((((((.(((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))...)))))))..))).
84	4	-92.80	(1) (((.((((((((((((((((((((((((((..(.((((((((((((((..(((.....)))..)))))))))))))).)..)))))))))))))))))..)))))))))..))).
85	1	-96.70	(3) (((.(((((((((((((((((((((((((((...((((((((((((((..((((...))))..))))))))))))))...)))))))))))))))))))..))))))))..))).
85	1	-95.90	(3) (((..((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..)))))))...))).
85	1	-96.10	(3) (((.((((((((.((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).))))))))))))))))))...))))))))..))).
85	5	-94.10	(1) .((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))...)).
```
The input parameter `-M SM` selects the attractive guiding potential sampling method. With `-e N` the default distortion is used, which 
means that both reference structures have the same probability as the mfe structure. The parameter `-i 10` sets 10 iterations (1 structure per iteration will be drawn).
In order to vary the distortion in each iteration you can choose other values for `-e`. This can produce more structures on the path between two references. The
number of produced samples is then the base pair distance between both references times the number of iterations (this would be the case with `-e S`).

The output contains the computed distortions for both references. The next lines are 5 columns with the following content: base pair distance to the first reference,
base pair distance to the second reference, free energy, number of structures with the same distances to both references, the structure.
In this example you can see that samples in the base pair vicinity of both references has been constructed. In contrast to that, pure Boltzmann sampling would
produce only structures which are similar to the second reference.


### retrieve local minima
Local minima via gradient walks can be produced with the option `-M RL`. Input is a file with the sequence in the first line and structures in the next lines.

```
cat sv11.fasta | RNAxplorer -M RL
```

The output consists of the sequence, the index of the corresponding structure in the input file, the local minima and the energies:
```
GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA
2 (((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))). -97.70
1 (((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....))))))))..... -66.00
```


You can adjust the output and speed up the computation by using additional command line parameters. All parameters are shown by
```
RNAxplorer --help
```

## Post-processing
A usual post processing of the sampled structures is clustering or the rate computation (e.g. via determining the partition funktions of the contact surfaces or 
via lowest saddles of paths between structures). The next step is the computation of of the RNA folding kinetics via numerical integration
of the rate matrix (treekin).

### recommended tools
- [clustering](https://github.com/ViennaRNA/RNAxplorer/scripts/pipeline/clusteralgorithms/diana.py) to reduce the number of representative structures (e.g. via the diana clustering script within RNAxplorer's scripts folder).
- [pourRNA](https://github.com/ViennaRNA/pourRNA/)
- [BHG](https://www.tbi.univie.ac.at/software/BHG/BHGbuilder.html#download)
- [FindPath](https://www.tbi.univie.ac.at/RNA/index.html#download)
  In the source folder of the viennaRNA package:
  ```
  viennarna/src/ViennaRNA/findpath.c
  ```
- [Treekin](https://www.tbi.univie.ac.at/RNA/Treekin/)



 
