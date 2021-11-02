# Contributing to the ViennaRNA Package

## Contents
- [General Remarks](#general-remarks)
- [Reporting Bugs](#reporting-bugs)
- [Pull Request Process](#pull-request-process)
- [Contributors License Agreement (CLA)](#contributors-license-agreement)

## General Remarks

The ViennaRNA Package is developed by humans and consequently may contain bugs that
prevent proper operation of the implemented algorithms. If you think you have found
any of those nasty animals, please help us to improve our software by [reporting the
bug](#reporting-bugs) to us.

The ViennaRNA Package also is open-source software, which means that everybody can
have a closer look into our implementations to understand and potentially extend
it's functionality. If you implemented any novel feature into the ViennaRNA
Package that might be of interest to a larger community, please don't hesitate to ask
for merging of your code into our official source tree. See the [Pull Request Process
section](#pull-request-process) below to find information on how to do that.

Please note that we have a code of conduct. Please follow it in all your interactions
with this project.

If you wish to contribute to this project, please first discuss any proposed changes
with the owners and main developers. You may do that either through making an issue
at [our official GitHub presence](https://github.com/ViennaRNA/ViennaRNA),
[by email](mailto:rna@tbi.univie.ac.at), or any other personal communication with the
core developer team.

More importantly, if you wish to contribute any files or software, you need to agree
to our ViennaRNA Package Contributors License Agreement (CLA)! Otherwise, your
contributions can't be merged into our source tree. [See below](#contributors-license-agreement)
for further information and the full CLA details.

## Reporting Bugs

1. Please make an issue at GitHub or notify us by emailing to rna@tbi.univie.ac.at
2. In your report, include as much information as possible, such that we are able
   to reproduce it. If possible, find a minimal example that triggers the bug.
3. Include the version number for the ViennaRNA Package you experience the bug with.
4. Include at least some minimal information regarding your operating system (Linux,
   Mac OS X, Windows, etc.)

## Pull Request Process

1. Ensure that you have not checked-in any files that are automatically build!
2. When contributing C source code, follow our code formatting guide lines. You
   may use the tool `uncrustify` together with our config located in `misc/uncrustify.cfg`
   to accomplish that.
3. Only expose symbols (functions, variables, etc.) to the libraries interface that are
   absolutely necessary! Hide all other symbols in the corresponding object file(s) by
   declaring them as `static`.
4. Use the prefixes `vrna_` for any symbol you add to the API of our library! Preprocessor
   macros in header files require the prefix in capital letters, i.e. `VRNA_`.
5. Use C-style comments at any place necessary to make sure your implementation can still
   be understood and followed in the future.
6. Add test cases for any new implementation! The test suite is located in the `tests`
   directory and is split into tests for the C-library, executable programs, and the
   individual scripting language interfaces.
7. Run `make check` to ensure that all other test suites still run properly with your
   applied changes!
8. When contributing via GitHub, make a personal fork of our project and create a separate
   branch for your changes. Then make a pull request to our `user-contrib` branch.
   Pull requests to the `master` branch will be rejected to keep its history clean.
9. Pull requests that have been successfully merged into the `user-contrib` branch
   usually find their way into the next release of the ViennaRNA Package. However, please
   note that the core developers may decide to include your changes in a later version.

## Contributors License Agreement

Thank you for your interest in contributing to the ViennaRNA Package ("We" or "Us").

Before contributing, please note that we adopted a standard Contributors License
Agreement (CLA) agreement provided by [Project Harmony](www.harmonyagreements.org),
a community-centered group focused on contributor agreements for free and open source
software (FOSS).

This contributor agreement ("Agreement") documents the rights granted by contributors
to Us. To make this document effective, please sign it and send it to Us by email to
rna@tbi.univie.ac.at.

The respective CLA PDF documents are available in the [doc/CLA directory](doc/CLA) of
the distribution tarball, and online at our
[official ViennaRNA Website](https://www.tbi.univie.ac.at/RNA/contributing.html).
