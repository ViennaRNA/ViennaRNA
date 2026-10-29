# Contributing to ViennaRNA

## Reporting Bugs

The ViennaRNA Package is developed by humans and consequently may contain bugs that
prevent proper operation of the implemented algorithms. If you believe you have found
any of those nasty animals, please help us to improve our software by sending comments
and bug reports either as issue at [GitHub](https://github.com/ViennaRNA/ViennaRNA/issues)
or by email to rna@tbi.univie.ac.at.

> If in doubt our program is right, nature is at fault.

**Adhere to the following rules while reporting any issues:**

* In your report, **include as much information as possible, such that we are able
  to reproduce it**. If possible, find a minimal example that triggers the bug.
* **Include the version number** for the ViennaRNA Package you experience the bug with.
* Include at least some minimal information regarding your operating system (Linux,
  Mac OS X, Windows, etc.)

## Pull Request Process

The ViennaRNA Package is open-source software, which means that everybody can
have a closer look into our implementations to understand and potentially extend
it's functionality. If you implemented any novel feature into the ViennaRNA
Package that might be of interest to a larger community, please don't hesitate to ask
for merging of your code into our official source tree.

::: {note}
If you wish to contribute to this project, please first discuss any proposed changes
with the owners and main developers. You may do that either through making an issue
at [our official GitHub presence](https://github.com/ViennaRNA/ViennaRNA),
[by email](mailto:rna@tbi.univie.ac.at), or any other personal communication with the
core developer team.
:::

::: {important}
More importantly, if you wish to contribute any files or software, you need to agree
to our ViennaRNA Package Contributors License Agreement (CLA)! Otherwise, your
contributions can't be merged into our source tree. [See below](#contributors-license-agreement)
for further information and the full CLA details.
:::

**Adhere to the following rules, otherwise your pull request can't be considered:**

* Ensure that you **do not check-in any files that are automatically build**!
* **Use descriptive commit messages**!
* Make **atomic commits**, i.e. one feature at a time! Do not commit too many changes at once
  but split them into individual commits.
* **Do not commit cosmetic changes!** Commits that merely beautify the code or fix
   spelling in the documentation will not be considered for a merge!
* When contributing C source code, **follow our code formatting style**. Have a look into
  our source code to get an idea how we format everything. Also, you may use the tool
  `uncrustify` together with our config located in `misc/uncrustify.cfg` to enforce our
  code formatting to your implementation.
* **Only expose symbols** (functions, variables, etc.) to the libraries interface **that are
  absolutely necessary**! Hide all other symbols in the corresponding object file(s) by
  declaring them as `static`.
* **Use the prefix `vrna_`** for any symbol you add to the API of our library! Preprocessor
  macros in header files require the prefix in capital letters, i.e. `VRNA_`.
* **Use C-style comments** at any place necessary to make sure your implementation can still
  be understood and followed in the future.
* **Add test cases** for any new implementation! The test suite is located in the `tests`
  directory and is split into tests for the C-library, executable programs, and the
  individual scripting language interfaces.
* Run `make check` to **ensure that all other test suites still run properly** with your
  applied changes!
* When contributing via GitHub, make a personal fork of our project and create a separate
  branch for your changes. Then make a **pull request to our `user-contrib` branch**.
  Pull requests to the `master` branch will be rejected to keep its history clean.
* Pull requests that have been successfully merged into the `user-contrib` branch
  usually find their way into the next release of the ViennaRNA Package. However, please
  note that the core developers may decide to include your changes in a later version.

## Contributors License Agreement

Thank you for your interest in contributing to the ViennaRNA Package ("We" or "Us").

Before contributing, please note that we adopted a standard Contributors License
Agreement (CLA) agreement provided by [Project Harmony](https://www.harmonyagreements.org),
a community-centered group focused on contributor agreements for free and open source
software (FOSS).

This contributor agreement ("Agreement") documents the rights granted by contributors
to Us. To make this document effective, please sign it and send it to Us by email to
rna@tbi.univie.ac.at.

The respective CLA PDF documents are available in the `doc/CLA` directory of
the distribution tarball, and online at our
[official ViennaRNA Website](https://www.tbi.univie.ac.at/RNA/contributors.html):

[ViennaRNA Package - Individual Contributor License Agreement](doc/CLA/ViennaRNA-CLA-Individual.pdf)

[ViennaRNA Package - Entity Contributor License Agreement](doc/CLA/ViennaRNA-CLA-Entity.pdf)



