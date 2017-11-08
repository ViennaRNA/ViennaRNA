# Contributing

If you wish to contribute to this project, please first discuss any proposed changes
with the owners and main developers. You may do that either through making an issue
at our official Github presence https://github.com/ViennaRNA/ViennaRNA, by email
(rna@tbi.univie.ac.at), or any other personal communication with the core developer
team.

Please note that we have a code of conduct. Please follow it in all your interactions
with this project.

## Reporting Bugs

1. Please make an issue at Github or notify us by emailing to rna@tbi.univie.ac.at
2. In your report, include as much information as possible, such that we are able
   to reproduce it. If possible, find a minimal example that triggers the bug.
3. Include the version number for the ViennaRNA Package you experience the bug with.
4. Include at least some minimal information regarding your operating system (Linux,
   MacOS X, Windows, etc.)

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
8. When contributing via Github, make a personal fork of our project and create a separate
   branch for your changes. Then make a pull request to our `user-contrib` branch.
   Pull requests to the `master` branch will be rejected to keep its histroy clean.
9. Pull requests that have been successfully merged into the `user-contrib` branch
   usually find their way into the next release of the ViennaRNA Package. However, please
   note that the core developers may decide to include your changes in a later version.
