doxy2swig (improved)
--------------------

Doxygen XML to SWIG docstring converter (improved version).

Converts Doxygen generated XML files into a file containing docstrings
for use by SWIG.


DEPRECATION NOTE
----------------

As of swig 4.0.0, direct support for parsing and conversion of doxygen documentation is available, and according to the swig documentation can be enabled via `swig -python -doxygen` (perhaps in conjunction with `%feature("autodoc", "1")`). It is still possible to generate and include docstrings via `doxy2swig` as described below, but this should no longer be required and is considered deprecated.


Usage
-----

1. Set doxygen to generate XML output by adding `GENERATE_XML = YES` to the project doxygen configuration file (or create a seperate one for this purpose).

2. Run doxygen to extract the documentation and output it in XML format.

3. Use doxy2swig.py to convert the XML documentation to doctrings for SWIG

        doxy2swig.py [options] index.xml output.i

        index.xml is your doxygen generated XML file and output.i is where the
        output will be written (the file will be clobbered).


        Options:
          -h, --help            show this help message and exit
          -f, --function-signature
                                include function signature in the documentation. This
                                is handy when not using swig auto-generated function
                                definitions %feature("autodoc", [0,1])
          -t, --type-info       include type information for arguments in function
                                signatures. This is similar to swig autodoc level 1
          -c, --constructor-list
                                generate a constructor list for class documentation.
                                Useful for target languages where the object
                                construction should be documented in the class
                                documentation.
          -a, --attribute-list  generate an attributes list for class documentation.
                                Useful for target languages where class attributes
                                should be documented in the class documentation.
          -o, --overloaded-functions
                                collect all documentation for overloaded functions.
                                Useful for target languages that have no concept of
                                overloaded functions, but also to avoid having to
                                attach the correct docstring to each function overload
                                manually
          -w W, --width=W       textwidth for wrapping (default: 80). Note that the
                                generated lines may include 2 additional spaces (for
                                markdown).
          -q, --quiet           be quiet and minimize output

4. Include the docstring file generated by doxy2swig.py at the top of your swig file using `%include output.i`.


Note
----

To attach docstrings to classes with SWIG for python and using the `-builtin` option,
a version of SWIG >= 3.0.7 is required.
Without this, the class documentation including constructor lists (`-c`) and attribute
lists (`-a`) will not be available from python.
