AC_DEFUN([RNA_CHECK_BUILD_REQUIREMENTS], [

## Check for any tool required to build source files
##
## These tools should only be required for maintainer builds,
## since the generated files are usually distributed along
## with the distribution tar ball
RNA_CHECK_SRC_BUILDERS

## Check for presence of (or capability to generate) the
## Postscript templates we use for RNA secondary structure-
## and Dot plots.
RNA_CHECK_POSTSCRIPT_TEMPLATES

## Check for presence of (or capability to generate) the
## SVG templates we use for RNA secondary structure-
## and Dot plots.
RNA_CHECK_SVG_TEMPLATES

RNA_CHECK_PARAMETER_FILES

RNA_CHECK_DLIB

RNA_CHECK_SWIG_SVM
])


AC_DEFUN([RNA_CHECK_SRC_BUILDERS], [
    AC_ARG_VAR([XXD], [the 'xxd' program to convert text files into C include header files])
    AC_PATH_PROG([XXD], [xxd], [no])
    AC_SUBST([XXD])
    AM_CONDITIONAL(VRNA_AM_SWITCH_HAS_XXD, test "x$XXD" != "xno")

    AC_ARG_VAR([GENGETOPT], [the 'gengetopt' program to generate command line argument parsers for executable programs])
    AC_PATH_PROG([GENGETOPT], [gengetopt], [no])
    AC_SUBST([GENGETOPT])
    AM_CONDITIONAL(VRNA_AM_SWITCH_HAS_GENGETOPT, test "x$GENGETOPT" != "xno")

    AC_ARG_VAR([HELP2MAN], [the 'help2man' script to generate man pages from command line options of our executable programs])
    AC_PATH_PROG([HELP2MAN], [help2man], [no])
    AC_SUBST([HELP2MAN])
    AM_CONDITIONAL(VRNA_AM_SWITCH_BUILD_MANPAGES, test "x$HELP2MAN" != "xno" && test "x$GENGETOPT" != "xno")
])


AC_DEFUN([RNA_CHECK_POSTSCRIPT_TEMPLATES], [
    STATIC_FILE_DIR="${srcdir}/src/ViennaRNA/static"
    POSTSCRIPT_FILE_LIST="$STATIC_FILE_DIR/ps_templates.txt"

    ## load list of postscript template files and replace '\n' by ' '
    PS_TEMPLATE_FILES=`cat $POSTSCRIPT_FILE_LIST | tr '\012' ' '`
    ## create list of hex postscript template files
    PS_TEMPLATE_FILES_HEX=`AS_ECHO("$PS_TEMPLATE_FILES") | sed 's/\.ps/\.hex/g'`

    if test "x$XXD" = "xno"
    then
        for template in $PS_TEMPLATE_FILES_HEX
        do
            AC_RNA_TEST_FILE($STATIC_FILE_DIR/$template,[],[
                AC_MSG_ERROR([
=================================================
Can't find the postscript hex template

${template}

Make sure you've installed the 'xxd' tool to
generate it from source!
=================================================
])
            ])
        done
    fi

    # prepare substitution string for
    # templates_postscript.h file
    PS_TEMPLATE_CONST=""
    for template in $PS_TEMPLATE_FILES_HEX
    do
      # remove the 'postscript/' prefix
      template_name=`AS_ECHO("$template") | sed 's/postscript\///g'`
      # remove the trailing .hex
      template_name=`AS_ECHO("$template_name") | sed 's/\.hex//g'`

      # create a C variable defintion for the template
      # note [[]] will turn into [] after M4 processed everythin
      PS_TEMPLATE_CONST="$PS_TEMPLATE_CONST
static const unsigned char PS_$template_name[[]] = {
#include \"$template\"
};
"
    done

    # Add templates_postscript.h to the files to be processed by
    # the configure script
    AC_CONFIG_FILES([src/ViennaRNA/static/templates_postscript.h])

    # substitute C variable definitions
    AC_SUBST(PS_TEMPLATE_CONST)
    # hack to avoid placing the multiline PS_TEMPLATE_CONST into any Makefile
    _AM_SUBST_NOTMAKE(PS_TEMPLATE_CONST)

    # substitute file list for static/Makefile.am
    AC_SUBST(PS_TEMPLATE_FILES)
    AC_SUBST(PS_TEMPLATE_FILES_HEX)
])

AC_DEFUN([RNA_CHECK_SVG_TEMPLATES], [
    STATIC_FILE_DIR="${srcdir}/src/ViennaRNA/static"
    SVG_FILE_LIST="$STATIC_FILE_DIR/svg_templates.txt"

    ## load list of svg template files and replace '\n' by ' '
    SVG_TEMPLATE_FILES=`cat $SVG_FILE_LIST | tr '\012' ' '`
    ## create list of hex postscript template files
    SVG_TEMPLATE_FILES_HEX=`AS_ECHO("$SVG_TEMPLATE_FILES") | sed 's/\.svg/\.hex/g'`

    if test "x$XXD" = "xno"
    then
        for template in $SVG_TEMPLATE_FILES_HEX
        do
            AC_RNA_TEST_FILE($STATIC_FILE_DIR/$template,[],[
                AC_MSG_ERROR([
=================================================
Can't find the svg hex template

${template}

Make sure you've installed the 'xxd' tool to
generate it from source!
=================================================
])
            ])
        done
    fi

    # prepare substitution string for
    # templates_svg.h file
    SVG_TEMPLATE_CONST=""
    for template in $SVG_TEMPLATE_FILES_HEX
    do
      # remove the 'svg/' prefix
      template_name=`AS_ECHO("$template") | sed 's/svg\///g'`
      # remove the trailing .hex
      template_name=`AS_ECHO("$template_name") | sed 's/\.hex//g'`

      # create a C variable defintion for the template
      # note [[]] will turn into [] after M4 processed everythin
      SVG_TEMPLATE_CONST="$SVG_TEMPLATE_CONST
static const char SVG_$template_name[[]] = {
#include \"$template\"
};
"
    done

    # Add templates_svg.h to the files to be processed by
    # the configure script
    AC_CONFIG_FILES([src/ViennaRNA/static/templates_svg.h])

    # substitute C variable definitions
    AC_SUBST(SVG_TEMPLATE_CONST)
    # hack to avoid placing the multiline SVG_TEMPLATE_CONST into any Makefile
    _AM_SUBST_NOTMAKE(SVG_TEMPLATE_CONST)

    # substitute file list for static/Makefile.am
    AC_SUBST(SVG_TEMPLATE_FILES)
    AC_SUBST(SVG_TEMPLATE_FILES_HEX)
])

AC_DEFUN([RNA_CHECK_PARAMETER_FILES], [
    STATIC_FILE_DIR="${srcdir}/src/ViennaRNA/static"
    PARAMETER_FILE_LIST="${srcdir}/misc/parameter_files.txt"

    ## load list of energy parameter files and replace '\n' by ' '
    PARAMETER_FILES=`cat $PARAMETER_FILE_LIST | sed 's/^/misc\//' | tr '\012' ' '`
    ## create list of hex energy parameter files
    PARAMETER_FILES_HEX=`AS_ECHO("$PARAMETER_FILES") | sed -E 's/\.par|\.json/\.hex/g'`

    if test "x$XXD" = "xno"
    then
        for parfile in $PARAMETER_FILES_HEX
        do
            AC_RNA_TEST_FILE($STATIC_FILE_DIR/$parfile,[],[
                AC_MSG_ERROR([
=================================================
Can't find the energy parameter hex file

${parfile}

Make sure you've installed the 'xxd' tool to
generate it from source!
=================================================
])
            ])
        done
    fi

    # prepare substitution string for
    # energy_parameter_sets.h file
    ENERGY_PARAMETER_CONST=""
    for parfile in $PARAMETER_FILES_HEX
    do
      # remove the 'parameters/' prefix
      parfile_name=`AS_ECHO("$parfile") | sed 's/misc\///g'`
      # remove the trailing .hex
      parfile_name=`AS_ECHO("$parfile_name") | sed 's/.hex//g'`

      # create a C variable defintion for the template
      # note [[]] will turn into [] after M4 processed everythin
      ENERGY_PARAMETER_CONST="$ENERGY_PARAMETER_CONST
static const unsigned char parameter_set_$parfile_name[[]] = {
#include \"$parfile\"
};
"
      # create a SWIG Python output typemap for the template
      # note [[]] will turn into [] after M4 processed everythin
      SWIG_ENERGY_PARAMETER_CONST_PYTHON="$SWIG_ENERGY_PARAMETER_CONST_PYTHON
%typemap(varout) const unsigned char parameter_set_$parfile_name[[]] {
  std::string str( parameter_set_$parfile_name, parameter_set_$parfile_name + sizeof (parameter_set_$parfile_name) / sizeof (parameter_set_$parfile_name[[0]]) );
  swig_result = PyUnicode_FromString(str.c_str());
}
"
      # create a SWIG Perl 5 output typemap for the template
      # note [[]] will turn into [] after M4 processed everythin
      SWIG_ENERGY_PARAMETER_CONST_PERL5="$SWIG_ENERGY_PARAMETER_CONST_PERL5
%typemap(varout) const unsigned char parameter_set_$parfile_name[[]] {
  std::string str( parameter_set_$parfile_name, parameter_set_$parfile_name + sizeof (parameter_set_$parfile_name) / sizeof (parameter_set_$parfile_name[[0]]) );
  sv_setpv(swig_result, str.c_str());
}
"
    done

    SWIG_ENERGY_PARAMETER_CONST_PYTHON=`AS_ECHO("$SWIG_ENERGY_PARAMETER_CONST_PYTHON") | sed 's/swig_result/$result/g'`
    SWIG_ENERGY_PARAMETER_CONST_PERL5=`AS_ECHO("$SWIG_ENERGY_PARAMETER_CONST_PERL5") | sed 's/swig_result/$result/g'`

    # Add templates_postscript.h to the files to be processed by
    # the configure script
    AC_CONFIG_FILES([src/ViennaRNA/static/energy_parameter_sets.h])
    AC_CONFIG_FILES([interfaces/parameter_sets.i])

    # substitute C variable definitions
    AC_SUBST(ENERGY_PARAMETER_CONST)
    AC_SUBST(SWIG_ENERGY_PARAMETER_CONST_PYTHON)
    AC_SUBST(SWIG_ENERGY_PARAMETER_CONST_PERL5)
    # hack to avoid placing the multiline ENERGY_PARAMETER_CONST into any Makefile
    _AM_SUBST_NOTMAKE(ENERGY_PARAMETER_CONST)
    _AM_SUBST_NOTMAKE(SWIG_ENERGY_PARAMETER_CONST_PYTHON)
    _AM_SUBST_NOTMAKE(SWIG_ENERGY_PARAMETER_CONST_PERL5)

    # substitute file list for static/Makefile.am
    AC_SUBST(PARAMETER_FILES)
    AC_SUBST(PARAMETER_FILES_HEX)
])


AC_DEFUN([RNA_CHECK_DLIB], [
  AX_CXX_COMPILE_STDCXX(17, noext, mandatory)
  DLIB_VERSION=19.24
  DLIB_DIR="dlib-${DLIB_VERSION}"
  DLIB_PATH="${srcdir}/src/${DLIB_DIR}"
  DLIB_SRC_FILE="${DLIB_PATH}/dlib/all/source.cpp"

  AC_RNA_TEST_FILE($DLIB_SRC_FILE,[],[
    AC_MSG_ERROR([
=================================================
Can't find dlib's source.cpp

Make sure you've unpacked 'src/dlib-${DLIB_VERSION}.tar.bz2'!

Usually, you only need to execute the following command:

tar -xjf src/dlib-${DLIB_VERSION}.tar.bz2 -C src/
=================================================
])
  ])

  AC_SUBST(DLIB_CPPFLAGS, "-I\$(top_srcdir)/src/${DLIB_DIR} -DDLIB_NO_GUI_SUPPORT")
  AC_SUBST(DLIB_DIR)

])

##
## Check whether the user deactivated libsvm support and still wants to
## build scripting language interfaces. In such cases, the user is required
## to have swig installed since the shipped swig wrappers default to libRNA
## with libsvm support
##
AC_DEFUN([RNA_CHECK_SWIG_SVM], [
  RNA_PACKAGE_IF_DISABLED([svm],[
      AS_IF([test "x$with_swig" == "xyes" && test "x$has_swig" != "xyes"], [
    AC_MSG_ERROR([
=================================================
Compilation requirements missing!

You deactivated SVM support but this requires the
scripting language interface wrappers (Perl 5, Python)
to be re-generated!

Please either install the 'swig' program to enable
re-generation of the respective wrapper files or:

a) deactivate the scripting language interfaces
   alltogether using the '--without-swig' configure
   option.
b) leave SVM support enabled by omitting the
   '--without-svm' configure option
=================================================
])        
      ])
    ])
  ])
])
