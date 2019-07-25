

#
# This macro checks for the programs
# pdflatex, bibtex, and makeindex and
# substitutes the following variables:
#
# LATEX_CMD
# BIBTEX_CMD
# MAKEINDEX_CMD
#
# If the corresponding programs can't be found,
# these variables are set to the string 'no'
#
AC_DEFUN([RNA_LATEX_ENVIRONMENT], [

    AC_PATH_PROG(pdflatex, [pdflatex], no)
    AC_PATH_PROG(bibtex,[bibtex],no)
    AC_PATH_PROG(makeindex,[makeindex],no)

    AC_SUBST([LATEX_CMD], [$pdflatex])
    AC_SUBST([BIBTEX_CMD], [$bibtex])
    AC_SUBST([MAKEINDEX_CMD], [$makeindex])
])


#
# RNA_LATEX_TEST_PACKAGES([pkg1 pkg2 ...], [ ACTION_IF_FOUND ], [ ACTION_IF_NOT_FOUND ])
#
# Checks for LaTeX packages and executes blocks of code if found,
# or not found, respectively
#
AC_DEFUN([RNA_LATEX_TEST_PACKAGES],[
    AC_REQUIRE([RNA_LATEX_ENVIRONMENT])

    AS_IF([test "x$LATEX_CMD" != "xno"], [
        AC_MSG_CHECKING([for LaTeX package(s) $1])
        ac_latex_package_test_success=no

        # create subdirectory to test LaTeX compilation with list of packages
        rm -rf .latextmp
        mkdir .latextmp
        cd .latextmp
        echo "\\documentclass{article}" > testfile.tex
        m4_foreach_w([latex_pkg], [$1], [
            echo "\\usepackage{latex_pkg}" >> testfile.tex
        ])
        echo "\\begin{document}\\end{document}" >> testfile.tex
        cat testfile.tex | $LATEX_CMD 2>&1 1>/dev/null && ac_latex_package_test_success=yes;
        cd ..
        rm -rf .latextmp
        AC_MSG_RESULT($ac_latex_package_test_success)

        # run user commands as supplied through arguments 2 and 3
        AS_IF([test "x$ac_latex_package_test_success" = "xyes"], [ $2 ], [ $3 ])
    ], [ $3 ])
])

