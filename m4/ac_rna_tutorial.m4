# RNA_ENABLE_TUTORIAL(PROJECT_NAME, [documentation-output-directory])
#
#
AC_DEFUN([RNA_GET_TUTORIAL],[
  ## collect the subpackages/programs we gonna build
  AS_IF([test "x$with_tutorial" != "xyes"], [ AC_RNA_APPEND_VAR_COMMA($1, [None]) ])
  AS_IF([test "x$with_tutorial_pdf" = "xyes"], [ AC_RNA_APPEND_VAR_COMMA($1, [PDF]) ])
  AS_IF([test "x$with_tutorial_html" = "xyes"], [ AC_RNA_APPEND_VAR_COMMA($1, [HTML]) ])
])

AC_DEFUN([RNA_ENABLE_TUTORIAL],[

RNA_ADD_PACKAGE([tutorial_pdf],
                [PDF tutorial],
                [yes])
RNA_ADD_PACKAGE([tutorial_html],
                [HTML tutorial],
                [no])
RNA_ADD_PACKAGE([tutorial],
                [ViennaRNA tutorial],
                [yes],
                [
                  with_tutorial=no
                  with_tutorial_pdf=no
                  with_tutorial_html=no
                ],
                [])


AC_PATH_PROG(pdflatex,[pdflatex],no)
AC_PATH_PROG(latex,[latex],no)
AC_PATH_PROG(dvipdf,[dvipdf],no)
AC_PATH_PROG(htlatex,[htlatex],no)

# setup everything in order to generate the doxygen configfile
RNA_PACKAGE_IF_ENABLED([tutorial],[

  AC_SUBST([TUTORIAL_PROJECT_NAME], [RNA-tutorial-$PACKAGE_VERSION])
  AC_SUBST([TUTORIAL_SRCDIR], [$srcdir])
  AC_SUBST([TUTORIAL_DOCDIR], [ifelse([$2], [], [RNA-Tutorial], [$2])])

  RNA_PACKAGE_IF_ENABLED([tutorial_pdf], [
    # check whether pdflatex or latex and dvipdf are available
    if test "x$pdflatex" = xno;
    then
      if test "x$latex" = xno;
      then
        AC_MSG_WARN([neither latex or pdflatex exists on your system!])
        AC_MSG_WARN([deactivating automatic (re)generation of tutorial!])
        tutorial_requirements_pdf=no
      else
        if test "x$dvipdf" = xno;
        then
          AC_MSG_WARN([dvipdf command is missing on your system!])
          AC_MSG_WARN([deactivating automatic (re)generation of tutorial!])
          tutorial_requirements_pdf=no
        else
          _latex_cmd=$latex
        fi
      fi
    else
      _latex_cmd=$pdflatex
    fi

    AC_SUBST([TUTORIAL_CMD_LATEX], [$_latex_cmd])
    AC_SUBST([TUTORIAL_CMD_DVIPDF], [$dvipdf])

    # check if a generated tutorial already exists
    if test "x$tutorial_requirements_pdf" = xno;
    then
      AC_RNA_TEST_FILE( [$TUTORIAL_DOCDIR/$TUTORIAL_PROJECT_NAME.pdf],
                        [with_tutorial_pdf=yes],
                        [with_tutorial_pdf=no])
    fi
  ])

  RNA_PACKAGE_IF_ENABLED([tutorial_html], [
    AC_MSG_WARN([

  ============================================================================
   Building the HTML tutorial is intended for the maintainers website only!

   YOU MOST PROBABLY DO NOT WANT TO DO THIS!

   Use the output at your own risk!
   Also note, that the HTML tutorial is ignored during the install step!
  ============================================================================
   ])

    if test "x$htlatex" = xno;
    then
        AC_MSG_WARN([htlatex does not exists on your system!])
        AC_MSG_WARN([deactivating automatic (re)generation of HTML tutorial!])
        with_tutorial_html=no
    fi

    AC_SUBST([TUTORIAL_CMD_HTLATEX], [$htlatex])
  ])

  ## disable tutorial in case PDF and HTML are disabled
  AS_IF([ test "x$with_tutorial_pdf" == "xno" && test "x$with_tutorial_html" == "xno"], [ with_tutorial=no ])
])

AC_SUBST([TUTORIAL_PDF_NAME], [ifelse([$with_tutorial],
                                      [no],
                                      [],
                                      [$TUTORIAL_PROJECT_NAME.pdf])])


AM_CONDITIONAL(WITH_TUTORIAL, test "x$with_tutorial" != xno)
AM_CONDITIONAL(WITH_TUTORIAL_BUILD, test "x$tutorial_requirements_pdf" != xno)
AM_CONDITIONAL(WITH_TUTORIAL_PDF, test "x$with_tutorial_pdf" != xno)
AM_CONDITIONAL(WITH_TUTORIAL_HTML, test "x$with_tutorial_html" != xno)
AM_CONDITIONAL(WITH_TUTORIAL_PDFLATEX, test "x$pdflatex" != xno)
])

