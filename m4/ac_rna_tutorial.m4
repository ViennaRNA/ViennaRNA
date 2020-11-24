# RNA_ENABLE_TUTORIAL(PROJECT_NAME, [documentation-output-directory])
#
#
AC_DEFUN([RNA_ENABLE_TUTORIAL],[

    AC_REQUIRE([RNA_LATEX_ENVIRONMENT])

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


AC_PATH_PROG(htlatex,[htlatex],no)

# setup everything in order to generate the doxygen configfile
RNA_PACKAGE_IF_ENABLED([tutorial],[

  AC_SUBST([TUTORIAL_PROJECT_NAME], [RNA-tutorial-$PACKAGE_VERSION])
  AC_SUBST([TUTORIAL_SRCDIR], [$srcdir])
  AC_SUBST([TUTORIAL_DOCDIR], [ifelse([$2], [], [RNA-Tutorial], [$2])])

  RNA_PACKAGE_IF_ENABLED([tutorial_pdf], [
    # check whether pdflatex is available
    if test "x$LATEX_CMD" = xno;
    then
      AC_MSG_WARN([Could not find pdflatex!])
      AC_MSG_WARN([deactivating automatic (re)generation of tutorial!])
      _tutorial_pdf_failed="pdflatex command is missing!"
      tutorial_requirements_pdf=no
    else
      RNA_LATEX_TEST_PACKAGES([url color caption amssymb amsmath upquote verbatim keystroke fancyvrb hyperref graphics pgf xcolor], [], [
        _tutorial_pdf_failed="Required LaTeX packages are missing!"
        tutorial_requirements_pdf=no
      ])
    fi

    AC_SUBST([TUTORIAL_CMD_LATEX], ["$LATEX_CMD -interaction=nonstopmode -halt-on-error"])

    # check if a generated tutorial already exists
    if test "x$tutorial_requirements_pdf" = xno;
    then
      AC_RNA_TEST_FILE( [$TUTORIAL_DOCDIR/$TUTORIAL_PROJECT_NAME.pdf],
                        [with_tutorial_pdf=yes],
                        [with_tutorial_pdf=no
                         tutorial_pdf_failed="($_tutorial_pdf_failed)"])
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
        tutorial_html_failed="(htlatex command is missing!)"
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
])

