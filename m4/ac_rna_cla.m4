# RNA_ENABLE_CLA(PROJECT_NAME, [documentation-output-directory])
#
#
AC_DEFUN([RNA_ENABLE_CLA],[

AC_REQUIRE([RNA_LATEX_ENVIRONMENT])

RNA_ADD_PACKAGE([cla_pdf],
                [PDF Contributors License Agreement],
                [yes])
RNA_ADD_PACKAGE([cla],
                [ViennaRNA Contributors License Agreement],
                [yes],
                [
                  with_cla=no
                  with_cla_pdf=no
                ],
                [])


# setup everything in order to generate the doxygen configfile
RNA_PACKAGE_IF_ENABLED([cla],[

  AC_SUBST([CLA_I_PROJECT_NAME], [ViennaRNA-CLA-Individual])
  AC_SUBST([CLA_E_PROJECT_NAME], [ViennaRNA-CLA-Entity])
  AC_SUBST([CLA_SRCDIR], [$srcdir])
  AC_SUBST([CLA_DOCDIR], [ifelse([$2], [], [doc/CLA], [$2])])

  RNA_PACKAGE_IF_ENABLED([cla_pdf], [
    # check whether pdflatex is available
    if test "x$LATEX_CMD" = xno;
    then
      AC_MSG_WARN([Could not find pdflatex!])
      AC_MSG_WARN([deactivating automatic (re)generation of CLA!])
      _cla_pdf_failed="pdflatex command is missing!"
      cla_requirements_pdf=no
    else
      RNA_LATEX_TEST_PACKAGES([babel fancyhdr geometry hyperref inputenc lastpage xcolor], [], [
        _cla_pdf_failed="Required LaTeX packages are missing!"
        cla_requirements_pdf=no
      ])
    fi

    AC_SUBST([CLA_CMD_LATEX], ["$LATEX_CMD -interaction=nonstopmode -halt-on-error"])

    # check if a generated CLA already exists
    if test "x$cla_requirements_pdf" = xno;
    then
      AC_RNA_TEST_FILE( [$CLA_DOCDIR/$CLA_I_PROJECT_NAME.pdf],
                        [with_cla_pdf=yes],
                        [with_cla_pdf=no
                         cla_pdf_failed="($_cla_pdf_failed)"])
      AC_RNA_TEST_FILE( [$CLA_DOCDIR/$CLA_E_PROJECT_NAME.pdf],
                        [with_cla_pdf=yes],
                        [with_cla_pdf=no
                         cla_pdf_failed="($_cla_pdf_failed)"])
    fi
  ])

  ## disable CLA in case PDF is disabled
  AS_IF([ test "x$with_cla_pdf" == "xno" ], [ with_cla=no ])
])

AC_SUBST([CLA_I_PDF_NAME], [ifelse([$with_cla],
                                      [no],
                                      [],
                                      [$CLA_I_PROJECT_NAME.pdf])])
AC_SUBST([CLA_E_PDF_NAME], [ifelse([$with_cla],
                                      [no],
                                      [],
                                      [$CLA_E_PROJECT_NAME.pdf])])


AM_CONDITIONAL(WITH_CLA, test "x$with_cla" != xno)
AM_CONDITIONAL(WITH_CLA_BUILD, test "x$cla_requirements_pdf" != xno)
AM_CONDITIONAL(WITH_CLA_PDF, test "x$with_cla_pdf" != xno)
])

