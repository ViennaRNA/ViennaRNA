# RNA_ENABLE_DOXYGEN_REFMAN(PROJECT_NAME, [config-file], [documentation-output-directory])
#
#




AC_DEFUN([RNA_ENABLE_DOXYGEN_REFMAN],[

    AC_REQUIRE([RNA_LATEX_ENVIRONMENT])

RNA_ADD_PACKAGE([doc_pdf],
                [PDF RNAlib reference manual],
                [yes])
RNA_ADD_PACKAGE([doc_html],
                [HTML RNAlib reference manual],
                [yes])
RNA_ADD_PACKAGE([doc],
                [RNAlib reference manual],
                [yes],
                [ with_doc=no
                  with_doc_pdf=no
                  with_doc_html=no],
                [])


AC_PATH_PROG(doxygen, [doxygen],no)
AC_PATH_PROG(dot,[dot],no)
AC_PATH_PROG(egrep,[egrep],no)
AC_PATH_PROG(perl,[perl],no)

# check whether we are able to generate the doxygen documentation
RNA_PACKAGE_IF_ENABLED([doc],[
  if test "x$doxygen" != xno;
  then

    # test for programs necessary in order to use doxygen
    if test "x$LATEX_CMD" = xno;
    then
      AC_MSG_WARN([Could not find pdflatex!])
      AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
      doxygen_failed="pdflatex command is missing!"
      doxygen=no
    else
      RNA_LATEX_TEST_PACKAGES(
        [fixltx2e calc graphicx makeidx multicol multirow textcomp xcolor fontenc helvet courier amssymb sectsty geometry fancyhdr natbib tocloft amsmath amsfonts newunicodechar hyperref caption etoc ], [], [
          doxygen_failed="Required LaTeX packages are missing!"
          doxygen=no
      ])
    fi

    if test "x$makeindex" = xno;
    then
      AC_MSG_WARN([makeindex command not found on your system!])
      AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
      doxygen_failed="makeindex command is missing!"
      doxygen=no
    fi

    if test "x$bibtex" = xno;
    then
      AC_MSG_WARN([bibtex command not found on your system!])
      AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
      doxygen_failed="bibtex command is missing!"
      doxygen=no
    fi

    if test "x$egrep" = xno;
    then
      AC_MSG_WARN([egrep command not found on your system!])
      AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
      doxygen_failed="egrep command is missing!"
      doxygen=no
    fi

    if test "x$dot" = xno;
    then
      AC_MSG_WARN([dot command not found on your system!])
      AC_MSG_WARN([deactivating graph output in reference manual!])
    fi

    if test "x$perl" = xno;
    then
      AC_MSG_WARN([perl command not found on your system!])
      AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
      doxygen_failed="perl command is missing!"
      doxygen=no
    fi
  else
    doxygen_failed="doxygen command is missing!"
  fi
])


# setup everything in order to generate the doxygen configfile

RNA_PACKAGE_IF_ENABLED([doc],[

  AC_SUBST([DOXYGEN_PROJECT_NAME], [$1-$PACKAGE_VERSION])
  AC_SUBST([DOXYGEN_SRCDIR], [$srcdir])
  AC_SUBST([DOXYGEN_DOCDIR], [ifelse([$3], [], [doc], [$3])])
  AC_SUBST([DOXYGEN_CONF], [ifelse([$2], [], [doxygen.conf], [$2])])


# prepare the config file for doxygen if we are able to generate a reference manual
  if test "x$doxygen" != xno;
  then

    AC_SUBST([DOXYGEN_CMD_LATEX], ["$LATEX_CMD -interaction=nonstopmode -halt-on-error"])
    AC_SUBST([DOXYGEN_CMD_BIBTEX], [$BIBTEX_CMD])
    AC_SUBST([DOXYGEN_CMD_MAKEINDEX], [$MAKEINDEX_CMD])
    AC_SUBST([DOXYGEN_HAVE_DOT],[ifelse([$dot], [no], [NO], [YES])])
    AC_SUBST([DOXYGEN_WITH_PDFLATEX], [YES])
    AC_SUBST([DOXYGEN_GENERATE_HTML], [ifelse([$with_doc_html], [no], [NO], [YES])])
    AC_SUBST([DOXYGEN_GENERATE_LATEX], [ifelse([$with_doc_pdf], [no], [NO], [YES])])

    AC_CONFIG_FILES([${DOXYGEN_DOCDIR}/${DOXYGEN_CONF}])
    AC_CONFIG_FILES([${DOXYGEN_DOCDIR}/refman.include/install.dox])

  else

# otherwise check if a generated reference manual already exists

    RNA_PACKAGE_IF_ENABLED([doc_pdf],[
      AC_RNA_TEST_FILE( [$DOXYGEN_DOCDIR/$DOXYGEN_PROJECT_NAME.pdf],
                        [with_doc_pdf=yes],
                        [with_doc_pdf=no
                         doc_pdf_failed="($doxygen_failed)"])])

    RNA_PACKAGE_IF_ENABLED([doc_html],[
      AC_RNA_TEST_FILE( [$DOXYGEN_DOCDIR/html/index.html],
                        [with_doc_html=yes],
                        [with_doc_html=no
                         doc_html_failed="($doxygen_failed)"])])

    if test "x$with_doc_pdf" = "x$with_doc_html";
    then
      if test "x$with_doc_pdf" = xno;
      then
        with_doc=no
      fi
    fi
  fi
])

AC_SUBST([REFERENCE_MANUAL_PDF_NAME], [ifelse([$with_doc_pdf],
                                              [no],
                                              [],
                                              [$DOXYGEN_PROJECT_NAME.pdf])])
AC_SUBST([REFERENCE_MANUAL_TAGFILE],  [ifelse([$doxygen],
                                              [no],
                                              [],
                                              [$DOXYGEN_PROJECT_NAME.tag])])


# setup variables used in Makefile.am

# Define ${htmldir} if the configure script was created with a version of
# autoconf older than 2.60
# Alternatively, if ${htmldir} is exactly '${docdir}', append a /html to
# separate html files from rest of doc.
# Otherwise, just append the PACKAGE_NAME to the htmldir
if test "x${htmldir}" = "x";
then
  AC_MSG_WARN([resetting htmldir])
  htmldir="${docdir}/html"
fi

if test "x${htmldir}" = 'x${docdir}';
then
  htmldir="${docdir}/html"
else
  htmldir=${htmldir}/${PACKAGE_NAME}
fi

AC_SUBST(htmldir, [${htmldir}])

#

AM_CONDITIONAL(WITH_REFERENCE_MANUAL, test "x$with_doc" != xno)
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_BUILD, test "x$doxygen" != xno)
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_PDF, test "x$with_doc_pdf" != xno)
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_HTML, test "x$with_doc_html" != xno)
])

