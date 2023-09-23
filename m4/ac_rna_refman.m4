# RNA_ENABLE_DOXYGEN_REFMAN(PROJECT_NAME, [config-file], [documentation-output-directory])
#
#

AC_DEFUN([RNA_ENABLE_DOXYGEN_REFMAN],[

RNA_ADD_PACKAGE([doc_pdf],
                [PDF RNAlib reference manual],
                [yes], [],
                [ with_doc=yes])
RNA_ADD_PACKAGE([doc_html],
                [HTML RNAlib reference manual],
                [yes], [],
                [ with_doc=yes])
RNA_ADD_PACKAGE([doc],
                [RNAlib reference manual],
                [yes],
                [ with_doc=no
                  with_doc_pdf=no
                  with_doc_html=no],
                [ with_doc_pdf=yes
                  with_doc_html=yes])

AC_REQUIRE([RNA_LATEX_ENVIRONMENT])
AC_PATH_PROG(doxygen, [doxygen],no)
AC_PATH_PROG(dot,[dot],no)
AC_PATH_PROG(grep,[grep],no)
AC_PATH_PROG(perl,[perl],no)
AC_CHECK_PROGS([SPHINXBUILD], [sphinx-build sphinx-build3], [no])

# check whether we are able to generate the doxygen documentation
RNA_PACKAGE_IF_ENABLED([doc],[
  doxygen_xml="yes"
  doc_build_requirements_html="yes"
  doc_build_requirements_pdf="yes"
  doc_build_pdf_failed=""
  doc_build_html_failed=""

  AC_SUBST([DOXYGEN_PROJECT_NAME], [$1-$PACKAGE_VERSION])
  AC_SUBST([DOXYGEN_SRCDIR], [$srcdir])
  AC_SUBST([DOXYGEN_DOCDIR], [doc/doxygen])
  AC_SUBST([DOXYGEN_CONF], [doxygen.conf])
  AC_SUBST([DOXYGEN_HAVE_DOT],[ifelse([$dot], [no], [NO], [YES])])
  AC_SUBST([SPHINX_SRCDIR], [doc/source])
  AC_SUBST([SPHINX_CMD_LATEX], ["$LATEX_CMD -interaction=nonstopmode -halt-on-error"])
  AC_SUBST([SPHINX_CMD_BIBTEX], [$BIBTEX_CMD])
  AC_SUBST([SPHINX_CMD_MAKEINDEX], [$MAKEINDEX_CMD])

  # First, check whether we are able to build doxygen XML data
  # or whether we at least have it already present
  if test "x$doxygen" != xno && test "x$perl" = xno;
  then
    AC_MSG_WARN([perl command required for automatic (re)generation of reference manual!])
    doxygen=no
  fi

  if test "x$doxygen" = xno;
  then
    AC_RNA_TEST_FILE( [$DOXYGEN_DOCDIR/xml/index.xml],
                      [doxygen_xml=yes],
                      [doxygen_xml=no
                       doc_build_requirements_html="no"
                       doc_build_requirements_pdf="no"
                       doc_build_pdf_failed="no doxygen xml"
                       doc_build_html_failed="no doxygen xml"])
  else
    doxygen_xml=yes
  fi

  AS_IF([test "x$SPHINXBUILD" = xno],
        [
          AC_MSG_WARN(sphinx-build is required to build the reference manual)
          doc_build_requirements_html="no"
          doc_build_requirements_pdf="no"
          AC_RNA_APPEND_VAR_COMMA(doc_build_pdf_failed, [no sphinx-build])
          AC_RNA_APPEND_VAR_COMMA(doc_build_html_failed, [no sphinx-build])
        ])

  RNA_PACKAGE_IF_DISABLED([python],
        [
          AC_MSG_WARN(Python module build required to build the reference manual)
          doc_build_requirements_html="no"
          doc_build_requirements_pdf="no"
          AC_RNA_APPEND_VAR_COMMA(doc_build_pdf_failed, [no python])
          AC_RNA_APPEND_VAR_COMMA(doc_build_html_failed, [no python])
        ])

  if test "x$LATEX_CMD" = xno;
  then
    AC_MSG_WARN([Could not find pdflatex!])
    AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
      doc_build_requirements_pdf="no"
      AC_RNA_APPEND_VAR_COMMA(doc_build_pdf_failed, [no pdflatex])
  else
    RNA_LATEX_TEST_PACKAGES(
      [inputenc cmap fontenc amsmath amssymb amstext babel tgtermes tgheros fncychap geometry hyperref hypcap], [], [
      doc_build_requirements_pdf="no"
      AC_RNA_APPEND_VAR_COMMA(doc_build_pdf_failed, [LaTeX packages missing])
    ])
  fi

  if test "x$makeindex" = xno;
  then
    AC_MSG_WARN([makeindex command not found on your system!])
    AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
    doc_build_requirements_pdf="no"
    AC_RNA_APPEND_VAR_COMMA(doc_build_pdf_failed, [no makeindex])
  fi

  if test "x$bibtex" = xno;
  then
    AC_MSG_WARN([bibtex command not found on your system!])
    AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
    doc_build_requirements_pdf="no"
    AC_RNA_APPEND_VAR_COMMA(doc_build_pdf_failed, [no bibtex])
  fi

  if test "x$grep" = xno;
  then
    AC_MSG_WARN([grep command not found on your system!])
    AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
    doc_build_requirements_pdf="no"
    AC_RNA_APPEND_VAR_COMMA(doc_build_pdf_failed, [no grep])
  fi

])


RNA_PACKAGE_IF_ENABLED([doc_pdf],[
  if test "x$doc_build_requirements_pdf" = xno;
  then
    # check if a generated reference manual already exists
    RNA_PACKAGE_IF_ENABLED([doc_pdf],[
      AC_RNA_TEST_FILE( [doc/$DOXYGEN_PROJECT_NAME.pdf],
                        [with_doc_pdf=yes],
                        [with_doc_pdf=no
                         doc_pdf_failed="($doc_build_pdf_failed)"])
    ])
  else
    doc_build_pdf_failed=""
  fi
])

RNA_PACKAGE_IF_ENABLED([doc_html],[
  if test "x$doc_build_requirements_html" = xno;
  then
    RNA_PACKAGE_IF_ENABLED([doc_html],[
      AC_RNA_TEST_FILE( [doc/build/html/index.html],
                        [with_doc_html=yes],
                        [with_doc_html=no
                         doc_html_failed="($doc_build_html_failed)"])
    ])
  else
    doc_build_html_failed=""
  fi
])

# setup everything in order to generate reference manual

RNA_PACKAGE_IF_ENABLED([doc],[


  AC_CONFIG_FILES([${DOXYGEN_DOCDIR}/${DOXYGEN_CONF}])
  AC_CONFIG_FILES([doc/conf.py])
  AC_CONFIG_FILES([${SPHINX_SRCDIR}/conf.py])
  AC_CONFIG_FILES([${SPHINX_SRCDIR}/install.rst])

  if test "x$with_doc_pdf" = "x$with_doc_html";
  then
    if test "x$with_doc_pdf" = "x$with_doc";
    then
        if test "x$with_doc_pdf" = "xno";
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
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_XML, test "x$doxygen_xml" != xno)
])

