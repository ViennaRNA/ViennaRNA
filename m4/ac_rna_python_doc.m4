AC_DEFUN([RNA_ENABLE_PYTHON_DOC],[

RNA_ADD_PACKAGE([pydoc],
                [python documentation],
                [no],
                [],
                [])

AC_CHECK_PROGS([SPHINXBUILD], [sphinx-build sphinx-build3], [no])

# check whether we are able to generate the doxygen documentation
RNA_PACKAGE_IF_ENABLED([pydoc],[
  AS_IF([test "x$SPHINXBUILD" = xno],
        [
          AC_MSG_WARN(sphinx-build is required to build python documentation)
          pydoc_enabled_but_failed="(sphinx-build unavailable)"
          with_pydoc=no
        ])
  RNA_PACKAGE_IF_DISABLED([python],
        [
          AC_MSG_WARN(Python module build required to build python documentation)
          #pydoc_enabled_but_failed="Python module not to be build"
          with_pydoc=no
        ])
  RNA_PACKAGE_IF_DISABLED([doc],
        [
          AC_MSG_WARN(Doxygen documentation is required to build python documentation)
          pydoc_enabled_but_failed="(doxygen documentation unavailable)"
          with_pydoc=no
        ])
])

# setup variables used in Makefile.am
pyhtmldir=${htmldir}/python

AC_SUBST(pyhtmldir, [${pyhtmldir}])

AM_CONDITIONAL(WITH_PYTHON_DOC, test "x$with_pydoc" != xno)

])
