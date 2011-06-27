# This file is part of Autoconf.                       -*- Autoconf -*-

# Copyright (C) 2011 Ronny Lorenz
#

##------------------##
## Default Options. ##
##------------------##

AC_ARG_WITH([documentation],
            [AS_HELP_STRING([--without-documentation],
              [disable installation of RNAlib reference manual])],
            [],
            [with_documentation=yes])

AC_ARG_WITH([documentation-pdf],
            [AS_HELP_STRING([--without-documentation-pdf],
              [disable installation of PDF RNAlib reference manual])],
            [],
            [with_documentation_pdf=yes])

AC_ARG_WITH([documentation-html],
            [AS_HELP_STRING([--without-documentation-html],
              [disable installation of HTML RNAlib reference manual])],
            [],
            [with_documentation_html=yes])



AC_DEFUN([AC_RNA_DOCUMENTATION_INIT],[

AC_CHECK_PROGS(doxygen, [doxygen],no)
AC_CHECK_PROGS(pdflatex,[pdflatex],no)
AC_CHECK_PROGS(latex,[latex],no)
AC_CHECK_PROGS(makeindex,[makeindex],no)
AC_CHECK_PROGS(dot,[dot],no)
AC_CHECK_PROGS(egrep,[dot],no)

])

