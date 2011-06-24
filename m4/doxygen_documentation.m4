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

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(ON)
DX_INIT_DOXYGEN(RNAlib, doxygen.conf, doc)
