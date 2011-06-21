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
