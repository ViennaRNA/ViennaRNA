
#
# Kinfold subpackage
#

AC_DEFUN([RNA_ENABLE_PKG_KINFOLD],[
  RNA_ADD_PACKAGE([kinfold],
                  [Kinfold program],
                  [yes],
                  [with_kinfold=no],
                  [with_kinfold=yes],
                  [${srcdir}/src/Kinfold/Makefile.am])

  RNA_PACKAGE_IF_ENABLED([kinfold],[
    AC_CONFIG_SUBDIRS([src/Kinfold])
  ])

  AM_CONDITIONAL(MAKE_KINFOLD, test "$with_kinfold" != "no")
])


#
# RNAforester subpackage
#
AC_DEFUN([RNA_ENABLE_PKG_FORESTER],[
  RNA_ADD_PACKAGE([forester],
                  [RNAforester program],
                  [yes],
                  [with_forester=no],
                  [with_forester=yes],
                  [${srcdir}/src/RNAforester/Makefile.am])

  RNA_PACKAGE_IF_ENABLED([forester],[
    AC_CONFIG_SUBDIRS([src/RNAforester])
  ])

  AM_CONDITIONAL(MAKE_FORESTER, test "$with_forester" != "no")
])


#
#
#

AC_DEFUN([RNA_ENABLE_PKG_CLUSTER],[
  RNA_ADD_PACKAGE([cluster],
                  [AnalyseSeqs and AnalyseDists],
                  [no],
                  [with_cluster=yes],
                  [with_cluster=no],
                  [${srcdir}/src/Cluster/Makefile.am])

  RNA_PACKAGE_IF_ENABLED([cluster],[
    AC_DEFINE([WITH_CLUSTER], [1], [Analyse{Dists,Seqs}])
    AC_SUBST([CLUSTER_DIR], [Cluster])
    AC_CONFIG_FILES([src/Cluster/Makefile])
  ])

  AM_CONDITIONAL(MAKE_CLUSTER, test "$with_cluster" = "yes")
])
