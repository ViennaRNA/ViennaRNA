ggo_verbose = $(ggo_verbose_@AM_V@)
ggo_verbose_ = $(ggo_verbose_@AM_DEFAULT_V@)

xxd_verbose = $(xxd_verbose_@AM_V@)
xxd_verbose_ = $(xxd_verbose_@AM_DEFAULT_V@)

doxy2swig_verbose = $(doxy2swig_verbose_@AM_V@)
doxy2swig_verbose_ = $(doxy2swig_verbose_@AM_DEFAULT_V@)

sphinx_verbose = $(sphinx_verbose_@AM_V@)
sphinx_verbose_ = $(sphinx_verbose_@AM_DEFAULT_V@)

doxygen_verbose = $(doxygen_verbose_@AM_V@)
doxygen_verbose_ = $(doxygen_verbose_@AM_DEFAULT_V@)

cla_verbose = $(cla_verbose_@AM_V@)
cla_verbose_ = $(cla_verbose_@AM_DEFAULT_V@)

swig_verbose = $(swig_verbose_@AM_V@)
swig_verbose_ = $(swig_verbose_@AM_DEFAULT_V@)

man2rst_verbose = $(man2rst_verbose_@AM_V@)
man2rst_verbose_ = $(man2rst_verbose_@AM_DEFAULT_V@)

manpages_verbose = $(manpages_verbose_@AM_V@)
manpages_verbose_ = $(manpages_verbose_@AM_DEFAULT_V@)

checkmk_verbose   = $(checkmk_verbose_@AM_V@)
checkmk_verbose_  = $(checkmk_verbose_@AM_DEFAULT_V@)

prepare_verbose = $(prepare_verbose_@AM_V@)
prepare_verbose_ = $(prepare_verbose_@AM_DEFAULT_V@)

pkg_verbose = $(pkg_verbose_@AM_V@)
pkg_verbose_ = $(pkg_verbose_@AM_DEFAULT_V@)

dmg_verbose = $(dmg_verbose_@AM_V@)
dmg_verbose_ = $(dmg_verbose_@AM_DEFAULT_V@)

if WITH_COLOR_MAKE

ggo_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mGGO\033[0m      %s\n" $@ || echo "  GGO      $@";
xxd_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mXXD\033[0m      %s\n" $@ || echo "  XXD      $@";
doxy2swig_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mDOX2SWIG\033[0m %s\n" $@ || echo "  DOX2SWIG $@";
sphinx_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mSPHINX\033[0m   %s\n" $@ || echo -e "  SPHINX    $@";
doxygen_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mDOXYGEN\033[0m  %s\n" $@ || echo "  DOXYGEN  $@";
cla_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mCLA\033[0m      %s\n" $@ || echo "  CLA      $@";
swig_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mSWIG\033[0m     %s\n" $@ || echo "  SWIG     $@";
man2rst_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mMAN2RST\033[0m  %s\n" $@ || echo "  MAN2RST  $@";
manpages_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mMAN\033[0m      %s\n" $@ || echo "  MAN      $@";
checkmk_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mCHECKMK\033[0m  %s\n" $@ || echo "  CHECKMK  $@";
prepare_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mPREPARE\033[0m  $@ || echo "  PREPARE  $@;
pkg_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mPKG-GEN\033[0m  $@ || echo "  PKG-GEN  $@;
dmg_verbose_0 = @[ -t 1 ] && printf "  \033[1;34mDMG-GEN\033[0m  $@ || echo "  DMG-GEN  $@;

else

ggo_verbose_0 = @echo "  GGO      $@";
xxd_verbose_0 = @echo "  XXD      $@";
doxy2swig_verbose_0 = @echo "  DOX2SWIG $@";
sphinx_verbose_0 = @echo -e "  SPHINX    $@";
doxygen_verbose_0 = @echo "  DOXYGEN  $@";
cla_verbose_0 = @echo "  CLA      $@";
swig_verbose_0 = @echo "  SWIG     $@";
man2rst_verbose_0 = @echo "  MAN2RST  $@";
manpages_verbose_0 = @echo "  MAN      $@";
checkmk_verbose_0 = @echo "  CHECKMK  $@";
prepare_verbose_0 = @echo "  PREPARE  $@;
pkg_verbose_0 = @echo "  PKG-GEN  $@;
dmg_verbose_0 = @echo "  DMG-GEN  $@;

endif
