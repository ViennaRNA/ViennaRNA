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

tutorial_verbose = $(tutorial_verbose_@AM_V@)
tutorial_verbose_ = $(tutorial_verbose_@AM_DEFAULT_V@)

if WITH_TTY_COLORS

#am__v_GEN_0 = @echo -e "  \e[1;34mGEN\e[0m     " $@;
#am__v_CC_0 = @echo -e "  \e[1;34mCC\e[0m      " $@;
#am__v_CCLD_0 = @echo -e "  \e[1;34mCCLD\e[0m    " $@;
#am__v_CXX_0 = @echo -e "  \e[1;34mCXX\e[0m     " $@;
#am__v_CXXLD_0 = @echo -e "  \e[1;34mCXXLD\e[0m   " $@;
ggo_verbose_0 = @echo -e "  \e[1;34mGGO\e[0m      $@";
xxd_verbose_0 = @echo -e "  \e[1;34mXXD\e[0m      $@";
doxy2swig_verbose_0 = @echo -e "  \e[1;34mDOXY2SWIG\e[0m $@";
sphinx_verbose_0 = @echo -e "  \e[1;34mSPHINX\e[0m   $@";
doxygen_verbose_0 = @echo -e "  \e[1;34mDOXYGEN\e[0m  $@";
cla_verbose_0 = @echo -e "  \e[1;34mCLA\e[0m      $@";
swig_verbose_0 = @echo -e "  \e[1;34mSWIG\e[0m      $@";
man2rst_verbose_0 = @echo -e "  \e[1;34mMAN2RST\e[0m  $@";
manpages_verbose_0 = @echo -e "  \e[1;34mMAN\e[0m      $@";
checkmk_verbose_0 = @echo -e "  \e[1;34mCHECKMK\e[0m  $@";
prepare_verbose_0 = @echo -e "  \e[1;34mPREPARE\e[0m  $@;
pkg_verbose_0 = @echo -e "  \e[1;34mPKG-GEN\e[0m  $@;
dmg_verbose_0 = @echo -e "  \e[1;34mDMG-GEN\e[0m  $@;
tutorial_verbose_0 = @echo -e "  \e[1;34mTUT\e[0m      $@";

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
tutorial_verbose_0 = @echo "  TUT      $@";

endif
