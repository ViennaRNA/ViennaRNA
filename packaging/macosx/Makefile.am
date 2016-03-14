prepare_verbose = $(prepare_verbose_@AM_V@)
prepare_verbose_ = $(prepare_verbose_@AM_DEFAULT_V@)
prepare_verbose_0 = @echo PREPARE $@;
pkg_verbose = $(pkg_verbose_@AM_V@)
pkg_verbose_ = $(pkg_verbose_@AM_DEFAULT_V@)
pkg_verbose_0 = @echo PKG-GEN $@;
dmg_verbose = $(dmg_verbose_@AM_V@)
dmg_verbose_ = $(dmg_verbose_@AM_DEFAULT_V@)
dmg_verbose_0 = @echo DMG-GEN $@;

EXTRA_DIST =  resources/welcome.txt \
              resources/conclusion.html \
              resources/ohm.png

srcdir  = @srcdir@
thisdir = packaging/macosx
abs_path_thisdir = $(abs_top_srcdir)/$(thisdir)

staging_core_dir = staging_core
staging_libRNA_dir = staging_libRNA
staging_perl_dir = staging_perl
staging_python_dir = staging_python
staging_kinfold_dir = staging_kinfold
staging_kinwalker_dir = staging_kinwalker
staging_forester_dir = staging_forester
staging_cluster_dir = staging_cluster

CORE_PKG_FILE="@PACKAGE_NAME@\ Package\ @PACKAGE_VERSION@\ Core".pkg
CLUSTER_PKG_FILE="@PACKAGE_NAME@\ Package\ @PACKAGE_VERSION@\ Cluster\ Analysis".pkg
KINFOLD_PKG_FILE="@PACKAGE_NAME@\ Package\ @PACKAGE_VERSION@\ Kinfold".pkg
KINWALKER_PKG_FILE="@PACKAGE_NAME@\ Package\ @PACKAGE_VERSION@\ Kinwalker".pkg
FORESTER_PKG_FILE="@PACKAGE_NAME@\ Package\ @PACKAGE_VERSION@\ RNAforester".pkg
DEV_PKG_FILE="@PACKAGE_NAME@\ Package\ @PACKAGE_VERSION@\ Development".pkg
PERL_PKG_FILE="@PACKAGE_NAME@\ Package\ @PACKAGE_VERSION@\ Perl".pkg
PYTHON_PKG_FILE="@PACKAGE_NAME@\ Package\ @PACKAGE_VERSION@\ Python".pkg
VRNA_PKG_FILE="@PACKAGE_NAME@\ Package\ @PACKAGE_VERSION@\ Installer".pkg
MACOSX_DMG_NAME="@PACKAGE_NAME@ @PACKAGE_VERSION@"
MACOSX_DMG_FILE="@PACKAGE_NAME@-@PACKAGE_VERSION@-MacOSX".dmg

remove_staging_core_dir = \
  { test ! -d $(staging_core_dir) \
    || { find $(staging_core_dir) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(staging_core_dir); }; }

remove_staging_libRNA_dir = \
  { test ! -d $(staging_libRNA_dir) \
    || { find $(staging_libRNA_dir) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(staging_libRNA_dir); }; }

remove_staging_kinfold_dir = \
  { test ! -d $(staging_kinfold_dir) \
    || { find $(staging_kinfold_dir) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(staging_kinfold_dir); }; }

remove_staging_cluster_dir = \
  { test ! -d $(staging_cluster_dir) \
    || { find $(staging_cluster_dir) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(staging_cluster_dir); }; }

remove_staging_kinwalker_dir = \
  { test ! -d $(staging_kinwalker_dir) \
    || { find $(staging_kinwalker_dir) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(staging_kinwalker_dir); }; }

remove_staging_forester_dir = \
  { test ! -d $(staging_forester_dir) \
    || { find $(staging_forester_dir) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(staging_forester_dir); }; }

remove_staging_perl_dir = \
  { test ! -d $(staging_perl_dir) \
    || { find $(staging_perl_dir) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(staging_perl_dir); }; }

remove_staging_python_dir = \
  { test ! -d $(staging_python_dir) \
    || { find $(staging_python_dir) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(staging_python_dir); }; }

if WITH_MACOSX_INSTALLER
DISK_IMAGE = "${MACOSX_DMG_FILE}"
else
DISK_IMAGE = 
endif

all: $(DISK_IMAGE)

##
## copy README and COPYING from top_srcdir to this ./resources/
##
resources:
	$(prepare_verbose)$(MKDIR_P) ./resources
	$(prepare_verbose)cp $(top_srcdir)/README $(top_srcdir)/COPYING ./resources/ 

##
## prepare Core staging directory
##
staging-core:
	$(prepare_verbose)$(remove_staging_core_dir)
	$(AM_V_at)$(MKDIR_P) $(staging_core_dir)$(prefix)/bin
	$(AM_V_at)cp uninstall-ViennaRNA-Package.sh $(staging_core_dir)$(prefix)/bin/
	$(AM_V_at)cd $(top_srcdir)/src/bin \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_core_dir)
	$(AM_V_at)cd $(top_srcdir)/src/Utils \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_core_dir)
	$(AM_V_at)cd $(top_srcdir)/man \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_core_dir)
	$(AM_V_at)cd $(top_srcdir)/misc \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_core_dir)

##
## prepare Cluster Analysis staging directory
##
staging-cluster:
	$(prepare_verbose)$(remove_staging_cluster_dir)
	$(AM_V_at)cd $(top_srcdir)/src/Cluster \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_cluster_dir)

##
## prepare Kinfold staging directory
##
staging-kinfold:
	$(prepare_verbose)$(remove_staging_kinfold_dir)
	$(AM_V_at)cd $(top_srcdir)/src/Kinfold \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_kinfold_dir)

##
## prepare RNAforester staging directory
##
staging-forester:
	$(prepare_verbose)$(remove_staging_forester_dir)
	$(AM_V_at)$(MKDIR_P) $(staging_forester_dir)$(prefix)/bin
	$(AM_V_at)cp $(top_srcdir)/src/RNAforester/src/RNAforester $(staging_forester_dir)$(prefix)/bin
	$(AM_V_at)cd $(top_srcdir)/src/RNAforester/man \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_forester_dir)

##
## prepare Kinwalker staging directory
##
staging-kinwalker:
	$(prepare_verbose)$(remove_staging_kinwalker_dir)
	$(AM_V_at)cd $(top_srcdir)/src/Kinwalker \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_kinwalker_dir)

##
## prepare Development staging directory
##
staging-libRNA:
	$(prepare_verbose)$(remove_staging_libRNA_dir)
	$(AM_V_at)cd $(top_srcdir)/doc \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_libRNA_dir)
	$(AM_V_at)cd $(top_srcdir)/src/ViennaRNA \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_libRNA_dir)
	$(AM_V_at)$(MKDIR_P) $(staging_libRNA_dir)$(prefix)/lib/pkgconfig
	$(AM_V_at)cp $(top_srcdir)/RNAlib2.pc $(staging_libRNA_dir)$(prefix)/lib/pkgconfig

##
## prepare Perl staging directory
##
staging-perl:
	$(prepare_verbose)$(remove_staging_perl_dir)
	$(AM_V_at)cd $(top_srcdir)/interfaces/Perl \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_perl_dir)

##
## prepare Python staging directory
##
staging-python:
	$(prepare_verbose)$(remove_staging_python_dir)
	$(AM_V_at)cd $(top_srcdir)/interfaces/Python \
      && $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_python_dir)

##
## Actual installer package targets below
##
"${CORE_PKG_FILE}": staging-core
	$(pkg_verbose)$(PKGBUILD) --identifier org.TBI.ViennaRNA.core --version @PACKAGE_VERSION@ \
      --root ${staging_core_dir} "${CORE_PKG_FILE}"

"${DEV_PKG_FILE}": staging-libRNA
	$(pkg_verbose)$(PKGBUILD) --identifier org.TBI.ViennaRNA.libRNA --version @PACKAGE_VERSION@ \
      --root ${staging_libRNA_dir} "${DEV_PKG_FILE}"

"${CLUSTER_PKG_FILE}": staging-cluster
	$(pkg_verbose)$(PKGBUILD) --identifier org.TBI.ViennaRNA.cluster --version @PACKAGE_VERSION@ \
      --root ${staging_cluster_dir} "${CLUSTER_PKG_FILE}"

"${KINFOLD_PKG_FILE}": staging-kinfold
	$(pkg_verbose)$(PKGBUILD) --identifier org.TBI.ViennaRNA.kinfold --version @PACKAGE_VERSION@ \
      --root ${staging_kinfold_dir} "${KINFOLD_PKG_FILE}"

"${KINWALKER_PKG_FILE}": staging-kinwalker
	$(pkg_verbose)$(PKGBUILD) --identifier org.TBI.ViennaRNA.kinwalker --version @PACKAGE_VERSION@ \
      --root ${staging_kinwalker_dir} "${KINWALKER_PKG_FILE}"

"${FORESTER_PKG_FILE}": staging-forester
	$(pkg_verbose)$(PKGBUILD) --identifier org.TBI.ViennaRNA.rnaforester --version @PACKAGE_VERSION@ \
      --root ${staging_forester_dir} "${FORESTER_PKG_FILE}"

"${PERL_PKG_FILE}": staging-perl
	$(pkg_verbose)$(PKGBUILD) --identifier org.TBI.ViennaRNA.perl --version @PACKAGE_VERSION@ \
      --root ${staging_perl_dir} "${PERL_PKG_FILE}"

"${PYTHON_PKG_FILE}": staging-python
	$(pkg_verbose)$(PKGBUILD) --identifier org.TBI.ViennaRNA.python --version @PACKAGE_VERSION@ \
      --root ${staging_python_dir} "${PYTHON_PKG_FILE}"

"${VRNA_PKG_FILE}":  resources \
                "${CORE_PKG_FILE}" \
                "${DEV_PKG_FILE}" \
                "${CLUSTER_PKG_FILE}" \
                "${KINFOLD_PKG_FILE}" \
                "${FORESTER_PKG_FILE}" \
                "${KINWALKER_PKG_FILE}" \
                "${PERL_PKG_FILE}" \
                "${PYTHON_PKG_FILE}"
	$(pkg_verbose)$(PRODUCTBUILD)  --distribution ./Distribution.xml \
      --package-path . \
      --resources ./resources \
      ./"${VRNA_PKG_FILE}"

disk_image:
	$(prepare_verbose)$(MKDIR_P) disk_image
	$(AM_V_at)cp "${VRNA_PKG_FILE}" disk_image
	$(AM_V_at)cp $(top_srcdir)/README $(top_srcdir)/NEWS $(top_srcdir)/ChangeLog $(top_srcdir)/AUTHORS $(top_srcdir)/THANKS disk_image/ 
	$(AM_V_at)cp $(top_srcdir)/COPYING disk_image/License

"${MACOSX_DMG_FILE}": "${VRNA_PKG_FILE}" disk_image
	$(dmg_verbose)$(HDIUTIL) create -ov -fs HFSX -layout SPUD -format UDZO \
    -volname ${MACOSX_DMG_NAME} -srcfolder disk_image \
    -imagekey zlib-level=9 -scrub "${MACOSX_DMG_FILE}"


clean-local:
	$(remove_staging_core_dir)
	$(remove_staging_libRNA_dir)
	$(remove_staging_kinfold_dir)
	$(remove_staging_cluster_dir)
	$(remove_staging_kinwalker_dir)
	$(remove_staging_forester_dir)
	$(remove_staging_perl_dir)
	$(remove_staging_python_dir)
	-rm -rf disk_image
	-rm -f  "${VRNA_PKG_FILE}" \
          "${PYTHON_PKG_FILE}" \
          "${PERL_PKG_FILE}" \
          "${FORESTER_PKG_FILE}" \
          "${KINWALKER_PKG_FILE}" \
          "${KINFOLD_PKG_FILE}" \
          "${CLUSTER_PKG_FILE}" \
          "${DEV_PKG_FILE}" \
          "${CORE_PKG_FILE}"


PHONY:  resources \
        staging-core \
        staging-libRNA \
        staging-cluster \
        staging-kinfold \
        staging-forester \
        staging-kinwalker \
        staging-perl \
        staging-python \
        disk_image
