
EXTRA_DIST = resources/welcome.txt

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
MACOSX_DMG_NAME=@PACKAGE_NAME@-@PACKAGE_VERSION@
MACOSX_DMG_FILE=${MACOSX_DMG_NAME}-macosx.dmg

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

##
## copy README and COPYING from top_srcdir to this ./resources/
##
prepare-pkg:
	$(MKDIR_P) ./resources
	cp $(top_srcdir)/README $(top_srcdir)/COPYING ./resources/ 

##
## prepare Core staging directory
##
prepare-staging-core:
	$(remove_staging_core_dir)
	$(MKDIR_P) $(staging_core_dir)$(prefix)/bin
	cd $(top_srcdir)/src/bin \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_core_dir)
	cd $(top_srcdir)/src/Utils \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_core_dir)
	cd $(top_srcdir)/man \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_core_dir)
	cd $(top_srcdir)/misc \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_core_dir)

##
## prepare Cluster Analysis staging directory
##
prepare-staging-cluster:
	$(remove_staging_cluster_dir)
	cd $(top_srcdir)/src/Cluster \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_cluster_dir)

##
## prepare Kinfold staging directory
##
prepare-staging-kinfold:
	$(remove_staging_kinfold_dir)
	cd $(top_srcdir)/src/Kinfold \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_kinfold_dir)

##
## prepare RNAforester staging directory
##
prepare-staging-forester:
	$(remove_staging_forester_dir)
	$(MKDIR_P) $(staging_forester_dir)$(prefix)/bin
	cp $(top_srcdir)/src/RNAforester/src/RNAforester $(staging_forester_dir)$(prefix)/bin
	cd $(top_srcdir)/src/RNAforester/man \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_forester_dir)

##
## prepare Kinwalker staging directory
##
prepare-staging-kinwalker:
	$(remove_staging_kinwalker_dir)
	cd $(top_srcdir)/src/Kinwalker \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_kinwalker_dir)

##
## prepare Development staging directory
##
prepare-staging-libRNA:
	$(remove_staging_libRNA_dir)
	cd $(top_srcdir)/doc \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_libRNA_dir)
	cd $(top_srcdir)/src/ViennaRNA \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_libRNA_dir)
	$(MKDIR_P) $(staging_core_dir)$(prefix)/lib/pkgconfig
	cp $(top_srcdir)/RNAlib2.pc $(staging_libRNA_dir)$(prefix)/lib/pkgconfig

##
## prepare Perl staging directory
##
prepare-staging-perl:
	$(remove_staging_perl_dir)
	cd $(top_srcdir)/interfaces/Perl \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_perl_dir)

##
## prepare Python staging directory
##
prepare-staging-python:
	$(remove_staging_python_dir)
	cd $(top_srcdir)/interfaces/Python \
	&& $(MAKE) install DESTDIR=$(abs_path_thisdir)/$(staging_python_dir)

"${CORE_PKG_FILE}": prepare-staging-core
	$(AM_V_GEN)pkgbuild --identifier org.TBI.ViennaRNA.core --version @PACKAGE_VERSION@ \
      --root ${staging_core_dir} "${CORE_PKG_FILE}"

"${DEV_PKG_FILE}": prepare-staging-libRNA
	$(AM_V_GEN)pkgbuild --identifier org.TBI.ViennaRNA.libRNA --version @PACKAGE_VERSION@ \
      --root ${staging_libRNA_dir} "${DEV_PKG_FILE}"

"${CLUSTER_PKG_FILE}": prepare-staging-cluster
	$(AM_V_GEN)pkgbuild --identifier org.TBI.ViennaRNA.cluster --version @PACKAGE_VERSION@ \
      --root ${staging_cluster_dir} "${CLUSTER_PKG_FILE}"

"${KINFOLD_PKG_FILE}": prepare-staging-kinfold
	$(AM_V_GEN)pkgbuild --identifier org.TBI.ViennaRNA.kinfold --version @PACKAGE_VERSION@ \
      --root ${staging_kinfold_dir} "${KINFOLD_PKG_FILE}"

"${KINWALKER_PKG_FILE}": prepare-staging-kinwalker
	$(AM_V_GEN)pkgbuild --identifier org.TBI.ViennaRNA.kinwalker --version @PACKAGE_VERSION@ \
      --root ${staging_kinwalker_dir} "${KINWALKER_PKG_FILE}"

"${FORESTER_PKG_FILE}": prepare-staging-forester
	$(AM_V_GEN)pkgbuild --identifier org.TBI.ViennaRNA.rnaforester --version @PACKAGE_VERSION@ \
      --root ${staging_forester_dir} "${FORESTER_PKG_FILE}"

"${PERL_PKG_FILE}": prepare-staging-perl
	$(AM_V_GEN)pkgbuild --identifier org.TBI.ViennaRNA.perl --version @PACKAGE_VERSION@ \
      --root ${staging_perl_dir} "${PERL_PKG_FILE}"

"${PYTHON_PKG_FILE}": prepare-staging-python
	$(AM_V_GEN)pkgbuild --identifier org.TBI.ViennaRNA.python --version @PACKAGE_VERSION@ \
      --root ${staging_python_dir} "${PYTHON_PKG_FILE}"

Installer.pkg:  prepare-pkg \
                "${CORE_PKG_FILE}" \
                "${DEV_PKG_FILE}" \
                "${CLUSTER_PKG_FILE}" \
                "${KINFOLD_PKG_FILE}" \
                "${FORESTER_PKG_FILE}" \
                "${KINWALKER_PKG_FILE}" \
                "${PERL_PKG_FILE}" \
                "${PYTHON_PKG_FILE}"
	$(AM_V_GEN)productbuild  --distribution ./Distribution.xml \
      --package-path . \
      --resources ./resources \
      ./Installer.pkg

"${DMG_FILE}": Installer.pkg
	$(MKDIR_P) disk_image
	cp Installer.pkg disk_image
	$(AM_V_GEN)hdiutil create -ov -fs HFSX -layout SPUD -format UDZO \
    -volname ${MACOSX_DMG_NAME} -srcfolder disk_image \
    -imagekey zlib-level=9 -scrub "${MACOSX_DMG_FILE}"

PHONY:  prepare-pkg \
        prepare-staging-core \
        prepare-staging-libRNA \
        prepare-staging-cluster \
        prepare-staging-kinfold \
        prepare-staging-forester \
        prepare-staging-kinwalker \
        prepare-staging-perl \
        prepare-staging-python
