include ../generic.mk

SOURCES = setup.py
EXTRA_DIST = \
  $(SWIG_wrapper) \
  setup.py \
  version_test.py \
  __init__.py

AM_CPPFLAGS = -I$(top_srcdir)/src

if WITH_LARGE_PF
  USE_LARGE_PF = -DLARGE_PF
else
  USE_LARGE_PF =
endif


##
## SWIG interface wrapper code generation
##

$(SWIG_wrapper) $(PERL_MODULE): $(SWIG_src)
	$(AM_V_GEN)swig $(INCLUDES) $(USE_LARGE_PF) \
                  -python -shadow \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)


##
## Python interface building
##

Python_interface: $(SWIG_wrapper) setup.py
	$(AM_V_GEN)@PythonCmd@ setup.py build

##
## Misc make rules
##

all-local: Python_interface
	

install-exec-hook:
	@PythonCmd@ setup.py install \
  `test -n "$(DESTDIR)" && echo --root $(DESTDIR)` \
  --record install.record

uninstall-hook:
	cat install.record \
  | awk '{print "$(DESTDIR)"$$1}' \
  | xargs -d\\n rm -vf --


distclean-local:
	-r -f install.record
  
mostlyclean-local: 
	@PythonCmd@ setup.py clean --all

# remove all generated content
maintainer-clean-local:
	-rm -f __init__.py RNA_wrap.c install.record

