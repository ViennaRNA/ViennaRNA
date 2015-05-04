if WITH_LARGE_PF
  USE_LARGE_PF = -DLARGE_PF
else
  USE_LARGE_PF =
endif

SOURCES = setup.py

EXTRA_DIST = RNA_wrap.c setup.py version_test.py __init__.py

RNA_wrap.c: ../RNA.i ../tmaps.i
	swig -python -shadow $(USE_LARGE_PF) -I../ RNA.i

all-local: RNA_wrap.c setup.py
	@PythonCmd@ setup.py build

install-exec-local: all-local
	@PythonCmd@ setup.py install --prefix=$(prefix) --record install.record

#check-local: all-local
#	$(MAKE) -f Makefile.perl test

mostlyclean-local: 
	@PythonCmd@ setup.py clean --all

# remove all generated content
maintainer-clean-local:
	-rm -f __init__.py RNA_wrap.c

