AUTOMAKE_OPTIONS = subdir-objects

xxd_verbose = $(xxd_verbose_@AM_V@)
xxd_verbose_ = $(xxd_verbose_@AM_DEFAULT_V@)
xxd_verbose_0 = @echo "  XXD      $@";

SUFFIXES = .par .ps .hex

POSTSCRIPT_FILES = $(PS_TEMPLATE_FILES)
ENERGY_PARAMETER_FILES = $(PARAMETER_FILES)

vpath %.par = ${top_srcdir}

##  Note that we only (re-)build .hex files for direct inclusion into our
##  library if the user has the xxd tool installed. Otherwise, any changes
##  to the postscript templates will NOT propagate into the binary files

if VRNA_AM_SWITCH_HAS_XXD

## The rule for generating .hex files from .ps
.ps.hex:
	$(xxd_verbose)$(XXD) -i < $< > $@ && echo ',0x00' >> $@

templates_postscript.h: $(PS_TEMPLATE_FILES_HEX)

.par.hex:
	$(xxd_verbose)$(MKDIR_P) misc && $(XXD) -i < $< > $@ && echo ',0x00' >> $@

energy_parameter_sets.h: $(PARAMETER_FILES_HEX)

endif


BUILT_SOURCES = \
    $(PS_TEMPLATE_FILES_HEX) \
    $(PARAMETER_FILES_HEX)


EXTRA_DIST = \
    ps_templates.txt \
    $(POSTSCRIPT_FILES) \
    templates_postscript.h \
    energy_parameter_sets.h \
    ${BUILT_SOURCES}


maintainer-clean-local:
	rm -rf $(BUILT_SOURCES)
