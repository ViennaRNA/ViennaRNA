if HAVE_HELP2MAN
  man1_MANS = $(ax_help2man_MANS)
  @AX_HELP2MAN_RULES@
endif

# remove all generated man pages
maintainer-clean-local:
	-rm -f $(man1_MANS)

distclean-local:
	-rm -f $(man1_MANS)

