include $(top_srcdir)/silent_rules.mk

AM_CFLAGS = $(RNA_CFLAGS) $(PTHREAD_CFLAGS)
AM_CXXFLAGS = $(RNA_CXXFLAGS) $(PTHREAD_CFLAGS)
AM_CPPFLAGS = $(RNA_CPPFLAGS) ${SVM_INC} -I$(top_srcdir)/src ${JSON_INC}
AM_LDFLAGS = $(RNA_LDFLAGS) $(PTHREAD_LIBS)

vrna_naviewdir         = $(pkgincludedir)/plotting/naview

# convenience library
noinst_LTLIBRARIES =  libRNA_naview.la

vrna_naview_HEADERS = naview.h

libRNA_naview_la_SOURCES = naview.c
