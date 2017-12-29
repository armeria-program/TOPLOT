#===============================================================================
# Makefile.am : for automake
# Copyright (C) 2006-2017 Jens Kleinjung
#===============================================================================

bin_PROGRAMS = toplot

AM_CPPFLAGS = $(INTI_CFLAGS)
AM_CFLAGS = -Wall

toplot_SOURCES = \
arg.c arg.h alphabet.c alphabet.h config.h contact.c contact.h \
getpdbs.c getpdbs.h structure.c structure.h toplot.c toplot.h

toplot_LDADD = $(INTI_LIBS)

EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

