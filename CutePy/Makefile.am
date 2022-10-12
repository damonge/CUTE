EXTRA_DIST = README.md

lib_LTLIBRARIES = libcute.la
libcute_la_SOURCES = src/utils.c src/cute_bins.c src/cute_angular.c
libcute_la_CPPFLAGS = $(OPENMP_CFLAGS) -I./src/
libcute_la_LDFLAGS = $(OPENMP_CFLAGS) -lm

include_HEADERS = src/cute.h
