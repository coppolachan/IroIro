SUBDIRS = lib docs

EXTRA_DIST = README 		\
	docs/doxy.cfg.in 	\
	docs/doxy.cfg.jap	\
	docs/keklogo-c.jpg 	\
	docs/Makefile.am


if HAVE_DOXYGEN
directory = $(top_srcdir)/docs/

doxyfile:
	(cd $(directory) && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1; 

endif 