SUBDIRS = lib tests docs

EXTRA_DIST = README 		\
	docs/doxy.cfg.in 	\
	docs/doxy.cfg.jp	\
	docs/keklogo-c.jpg 	\
	docs/Makefile.am

if DOXYGEN_DOC
directory = $(top_srcdir)/docs/

doxyfile:
	(cd $(directory) && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1; 

endif 