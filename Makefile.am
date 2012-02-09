SUBDIRS = lib tests docs

EXTRA_DIST = 	m4	\
		README 	\
		todo

if DOXYGEN_DOC
directory = $(top_srcdir)/docs/

doxyfile:
	(cd $(directory) && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1; 

endif 