EXTRA_DIST =	doxy.cfg.in 	\
		doxy.cfg.jp	\
		keklogo-c.jpg 	\
		Makefile.am	\
		footer.html	\
		header.html	\
		IroIroCodeBanner250.png \
		JICFUSsymbolmark170px.png
		customdoxygen.css


doxyfile:
	echo Timestamp > doxyfile.stamp
	$(DOXYGEN) doxy.cfg
	cp *.png ./doxy*/html/


