testdir = $(bindir)/Test

LDADD = -L../lib/ -lcommcode
AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/lib/
AM_CFLAGS  =-I$(top_srcdir) -I$(top_srcdir)/lib/


EXTRA_DIST = 	test_HMC.xml	 \
		test_RNG.xml	 \
		test_Wilson.xml	 \
		test_ODWF.xml 	\
		test_ResMass.xml \
		test_IO.xml

test_PROGRAMS = test_HMC \
		test_RNG \
		test_Wilson \
		test_ODWF \
		test_ResMass \
		test_IO

test_HMC_SOURCES = 	test_HMC.cpp \
		  	testerHMC.cpp

test_RNG_SOURCES = test_RNG.cpp

test_Wilson_SOURCES =	test_wilson.cpp \
			testerWilson.cpp

test_ODWF_SOURCES = test_DomainWall.cpp \
			testerODWF.cpp

test_ResMass_SOURCES = test_ResidualMass.cpp \
			testerResMass.cpp

test_IO_SOURCES = test_IO.cpp 	\
		testerIO.cpp

include_HEADERS =	test_HMC.hpp \
			test_wilson.hpp \
			test_DomainWall.hpp \
			test_ResidualMass.hpp \
			test_IO.hpp