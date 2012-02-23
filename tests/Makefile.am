# Time-stamp: <2012-02-23 10:15:19 neo>

testdir = $(bindir)/Test

LDADD = -L../lib/ -lcommcode
AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/lib/
AM_CFLAGS  =-I$(top_srcdir) -I$(top_srcdir)/lib/

EXTRA_DIST = 	testResults.info	\
		test_HMC.xml	 	\
		test_RNG.xml	 	\
		test_Wilson.xml	 	\
		test_Wilson_EvenOdd.xml	\
		test_DomainWall.xml 	\
		test_ResMass.xml 	\
		test_IO.xml		\
		test_Gauge.xml		\
		test_Overlap.xml	\
		test_HMC_EvenOdd.xml	\
		test_HMC_DomainWall.xml	\
		test_MultiShiftSolver.xml

test_PROGRAMS = test_HMC_nofact 	\
		test_Gauge		\
		test_HMC_smear_noF	\
		test_IO


test_Gauge_SOURCES 	= 	test_Gauge.cpp			\
				testerGauge.cpp

test_HMC_nofact_SOURCES = 	test_HMC_noFactories.cpp 	\
			  	testerHMC.cpp

test_HMC_smear_noF_SOURCES = 	test_HMC_Smear_noFact.cpp 	\
			  	testerHMC.cpp

test_IO_SOURCES 	= 	test_IO.cpp 			\
				testerIO.cpp

include_HEADERS 	=	tests.hpp			\
				tests_all.hpp			\
				test_HMC.hpp 			\
				test_HMC_EvenOdd.hpp		\
				test_HMC_DomainWall.hpp		\
				test_wilson.hpp 		\
				test_wilson_EvenOdd.hpp 	\
				test_DomainWall.hpp 		\
				test_ResidualMass.hpp 		\
				test_IO.hpp			\
				test_EigenModes_IRL.hpp		\
				test_smear.hpp			\
				test_Gauge.hpp			\
				test_overlap.hpp		\
				test_MultiShiftSolver.hpp
