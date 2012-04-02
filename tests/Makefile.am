# Time-stamp: <2012-03-29 17:48:37 neo>

testdir = $(bindir)/Test

LDADD = -L../lib/ -lcommcode
AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/lib/
AM_CFLAGS  =-I$(top_srcdir) -I$(top_srcdir)/lib/

EXTRA_DIST = 	testResults.info	\
		test_HMC.xml	 	\
		test_HMC.clover.xml	\
		test_RNG.xml	 	\
		test_Wilson.xml	 	\
		test_Wilson_EvenOdd.xml	\
		test_DomainWall.xml 	\
		test_QpropMom.xml 	\
		test_ResMass_full.xml 	\
		test_ResMass_eo.xml 	\
		test_IO.xml		\
		test_Gauge.xml		\
		test_Overlap.xml	\
		test_HMC_EvenOdd.xml	\
		test_HMC.DWF.zolo.xml	\
		test_MultiShiftSolver.xml
		test_RationalApprox.xml	\
		input.xml

test_PROGRAMS = test_HMC 		\
		test_HMC_noFactories 	\
		test_RHMC_noFactories 	\
		test_HMC_EvenOdd	\
		test_HMC_smear_noF	\
		test_HMC_DWall_noF	\
		test_Gauge		\
		test_Wilson 		\
		test_Clover 		\
		test_Wilson_EvenOdd 	\
		test_DomainWall		\
		test_QpropMom		\
		test_ResMass 		\
		test_EigenModes		\
		test_Smearing		\
		test_Overlap		\
		test_MultiShiftSolver	\
		test_RNG 		\
		test_IO			\
		test_RationalApprox


test_HMC_SOURCES 	= 	test_HMC.cpp 			\
		  		testerHMC.cpp

test_HMC_noFactories_SOURCES = 	test_HMC_noFactories.cpp 	\
			  	testerHMC.cpp

test_RHMC_noFactories_SOURCES = test_RHMC_noFactories.cpp 	\
			  	testerHMC.cpp

test_HMC_smear_noF_SOURCES = 	test_HMC_Smear_noFact.cpp 	\
			  	testerHMC.cpp

test_HMC_EvenOdd_SOURCES  =	test_HMC_EvenOdd.cpp		\
				testerHMCeven-odd.cpp

test_HMC_DWall_noF_SOURCES =	test_HMC_DWall_noFact.cpp	\
				testerHMC_DWall_noFact.cpp

test_Gauge_SOURCES 	= 	test_Gauge.cpp			\
				testerGauge.cpp

test_Wilson_SOURCES 	=	test_wilson.cpp 		\
				testerWilson.cpp

test_Clover_SOURCES 	=	test_clover.cpp 		\
				testerClover.cpp

test_Wilson_EvenOdd_SOURCES =	test_wilson_EvenOdd.cpp 	\
		        	testerWilson_EvenOdd.cpp

test_DomainWall_SOURCES = 	test_DomainWall.cpp 		\
				testerDWF.cpp

test_QpropMom_SOURCES 	= 	test_QpropMom.cpp 		\
				testerQpropMom.cpp

test_ResMass_SOURCES 	= 	test_ResidualMass.cpp 		\
				testerResMass.cpp

test_EigenModes_SOURCES = 	test_EigenModes_IRL.cpp		\
				testerEigenModes.cpp

test_Overlap_SOURCES 	= 	test_Overlap.cpp		\
				testerOverlap.cpp

test_MultiShiftSolver_SOURCES =	test_MultiShiftSolver.cpp	\
				testerMultiShiftSolver.cpp

test_Smearing_SOURCES   =       test_smear.cpp                  \
                                testerSmear.cpp

test_RNG_SOURCES 	= 	test_RNG.cpp

test_IO_SOURCES 	= 	test_IO.cpp 			\
				testerIO.cpp

test_RationalApprox_SOURCES 	= 	test_RationalApprox.cpp	\
				testerRationalApprox.cpp


include_HEADERS 	=	tests.hpp			\
				tests_all.hpp			\
				test_HMC.hpp 			\
				test_HMC_EvenOdd.hpp		\
				test_HMC_DWall_noFact.hpp	\
				test_wilson.hpp 		\
				test_clover.hpp 		\
				test_wilson_EvenOdd.hpp 	\
				test_DomainWall.hpp 		\
				test_QpropMom.hpp 		\
				test_ResidualMass.hpp 		\
				test_IO.hpp			\
				test_EigenModes_IRL.hpp		\
				test_smear.hpp			\
				test_Gauge.hpp			\
				test_Overlap.hpp		\
				test_MultiShiftSolver.hpp	\
				test_RationalApprox.hpp
