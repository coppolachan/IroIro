testdir = $(bindir)/Test

LDADD = -L../lib/ -lcommcode
AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/lib/
AM_CFLAGS  =-I$(top_srcdir) -I$(top_srcdir)/lib/

EXTRA_DIST = 	test_HMC.xml	 \
		test_RNG.xml	 \
		test_Wilson.xml	 \
		test_DomainWall.xml \
		test_ResMass.xml \
		test_IO.xml	

test_PROGRAMS = test_HMC 	\
		test_HMC_nofact \
		test_HMC_smear_noF\
		test_RNG 	\
		test_Wilson 	\
		test_DomainWall	\
		test_ResMass 	\
		test_IO   	\
		test_EigenModes	\
		test_Smearing

test_HMC_SOURCES = 	test_HMC.cpp \
		  	testerHMC.cpp

test_HMC_nofact_SOURCES = 	test_HMC_noFactories.cpp \
			  	testerHMC.cpp

test_HMC_smear_noF_SOURCES = 	test_HMC_Smear_noFact.cpp \
			  	testerHMC.cpp

test_RNG_SOURCES = test_RNG.cpp

test_Wilson_SOURCES =	test_wilson.cpp \
			testerWilson.cpp

test_DomainWall_SOURCES = test_DomainWall.cpp \
			testerODWF.cpp

test_ResMass_SOURCES = test_ResidualMass.cpp \
			testerResMass.cpp

test_IO_SOURCES = test_IO.cpp 	\
		testerIO.cpp

test_EigenModes_SOURCES = test_EigenModes_IRL.cpp\
			testerEigenModes.cpp

test_Smearing_SOURCES = test_smear.cpp	\
			testerSmear.cpp

include_HEADERS =	test_HMC.hpp 		\
			test_wilson.hpp 	\
			test_DomainWall.hpp 	\
			test_ResidualMass.hpp 	\
			test_IO.hpp		\
			test_EigenModes_IRL.hpp	\
			test_smear.hpp
