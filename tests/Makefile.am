#                   -*- Automake -*-
# Process this file with automake to produce a makefile.
#
# Project IROIRO JLQCD package (KEK) 
# Lattice QCD MonteCarlo simulations Test files
#
# Time-stamp: <2013-08-09 11:55:44 cossu>

testdir = $(bindir)/Test
OPTLIB =
BMF_TESTPROGS = 
AM_CXXFLAGS =
AM_LDFLAGS = 

if USE_BFM
BFM_HACK=/bgsys/drivers/V1R1M2/ppc64/
BFM_INCLUDE= -I$(BFM_HACK) -I/dirac1/work/x03/SFW/IROIRO/V1.0/install -I$(BFM_HACK)/spi/include/kernel/cnk/
BFM_LIBS= -L/dirac1/work/x03/SFW/IROIRO/V1.0/install/lib 
BFM_LINKAGE= -lbfm -lbagel
BMF_TESTPROGS += test_Solver_BFM
BMF_TESTPROGS += test_BFM_Operator
BMF_TESTPROGS += test_BFM_Moebius
BMF_TESTPROGS += test_BFM_Inverter
BMF_TESTPROGS += test_BFM_Inverter_wrapper
endif

if IBM_BGQ
#OPTLIB += -lkek -lbgnet -lpthread -lmpich -lopa -lmpl -lpami -lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++
OPTLIB += -lkek -lbgnet -lpthread -lmpich-xl -lopa-xl -lmpl-xl -lpami-gcc -lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++ 

if IBM_BGQWILSON_SMP
AM_CXXFLAGS += -qsmp=omp:auto
OPTLIB += -lbgqwilson3x4smp -lbgqthread
else
OPTLIB += -lbgqwilson3x4
endif

endif

AM_CXXFLAGS += $(BFM_INCLUDE) 
AM_LDFLAGS += $(BFM_LIBS)

LDADD = -L../lib/ -lcommcode  $(BFM_LINKAGE) $(OPTLIB)
AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/lib/
AM_CFLAGS  =-I$(top_srcdir) -I$(top_srcdir)/lib/

EXTRA_DIST = 	testResults.info	\
		test_HMC.clover.xml	\
		test_HMC.wilson.xml	\
		test_HMC.wilson.smear.xml \
		test_HMC.wilson.rational.xml	\
		test_HMC.wilson.ratio.rational.xml	\
		test_HMC.wilson.ratio.xml \
		test_HMC.wilson.ratio.smeared.xml \
		test_HMC.wilson.rational.smeared.xml \
		test_HMC.DWF.rational.xml	\
		test_RNG.xml	 	\
		test_Wilson.xml	 	\
		test_Wilson_EvenOdd.xml	\
		test_Wilson_Brillouin.xml	\
		test_DomainWall.xml 	\
		test_QpropMom.xml 	\
		test_EigenModesSolver.xml\
		test_GeneralMeas.xml 	\
		test_PolyakovLoop.xml	\
		test_WilsonFlow.xml	\
		test_MesonSpectrum_Nf2p1.xml 	\
		test_BaryonSpectrum_Nf2p1.xml 	\
		test_HadronSpectrum_Nf2p1.xml 	\
		test_ResMass_full.xml 	\
		test_ResMass_eo.xml 	\
		test_IO.xml		\
		test_Gauge.xml		\
		test_NoMeas.xml	        \
		test_Overlap.xml	\
		test_HMC_DWF.xml	\
		test_MultiShiftSolver.xml	\
		test_RationalApprox.xml	\
		test_Solver_Perf.xml	\
		test_sunRep.xml		\
		test_noFact.xml

test_PROGRAMS = test_HMC 		\
		test_HMC_noFactories 	\
		test_HMC_smear_noF	\
		test_HMC_DWall_noF	\
		test_Gauge		\
		test_NoMeas		\
		test_PolyakovLoop	\
		test_WilsonFlow		\
		test_Wilson 		\
		test_Clover 		\
		test_Wilson_EvenOdd 	\
		test_Wilson_Brillouin 	\
		test_Wilson_Brillouin_Imp\
		test_QpropMom		\
		test_ResMass 		\
		test_MesonSpectrum	\
		test_MesonSpectrum_SpMom\
		test_MesonSpectrum_Nf2p1\
		test_BaryonSpectrum_Nf2p1\
		test_HadronSpectrum_Nf2p1\
		test_EigenModesSolver	\
		test_Smearing		\
		test_Overlap		\
		test_MultiShiftSolver	\
		test_Solver_Perf	\
		$(BMF_TESTPROGS)        \
		test_RNG 		\
		test_IO			\
		test_RationalApprox	\
		test_RHMC_noFactories 	\
		test_SUNRepresentation	\
		test_HeatBath           \
                test_DWF              	

test_HMC_SOURCES 	= 	test_HMC.cpp 			\
		  		testerHMC.cpp

test_HMC_noFactories_SOURCES = 	test_HMC_noFactories.cpp 	\
			  	testerHMC.cpp

test_RHMC_noFactories_SOURCES = test_RHMC_noFactories.cpp 	\
			  	testerHMC.cpp

test_HMC_smear_noF_SOURCES = 	test_HMC_Smear_noFact.cpp 	\
			  	testerHMC.cpp

test_HMC_DWall_noF_SOURCES =	test_HMC_DWall_noFact.cpp	\
				testerHMC_DWall_noFact.cpp

test_NoMeas_SOURCES 	= 	testerNoMeas.cpp

test_Gauge_SOURCES 	= 	test_Gauge.cpp			\
				testerGauge.cpp

test_Wilson_SOURCES 	=	test_wilson.cpp 		\
				testerWilson.cpp

test_Clover_SOURCES 	=	test_clover.cpp 		\
				testerClover.cpp

test_Wilson_EvenOdd_SOURCES =	test_wilson_EvenOdd.cpp 	\
		        	testerWilson_EvenOdd.cpp

test_Wilson_Brillouin_SOURCES =	test_wilson_Brillouin.cpp 	\
		        	testerWilsonBrillouin.cpp

test_Wilson_Brillouin_Imp_SOURCES = test_wilson_Brillouin_Imp.cpp 	\
		        	testerWilsonBrillouin_Imp.cpp

test_QpropMom_SOURCES 	= 	test_QpropMom.cpp 		\
				testerQpropMom.cpp

test_MesonSpectrum_SOURCES= 	test_MesonSpectrum.cpp 		\
				testerMesonSpectrum.cpp

test_MesonSpectrum_SpMom_SOURCES= test_MesonSpectrum_SpMom.cpp 	\
				testerMesonSpectrumSpMom.cpp

test_MesonSpectrum_Nf2p1_SOURCES= test_MesonSpectrum_Nf2p1.cpp 	\
				testerMesonSpectrumNf2p1.cpp

test_BaryonSpectrum_Nf2p1_SOURCES= test_BaryonSpectrum_Nf2p1.cpp 	\
				testerBaryonSpectrumNf2p1.cpp

test_HadronSpectrum_Nf2p1_SOURCES= test_HadronSpectrum_Nf2p1.cpp 	\
				testerHadronSpectrumNf2p1.cpp

test_Solver_Perf_SOURCES = 	test_Solver_Perf.cpp		\
				testerSolver_Perf.cpp

test_ResMass_SOURCES 	= 	test_ResidualMass.cpp 		\
				testerResMass.cpp

test_EigenModesSolver_SOURCES =	test_EigenModesSolver.cpp	\
				testerEigenModesSolver.cpp

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

test_WilsonFlow_SOURCES 	= 	test_WilsonFlow.cpp	\
					testerWilsonFlow.cpp
test_PolyakovLoop_SOURCES 	= 	test_PolyakovLoop.cpp	\
	 				testerPolyakovLoop.cpp

test_SUNRepresentation_SOURCES  =	test_sunRepres.cpp	\
					testerSunRep.cpp	

test_HeatBath_SOURCES 		=	test_HeatBath.cpp	\
					testerHeatBath.cpp

test_DWF_SOURCES 		=	test_DWF.cpp	\
					testerDWF.cpp

## Put all BFM related tests here
if USE_BFM
test_Solver_BFM_SOURCES =       test_Solver_BFM.cpp             \
                                testerSolver_BFM.cpp            

test_BFM_Operator_SOURCES =     test_BFM_Operator.cpp             \
                                testerSolver_BFM.cpp            
test_BFM_Moebius_SOURCES =     test_BFM_Moebius.cpp             \
                                testerSolver_BFM.cpp            
test_BFM_Inverter_SOURCES =     test_BFM_Inverter.cpp 		\
                                testerSolver_BFM.cpp            

test_BFM_Inverter_wrapper_SOURCES = test_BFM_Inverter_wrapper.cpp	\
                                testerSolver_BFM.cpp            

endif


include_HEADERS 	=	tests.hpp			\
				tests_all.hpp			\
				test_HMC.hpp 			\
				test_HMC_DWall_noFact.hpp	\
				test_wilson.hpp 		\
				test_clover.hpp 		\
				test_wilson_EvenOdd.hpp 	\
				test_wilson_Brillouin.hpp 	\
				test_wilson_Brillouin_Imp.hpp 	\
				test_WilsonFlow.hpp     	\
				test_PolyakovLoop.hpp     	\
				test_NoMeas.hpp     	        \
				test_QpropMom.hpp 		\
				test_ResidualMass.hpp 		\
				test_IO.hpp			\
				test_EigenModesSolver.hpp	\
				test_smear.hpp			\
				test_Gauge.hpp			\
				test_MesonSpectrum.hpp		\
				test_MesonSpectrum_Nf2p1.hpp	\
	                        test_MesonSpectrum_SpMom.hpp    \
				test_BaryonSpectrum_Nf2p1.hpp	\
				test_HadronSpectrum_Nf2p1.hpp	\
				test_Overlap.hpp		\
				test_MultiShiftSolver.hpp	\
				test_RationalApprox.hpp		\
	                        test_Solver_Perf.hpp		\
				test_sunRepres.hpp		\
				test_HeatBath.hpp               \
                                test_DWF.hpp
