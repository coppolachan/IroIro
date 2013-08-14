#                         -*- Automake -*-
# Process this file with automake to produce a makefile.
#
# Project IROIRO JLQCD package (KEK) 
# Lattice QCD MonteCarlo simulations
#
# Time-stamp: <2013-08-09 14:32:17 cossu>

###########################################
# Initial definitions 
SUBDIRS = . 

AM_CXXFLAGS = 
AM_LDFLAGS = 
OPTLIB =

DWF_BGQ_SRC = 
SMEAR_SRCS =
###########################################
# Configure driven path selector
if USE_MPI

if IBM_BGQ
# NB: Libraries ordering is essential
#OPTLIB += -lkek -lbgnet -lpthread -lmpich -lopa -lmpl -lpami -lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++ 
OPTLIB += -lkek -lbgnet -lpthread -lmpich-xl -lopa-xl -lmpl-xl -lpami-gcc -lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++ 

COMMUNICATOR_SRCS = Communicator/ArchitectureDep/communicator_mpi_bgq.cpp

if IBM_BGQWILSON_SMP
OPTLIB += -lbgqwilson3x4smp -lbgqthread
else
OPTLIB += -lbgqwilson3x4
endif #IBM_BGQWILSON_SMP

else #not using BGQ routines but MPI

if HITACHI	
AM_CXXFLAGS += -q64 -qarch=pwr7 -qtune=pwr7 -qmaxmem=-1 
AM_LDFLAGS += -q64
endif

COMMUNICATOR_SRCS = Communicator/communicator_mpi.cpp
endif #IBM_BGQ

else
# Fall-back to the single core case
COMMUNICATOR_SRCS = Communicator/communicator_single.cpp	
endif #USE_MPI

############################################
# Sources by architecture
if IBM_BGQKERNELS
DWF_BGQ_SRC += Solver/Architecture_Optimized/rationalSolver_DWF_BGQ.cpp
DWF_BGQ_SRC += Solver/Architecture_Optimized/solver_CG_DWF_BGQ.cpp
ACTION_SRCS = Action/Architecture_Optimized/action_gauge_rect.bgq.cpp 
ACTION_SRCS += Action/Architecture_Optimized/action_gauge_wilson_force.bgq.cpp 	
SMEAR_SRCS += Smearing/Architecture_Optimized/stoutSmear_core_BGQ.cpp
else
ACTION_SRCS = Action/action_gauge_rect.cpp
ACTION_SRCS += Action/action_gauge_wilson_force.cpp
SMEAR_SRCS += Smearing/stoutSmear_core.cpp
endif 
###########################################
# Site indexing (Lexicographic/Even-Odd)
if USE_EVENODD
SITE_INDEXING_SRCS = Geometry/siteIndexEO.cpp
else
SITE_INDEXING_SRCS = Geometry/siteIndex.cpp
endif

###########################################
# Optional DCMT Random number generators
if USE_DCMT
DCMT_SRCS = Tools/RandomNumGen/dcmt_wrapper.cpp
else
DCMT_SRCS = 
endif
###########################################
# Optional for Peter's Bagel/BFM libs
if USE_BFM
BFM_HACK=/bgsys/drivers/V1R1M2/ppc64/
BFM_INCLUDE= -I$(BFM_HACK) -I/dirac1/work/x03/SFW/IROIRO/V1.0/install -I$(BFM_HACK)/spi/include/kernel/cnk/
BFM_LIBS= -L/dirac1/work/x03/SFW/IROIRO/V1.0/install/lib 
BFM_SRCS = Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper.cpp 
BFM_SRCS += Tools/Bagel/bfm_storage.cpp 
BFM_SRCS += Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper_factory.cpp
else
BFM_SRCS = 
endif


AM_CPPFLAGS =-I$(top_srcdir)
AM_CPPFLAGS += $(BFM_INCLUDE)
AM_CFLAGS   =-I$(top_srcdir) 
AM_LDFLAGS  += $(BFM_LIBS)

bin_PROGRAMS = main_hmc
lib_LIBRARIES = libcommcode.a 


EXTRA_DIST =    HMC/mdExec_leapfrog_eigen.cpp	\
		Main/documentation_pages.h	\
                Dirac_ops/dirac_wilson_improved.code    \
                Dirac_ops/dirac_staggered_EvenOdd_improved.code    \
                Dirac_ops/dirac_staggered_EvenOdd_Adjoint_improved.code    \
                Dirac_ops/dirac_wilson_standard.code    	\
                Dirac_ops/dirac_DomainWall_BGQ.code \
                Dirac_ops/domainWallSolver_BGQ.code \
		Fields/field.inl

libcommcode_a_SOURCES = Main/commonPrms.cpp		\
		$(COMMUNICATOR_SRCS)     	        \
		$(SITE_INDEXING_SRCS)			\
		Action/action.cpp			\
		Action/action_gauge_wilson.cpp		\
	 	Action/action_gauge_wilson_adjoint.cpp	\
		$(ACTION_SRCS)				\
	 	Action/action_gaugetype_factory_creator.cpp\
		Action/action_gaugetype_factory.cpp     \
		Action/action_Nf2.cpp			\
		Action/action_Nf2_ratio.cpp		\
		Action/action_Nf.cpp			\
		Action/action_Nf_ratio.cpp		\
		Action/action_staggered.cpp		\
		Action/action_staggered_ratio.cpp	\
		Action/action_Nf2_DomainWall.cpp	\
                Action/action_fermiontype_factory_creator.cpp\
		Action/action_fermiontype_factory.cpp	\
		Communicator/comm_io.cpp		\
		IO/fields_io.cpp			\
		IO/text_header.cpp			\
		IO/binary_reader.cpp			\
		IO/NERSC_reader.cpp			\
		IO/ILDG_reader.cpp			\
		IO/JLQCDLegacy_reader.cpp		\
		Dirac_ops/dirac_Operator_FactoryCreator.cpp \
		Dirac_ops/dirac_Operator_Factory.cpp	\
		Dirac_ops/dirac.cpp     		\
		Dirac_ops/dirac_WilsonLike.cpp		\
		Dirac_ops/dirac_staggeredLike.cpp	\
		Dirac_ops/dirac_wilson.cpp		\
		Dirac_ops/dirac_wilson_EvenOdd.cpp	\
		Dirac_ops/dirac_staggered_EvenOdd.cpp	\
		Dirac_ops/dirac_staggered_EvenOdd_Adjoint.cpp	\
		Dirac_ops/dirac_clover.cpp			\
		Dirac_ops/dirac_Mobius.cpp			\
		Dirac_ops/fopr_signH_Zolotarev.cpp		\
		Dirac_ops/dirac_overlap_Zolotarev.cpp		\
		Dirac_ops/dirac_DWoverlap.cpp    		\
		Dirac_ops/dirac_DomainWall.cpp   		\
		$(DWF_BGQ_SRC)					\
		Dirac_ops/dirac_DomainWall_EvenOdd.cpp  	\
		Dirac_ops/dirac_DomainWall_4D_fullSolv.cpp 	\
		Dirac_ops/dirac_DomainWall_4D_eoSolv.cpp 	\
		Dirac_ops/dirac_wilson_Brillouin.cpp		\
		Dirac_ops/eoUtils.cpp           		\
		Dirac_ops/dirac_Operator_Factory.cpp		\
		Dirac_ops/BoundaryConditions/boundaryCond.cpp	\
		Dirac_ops/subSpaceProjector.cpp	  	        \
		Dirac_ops/dirac_LowModeDeflation_ExactEigen.cpp \
		Dirac_ops/dirac_LowModeDeflation_Approx.cpp     \
		$(BFM_SRCS)					\
		EigenModes/eigenCalcGeneral.cpp			\
		EigenModes/eigenModesSolver_Factory.cpp		\
		EigenModes/eigenModesSolver_IRL.cpp		\
		EigenModes/eigenProc_Zolotarev.cpp		\
		EigenModes/eigenSorter.cpp			\
		EigenModes/findminmax.cpp			\
		Measurements/measGeneral.cpp			\
		Measurements/GaugeM/staples.cpp			\
		Measurements/GaugeM/wilsonFlow.cpp		\
		Measurements/GaugeM/polyakovLoop.cpp			\
		Measurements/GaugeM/gaugeFixingStep.cpp			\
		Measurements/GaugeM/gaugeFixing_Coulomb.cpp		\
		Measurements/GaugeM/gaugeFixing_Landau.cpp		\
		Measurements/GaugeM/gfixFactory.cpp			\
		Measurements/GaugeM/topologyGeom.cpp			\
		Measurements/FermionicM/qprop_MultiShift.cpp 		\
		Measurements/FermionicM/qprop.cpp 			\
		Measurements/FermionicM/qprop_mom.cpp 			\
		Measurements/FermionicM/qprop_EvenOdd.cpp		\
		Measurements/FermionicM/qprop_DomainWall.cpp		\
		Measurements/FermionicM/utils_DWF4d.cpp          	\
		Measurements/FermionicM/quark_prop_meas_factory.cpp 	\
		Measurements/FermionicM/baryonCorrelator.cpp	        \
		Geometry/mapping.cpp			\
		Geometry/mapping_improved.cpp		\
		Geometry/siteIndex_EvenOdd.cpp		\
		Geometry/siteMap.cpp			\
		Geometry/shiftField.cpp			\
		Main/gaugeConf.cpp			\
		Main/gaugeGlobal.cpp			\
		Main/commandline.cpp			\
		HMC/hmcGeneral.cpp                      \
		HMC/mdExec.cpp		                \
		HMC/mdExec_leapfrog.cpp                 \
		HMC/mdExec_2MN.cpp	                \
		HMC/mdExec_factory.cpp                  \
		PugiXML/pugixml.cpp			\
		PugiXML/pugi_interface.cpp		\
		Solver/solver_BiCGStab.cpp		\
		Solver/solver_CG.cpp			\
		Solver/solver_Factory.cpp		\
		Solver/multiShiftSolver_CG.cpp		\
		Solver/rationalSolver_CG.cpp		\
		Solver/solver.cpp			\
		Smearing/APEsmear.cpp			\
		Smearing/stoutSmear_base.cpp		\
		$(SMEAR_SRCS)				\
		Smearing/smartConf.cpp			\
		Smearing/smearingFactories.cpp		\
		Tools/fieldUtils.cpp			\
		Tools/jobUtils.cpp			\
		Tools/sunMatUtils.cpp			\
		Tools/sunAdjMatUtils.cpp		\
		Tools/randNum.cpp			\
		Tools/randNum_MP.cpp			\
		Tools/randNum_Factory.cpp		\
		Tools/gammaMatrices.cpp			\
		Tools/sunAlgUtils.cpp			\
		Tools/complexMatrix.cpp		        \
		Tools/Remez/alg_remez.cpp		\
		Tools/RationalApprox/rationalapprox.cpp	\
		Tools/ErrorsHandling/xml_errors.cpp	\
		Tools/ErrorsHandling/allocation_errors.cpp\
		Tools/ErrorsHandling/base_errors.cpp	\
		Tools/ErrorsHandling/io_errors.cpp	\
		Tools/ErrorsHandling/convergence_errors.cpp\
		Tools/ErrorsHandling/parameter_errors.cpp\
		$(DCMT_SRCS)				\
		Tools/RandomNumGen/randNum_MT19937.cpp	\
		Tools/Timings/timings.cpp		\
		UpdatingAlg/HeatBath/heatbath.cpp	\
		UpdatingAlg/HeatBath/HBActions_SUN.cpp	


main_hmc_SOURCES = Main/main_hmc.cpp 				
main_hmc_LDADD   = libcommcode.a $(OPTLIB)

include_HEADERS = 	$(top_srcdir)/include/iroiro_config.h		\
			$(top_srcdir)/include/commonPrms.hpp		\
                        $(top_srcdir)/include/commandline.hpp           \
			$(top_srcdir)/include/iroiro_code.hpp		\
			$(top_srcdir)/include/common_fields.hpp		\
			$(top_srcdir)/include/config.h			\
			$(top_srcdir)/include/geometry.hpp		\
			$(top_srcdir)/include/field.h			\
			$(top_srcdir)/include/factories.hpp		\
			$(top_srcdir)/include/fopr_Chebyshev.h	        \
			$(top_srcdir)/include/fopr_DdagDLinear.h        \
			$(top_srcdir)/include/fopr_Dhat.h		\
			$(top_srcdir)/include/fopr_eo.h			\
			$(top_srcdir)/include/fopr_eo_wilson.h		\
			$(top_srcdir)/include/fopr.h			\
			$(top_srcdir)/include/fopr_signH_Zolotarev.h	\
			$(top_srcdir)/include/foprHermFactory.hpp	\
			$(top_srcdir)/include/format_A.h		\
			$(top_srcdir)/include/format_F.h		\
			$(top_srcdir)/include/format_G.h		\
			$(top_srcdir)/include/format_Ga.h		\
			$(top_srcdir)/include/format_S.h		\
			$(top_srcdir)/include/numerical_const.hpp	\
			$(top_srcdir)/include/macros.hpp		\
			$(top_srcdir)/include/typeDefs.h		\
			$(top_srcdir)/include/pugi_interface.h		\
			$(top_srcdir)/include/singleton.h		\
			$(top_srcdir)/include/observer.hpp		\
			$(top_srcdir)/include/messages_macros.hpp	\
			$(top_srcdir)/include/timings.hpp		\
			$(top_srcdir)/include/errors.hpp		\
			$(top_srcdir)/include/inputConfig.hpp		\
			PugiXML/pugixml.hpp				\
			PugiXML/pugiconfig.hpp				\
			Action/action_gauge_wilson.hpp			\
	 		Action/action_gauge_wilson_adjoint.hpp		\
			Action/action_gauge_rect.hpp                    \
			Action/action.hpp				\
			Action/action_gaugetype_factory.hpp		\
			Action/action_gaugetype_factory_creator.hpp	\
			Action/action_gaugetype_factory_abs.hpp		\
			Action/action_fermiontype_factory_abs.hpp	\
			Action/action_fermiontype_factory.hpp		\
			Action/action_fermiontype_factory_creator.hpp	\
			Action/action_factory.hpp			\
			Action/action_Nf2.hpp				\
			Action/action_Nf2_ratio.hpp			\
			Action/action_Nf.hpp				\
			Action/action_Nf_ratio.hpp			\
			Action/action_staggered.hpp	  		\
			Action/action_staggered_ratio.hpp		\
			Action/action_Nf2_DomainWall.hpp		\
			Measurements/measGeneral.hpp			\
			Measurements/GaugeM/staples.hpp			\
			Measurements/GaugeM/gaugeFixing.hpp		\
			Measurements/GaugeM/gaugeFixingStep.hpp		\
			Measurements/GaugeM/gaugeFixing_Coulomb.hpp	\
			Measurements/GaugeM/gaugeFixing_Landau.hpp	\
			Measurements/GaugeM/gaugeFixing_Free.hpp	\
			Measurements/GaugeM/gfixFactory.hpp		\
			Measurements/GaugeM/wilsonFlow.hpp		\
			Measurements/GaugeM/polyakovLoop.hpp		\
			Measurements/GaugeM/topologyGeom.hpp		\
			Measurements/FermionicM/source.hpp		\
			Measurements/FermionicM/source_types.hpp	\
			Measurements/FermionicM/qprop_MultiShift.hpp 	\
			Measurements/FermionicM/qprop.hpp		\
			Measurements/FermionicM/qprop_mom.hpp		\
			Measurements/FermionicM/qprop_EvenOdd.hpp	\
			Measurements/FermionicM/qprop_DomainWall.hpp 	\
			Measurements/FermionicM/quark_propagators.hpp	\
			Measurements/FermionicM/sources_factory.hpp	\
			Measurements/FermionicM/meson_correlator.hpp	\
			Measurements/FermionicM/baryonCorrelator.hpp	\
			Measurements/FermionicM/quark_prop_meas_factory.hpp \
                        Measurements/FermionicM/utils_DWF4d.hpp \
			Communicator/communicator.hpp		\
			Communicator/comm_io.hpp		\
			IO/fields_io.hpp			\
			IO/binary_reader.hpp			\
			IO/general_reader.hpp			\
			IO/generic_header.hpp			\
			IO/ILDG_header.hpp			\
			IO/ILDG_reader.hpp			\
			IO/JLQCDLegacy_reader.hpp		\
			IO/NERSC_reader.hpp			\
			IO/text_header.hpp			\
			Geometry/siteIndex.hpp     		\
			Geometry/siteIndex_EvenOdd.hpp		\
			Geometry/siteMap.hpp			\
			Geometry/mapping.hpp			\
			Geometry/shiftField.hpp			\
			Main/gaugeConf.hpp			\
			Main/gaugeGlobal.hpp			\
		        Dirac_ops/dirac_Operator_FactoryCreator.hpp  \
			Dirac_ops/dirac_Operator_Factory.hpp	\
			Dirac_ops/dirac_DWoverlap.hpp 		\
			Dirac_ops/dirac_overlap_Zolotarev.hpp	\
			Dirac_ops/dirac.hpp			\
			Dirac_ops/dirac_WilsonLike.hpp		\
			Dirac_ops/dirac_staggeredLike.hpp	\
			Dirac_ops/dirac_wilson.hpp		\
			Dirac_ops/dirac_wilson_EvenOdd.hpp	\
			Dirac_ops/dirac_wilson_Brillouin.hpp	\
			Dirac_ops/dirac_staggered_EvenOdd.hpp	\
			Dirac_ops/dirac_staggered_EvenOdd_Adjoint.hpp	\
			Dirac_ops/dirac_clover.hpp		\
			Dirac_ops/dirac_Mobius.hpp		\
			Dirac_ops/dirac_DomainWall.hpp   	\
			Dirac_ops/dirac_DomainWall_EvenOdd.hpp	\
			Dirac_ops/dirac_DomainWall_4D_fullSolv.hpp \
			Dirac_ops/dirac_DomainWall_4D_eoSolv.hpp \
			Dirac_ops/dirac_LowModeDeflation.hpp    \
			Dirac_ops/dirac_LowModeDeflation_ExactEigen.hpp \
			Dirac_ops/dirac_LowModeDeflation_Approx.hpp     \
			Dirac_ops/eoUtils.hpp            	\
			Dirac_ops/BoundaryConditions/boundaryCond.hpp      	\
			Dirac_ops/BoundaryConditions/antiPeriodicBC.hpp      	\
			Dirac_ops/subSpaceProjector.hpp		\
			Dirac_ops/utils_BGQ.hpp	        	\
			Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper.hpp\
			Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper_factory.hpp\
                        EigenModes/eigenCalcGeneral.hpp  	\
			EigenModes/eigenModesSolver_Factory.hpp	\
			EigenModes/eigenModesSolver.hpp  	\
			EigenModes/eigenModesSolver_IRL.hpp	\
			EigenModes/eigenSorter_Factory.hpp	\
			EigenModes/eigenSorter.hpp		\
			EigenModes/foprHermFactory_ChebyshevDdagDLin.hpp\
			EigenModes/eigenProc_Zolotarev.h	\
			EigenModes/eigenModes.hpp		\
			EigenModes/findminmax.hpp		\
			Fields/field_expressions.hpp		\
			HMC/hmcGeneral.hpp                      \
			HMC/mdExec.hpp		                \
			HMC/mdExec_leapfrog.hpp			\
			HMC/mdExec_2MN.hpp	                \
			HMC/mdExec_leapfrog_eigen.hpp		\
			HMC/mdExec_factory.hpp			\
			HMC/mdExec_factory_abs.hpp		\
			Solver/solver_CG_params.hpp		\
			Solver/solver_Factory.hpp		\
			Solver/solver_BiCGStab.hpp		\
			Solver/solver.hpp			\
			Solver/solver_CG.hpp			\
			Solver/multiShiftSolver_CG.hpp		\
			Solver/multiShiftSolver.hpp		\
			Solver/rationalSolver.hpp		\
			Solver/rationalSolver_CG.hpp		\
			Solver/Architecture_Optimized/rationalSolver_DWF_BGQ.hpp\
			Solver/Architecture_Optimized/solver_CG_DWF_BGQ.hpp\
			Smearing/APEsmear.hpp			\
			Smearing/baseSmear.hpp			\
			Smearing/stoutSmear.hpp			\
			Smearing/smartConf.hpp			\
			Smearing/smearingFactories.hpp		\
			Tools/jobUtils.hpp			\
			Tools/sunMat.hpp			\
			Tools/sunMatUtils.hpp   		\
			Tools/sunAdjMat.hpp			\
			Tools/sunAdjVec.hpp			\
			Tools/sunAdjMatUtils.hpp		\
			Tools/sunRepresentations.hpp		\
			Tools/sunVec.hpp			\
			Tools/sunAlgUtils.hpp			\
			Tools/colorEpsilon.hpp	                \
			Tools/randNum_MP.h			\
			Tools/randNum.h				\
			Tools/randNum_Factory.h			\
			Tools/fieldUtils.hpp			\
			Tools/RAIIFactory.hpp			\
			Tools/EnumToString.hpp			\
			Tools/byteswap.hpp			\
			Tools/gammaMatrices.hpp			\
			Tools/complexMatrix.hpp 		\
			Tools/Remez/alg_remez.hpp		\
			Tools/Remez/bigfloat.hpp		\
			Tools/RationalApprox/rationalapprox.hpp	\
			Tools/Bagel/bfm_storage.hpp		\
			Tools/RandomNumGen/dcmt_wrapper.hpp	\
                        Tools/RandomNumGen/randNum_MT19937.h    \
			UpdatingAlg/HeatBath/heatbath.hpp	\
			UpdatingAlg/HeatBath/HBActions.hpp	\
			UpdatingAlg/HeatBath/HBActions_SUN.hpp	





