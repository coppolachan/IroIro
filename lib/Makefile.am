SUBDIRS = . 

if USE_MPI
COMMUNICATOR_SRCS = Communicator/communicator_mpi.cpp
if HITACHI	
CXXFLAGS += -q64 -qarch=pwr7 -qtune=pwr7 -qmaxmem=-1 
LDFLAGS += -q64
else
CXXFLAGS += 
LDFLAGS += 
endif
else
COMMUNICATOR_SRCS = Communicator/communicator_dummy.cpp	
endif

AM_CPPFLAGS=-I$(top_srcdir)
AM_CFLAGS  =-I$(top_srcdir)
AM_LDFLAGS =

bin_PROGRAMS = main_hmc

lib_LIBRARIES = libcommcode.a

EXTRA_DIST =    HMC/mdExec_leapfrog_eigen.cpp	\
		../tests/test_HMC_noFactories.cpp \
		Main/documentation_pages.h\
		Fields/field.inl


libcommcode_a_SOURCES = Main/commonPrms.cpp		\
		$(COMMUNICATOR_SRCS)     	        \
		Action/action.cpp			\
	 	Action/action_gauge_wilson.cpp 		\
	 	Action/action_gauge_rect.cpp 		\
		Action/action_Nf2.cpp			\
		Action/action_Nf2_ratio.cpp		\
		Action/action_gaugetype_factory.cpp	\
		Action/action_fermiontype_factory.cpp	\
		Communicator/comm_io.cpp		\
		Communicator/fields_io.cpp		\
		Dirac_ops/dirac_wilson.cpp		\
		Dirac_ops/dirac_clover.cpp		\
		Dirac_ops/dirac_wilson_EvenOdd.cpp	\
		Dirac_ops/fopr_signH_Zolotarev.cpp	\
		Dirac_ops/dirac_overlap_Zolotarev.cpp	\
		Dirac_ops/dirac_Operator_Factory.cpp	\
		Dirac_ops/dirac_DomainWall.cpp   	\
		Dirac_ops/dirac_DomainWall_EvenOdd.cpp  \
		Dirac_ops/dirac_DomainWall_4D.cpp	\
		EigenModes/eigenModes_IRL.cpp		\
		EigenModes/eigenProc_Zolotarev.cpp	\
		EigenModes/sortEigen.cpp		\
		EigenModes/fopr_chebyshev_DdagD.cpp	\
		Solver/solver_BiCGStab.cpp		\
		Solver/solver_CG.cpp			\
		Solver/solver_Factory.cpp		\
		Solver/multiShiftSolver_CG.cpp		\
		Solver/solver.cpp			\
		Smearing/APEsmear.cpp			\
		Smearing/StoutSmear.cpp			\
		Smearing/SmartConf.cpp			\
		Smearing/Smearing_Factories.cpp		\
		HMC/hmcGeneral.cpp			\
		HMC/mdExec.cpp		        	\
		HMC/mdExec_leapfrog.cpp			\
		HMC/mdExec_Factory.cpp			\
		Measurements/GaugeM/staples.cpp		\
		Measurements/FermionicM/qprop_MultiShift.cpp 		\
		Measurements/FermionicM/qprop.cpp 			\
		Measurements/FermionicM/qprop_EvenOdd.cpp		\
		Measurements/FermionicM/qprop_DomainWall.cpp		\
		Measurements/FermionicM/fermion_meas_factory.cpp 	\
		Main/Geometry/siteIndex.cpp		\
		Main/Geometry/siteIndex_eo.cpp		\
		Main/gaugeConf.cpp			\
		Main/commandline.cpp			\
		PugiXML/pugixml.cpp			\
		PugiXML/pugi_interface.cpp		\
		Tools/sunMatUtils.cpp			\
		Tools/randNum_MT19937.cpp		\
		Tools/randNum_Factory.cpp		\
		Tools/gammaMatrices.cpp			\
		../tests/test_MultiShiftSolver.cpp	\
		../tests/test_overlap.cpp		\
		../tests/test_DomainWall.cpp		\
		../tests/test_HMC_DomainWall.cpp	\
                ../tests/test_HMC_EvenOdd.cpp           \
		$(top_srcdir)/include/common_fields.cpp

main_hmc_SOURCES = Main/main_hmc.cpp 				
main_hmc_LDADD   = libcommcode.a

include_HEADERS = 	$(top_srcdir)/include/commonPrms.h		\
                        $(top_srcdir)/include/commandline.hpp           \
			$(top_srcdir)/include/common_code.hpp		\
			$(top_srcdir)/include/common_fields.hpp		\
			$(top_srcdir)/include/config.h			\
			$(top_srcdir)/include/geometry.hpp		\
			$(top_srcdir)/include/field.h			\
			$(top_srcdir)/include/factories.hpp		\
			$(top_srcdir)/include/fopr_chebyshev_DdagD.h	\
			$(top_srcdir)/include/fopr_Dhat.h		\
			$(top_srcdir)/include/fopr_eo.h			\
			$(top_srcdir)/include/fopr_eo_wilson.h		\
			$(top_srcdir)/include/fopr.h			\
			$(top_srcdir)/include/fopr_signH_Zolotarev.h	\
			$(top_srcdir)/include/format_A.h		\
			$(top_srcdir)/include/format_F.h		\
			$(top_srcdir)/include/format_G.h		\
			$(top_srcdir)/include/numerical_const.hpp	\
			$(top_srcdir)/include/macros.hpp		\
			$(top_srcdir)/include/typeDefs.h		\
			$(top_srcdir)/include/pugi_interface.h		\
			$(top_srcdir)/include/singleton.h		\
			$(top_srcdir)/include/observer.hpp		\
			PugiXML/pugixml.hpp				\
			PugiXML/pugiconfig.hpp				\
			Action/action_gauge_wilson.hpp			\
			Action/action_gauge_rect.hpp			\
			Action/action.hpp				\
			Action/action_gaugetype_factory.hpp		\
			Action/action_fermiontype_factory.hpp		\
			Action/action_Factory.hpp			\
			Action/action_Nf2.hpp				\
			Action/action_Nf2_ratio.hpp			\
			Action/action_Nf2_DomainWall.hpp		\
			EigenModes/eigenModes_IRL.hpp			\
			EigenModes/sortEigen.h				\
			EigenModes/eigenProc_Zolotarev.h		\
			Measurements/GaugeM/staples.h				\
			Measurements/FermionicM/source.hpp			\
			Measurements/FermionicM/source_types.hpp		\
			Measurements/FermionicM/qprop_MultiShift.hpp 		\
			Measurements/FermionicM/qprop.h				\
			Measurements/FermionicM/qprop_EvenOdd.h			\
			Measurements/FermionicM/qprop_DomainWall.hpp 		\
			Measurements/FermionicM/fermion_meas_factory.hpp 	\
			Measurements/FermionicM/quark_propagators.hpp		\
			Measurements/FermionicM/sources_factory.hpp		\
			Measurements/FermionicM/meson_correlator.hpp		\
			Communicator/communicator.h		\
			Communicator/comm_io.hpp		\
			Communicator/fields_io.hpp		\
			Main/Geometry/siteIndex.h		\
			Main/Geometry/shiftField.h		\
			Main/Geometry/shiftField_eo.h		\
			Main/Geometry/siteIndex_eo.h		\
			Main/gaugeConf.hpp			\
			Main/commandline.cpp			\
			Dirac_ops/dirac_Operator_Factory.hpp	\
			Dirac_ops/dirac_overlap_Zolotarev.h	\
			Dirac_ops/dirac_wilson.h		\
			Dirac_ops/dirac_clover.hpp		\
			Dirac_ops/dirac_wilson_EvenOdd.h	\
			Dirac_ops/dirac.h			\
			Dirac_ops/dirac_DomainWall.hpp   	\
			Dirac_ops/dirac_DomainWall_EvenOdd.hpp	\
			Dirac_ops/dirac_DomainWall_4D.hpp	\
			Dirac_ops/dirac_Preconditioners.hpp	\
			Fields/field_expressions.hpp		\
			Solver/solver_Factory.hpp		\
			Solver/solver_BiCGStab.h		\
			Solver/solver.hpp			\
			Solver/solver_CG.h			\
			Solver/multiShiftSolver_CG.hpp		\
			Solver/multiShiftSolver.hpp		\
			Smearing/APEsmear.hpp			\
			Smearing/BaseSmear.hpp			\
			Smearing/StoutSmear.hpp			\
			Smearing/SmartConf.hpp			\
			Smearing/Smearing_Factories.hpp		\
			HMC/hmcGeneral.hpp			\
			HMC/mdExec.h				\
			HMC/mdExec_leapfrog.hpp			\
			HMC/mdExec_leapfrog_eigen.h		\
			HMC/mdExec_Factory.hpp			\
			Tools/sunMat.h				\
			Tools/sunMatUtils.hpp			\
			Tools/randNum_MP.h			\
			Tools/randNum.h				\
			Tools/sunVec.h				\
			Tools/randNum_MT19937.h			\
			Tools/randNum_Factory.h			\
			Tools/RAIIFactory.hpp			\
			Tools/EnumToString.hpp			\
			Tools/byteswap.hpp			\
			Tools/gammaMatrices.hpp			\
			../tests/tests_all.hpp			\
			../tests/test_MultiShiftSolver.h	\
			../tests/test_overlap.hpp		\
                        ../tests/test_HMC_DomainWall.hpp        \
                        ../tests/test_HMC_EvenOdd.hpp           \
			../tests/tests.hpp                      




