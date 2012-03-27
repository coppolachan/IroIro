SUBDIRS = . 

CXXFLAGS += 
LDFLAGS += 
OPTLIB =

if USE_MPI
COMMUNICATOR_SRCS = Communicator/communicator_mpi.cpp
MAPPER_SRCS       = Main/Geometry/map_parallelscalar.cpp
MAPPER_HEAD       = Main/Geometry/map_parallelscalar.hpp

if HITACHI      
CXXFLAGS += -q64 -qarch=pwr7 -qtune=pwr7 -qmaxmem=-1 
LDFLAGS += -q64
endif

if IBM_WATSON
DRIVER=/bgsys/drivers/ppcfloor
CXXFLAGS += -q64 -qarch=qp -qtune=qp -I$(DRIVER) -I$(DRIVER)/spi/include/kernel/cnk \
	-I$(DRIVER)/spi/include/mu/default -I$(DRIVER)/comm/xl/include 
LDFLAGS += -L$(DRIVER)/comm/xl/lib -L$(DRIVER)/comm/sys/lib -L$(DRIVER)/spi/lib
OPTLIB += -lpthread -lmpich -lopa -lmpl -lpami -lSPI_cnk -lSPI -lSPI_l1p -lrt -lstdc++
endif

else
COMMUNICATOR_SRCS = Communicator/communicator_dummy.cpp 
MAPPER_SRCS       = Main/Geometry/map_scalar.cpp
MAPPER_HEAD       = Main/Geometry/map_scalar.hpp
endif

AM_CPPFLAGS=-I$(top_srcdir)
AM_CFLAGS  =-I$(top_srcdir)
AM_LDFLAGS =

bin_PROGRAMS = main_hmc

lib_LIBRARIES = libcommcode.a

EXTRA_DIST =    HMC/mdExec_leapfrog_eigen.cpp		\
		Main/documentation_pages.h		\
		Dirac_ops/dirac_wilson_improved.code	\
		Dirac_ops/dirac_wilson_improved2.code	\
		Dirac_ops/dirac_wilson_improved3.code	\
		Dirac_ops/dirac_wilson_standard.code	\
		Fields/field.inl


libcommcode_a_SOURCES = Main/commonPrms.cpp		\
		Main/common_fields.cpp			\
		$(COMMUNICATOR_SRCS)     	        \
		Action/action.cpp			\
	 	Action/action_gauge_wilson.cpp 		\
		Action/action_gauge_rect.cpp            \
		Action/action_gaugetype_factory.cpp     \
		Action/action_Nf2.cpp			\
		Action/action_Nf2_ratio.cpp		\
		Action/action_Nf2_DomainWall.cpp	\
		Action/action_fermiontype_factory.cpp	\
		Communicator/comm_io.cpp		\
		Communicator/fields_io.cpp		\
		Dirac_ops/dirac_wilson.cpp		\
		Dirac_ops/dirac_wilson_EvenOdd.cpp	\
		Dirac_ops/dirac_clover.cpp		\
		Dirac_ops/fopr_signH_Zolotarev.cpp	\
		Dirac_ops/dirac_overlap_Zolotarev.cpp	\
		Dirac_ops/dirac_DomainWall.cpp   	\
		Dirac_ops/dirac_DomainWall_EvenOdd.cpp  \
		Dirac_ops/dirac_DomainWall_4D_fullSolv.cpp \
		Dirac_ops/dirac_DomainWall_4D_eoSolv.cpp \
		Dirac_ops/eoUtils.cpp           	\
		Dirac_ops/dirac_Operator_Factory.cpp	\
		EigenModes/eigenModes_IRL.cpp		\
		EigenModes/eigenProc_Zolotarev.cpp	\
		EigenModes/sortEigen.cpp		\
		EigenModes/fopr_chebyshev_DdagD.cpp	\
		Measurements/GaugeM/staples.cpp		\
		Main/Geometry/siteIndex.cpp		\
		Main/Geometry/siteIndex_eo.cpp		\
		Main/Geometry/mapper.cpp		\
		$(MAPPER_SRCS)		     	        \
		Main/gaugeConf.cpp			\
		Main/commandline.cpp			\
		Measurements/FermionicM/qprop_MultiShift.cpp 		\
		Measurements/FermionicM/qprop.cpp 			\
		Measurements/FermionicM/qprop_EvenOdd.cpp		\
		Measurements/FermionicM/qprop_DomainWall.cpp		\
		Measurements/FermionicM/qprop_mom.cpp 			\
		Measurements/FermionicM/utils_DWF4d.cpp			\
		Measurements/FermionicM/quark_prop_meas_factory.cpp 	\
		HMC/hmcGeneral.cpp                      \
		HMC/mdExec.cpp		                \
		HMC/mdExec_leapfrog.cpp                 \
		HMC/mdExec_factory.cpp                  \
		PugiXML/pugixml.cpp			\
		PugiXML/pugi_interface.cpp		\
		Solver/solver_BiCGStab.cpp		\
		Solver/solver_CG.cpp			\
		Solver/solver_Factory.cpp		\
		Solver/multiShiftSolver_CG.cpp		\
		Solver/solver.cpp			\
		Smearing/APEsmear.cpp			\
		Smearing/StoutSmear.cpp			\
		Smearing/SmartConf.cpp			\
		Smearing/Smearing_Factories.cpp		\
		Tools/sunMatUtils.cpp			\
		Tools/randNum_MT19937.cpp		\
		Tools/randNum_Factory.cpp		\
		Tools/gammaMatrices.cpp			\
		Tools/Remez/alg_remez.cpp


main_hmc_SOURCES = Main/main_hmc.cpp 				
main_hmc_LDADD   = libcommcode.a $(OPTLIB)

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
			$(top_srcdir)/include/messages_macros.hpp	\
			PugiXML/pugixml.hpp				\
			PugiXML/pugiconfig.hpp				\
			Action/action_gauge_wilson.hpp			\
			Action/action_gauge_rect.hpp                    \
			Action/action.hpp				\
			Action/action_gaugetype_factory.hpp		\
			Action/action_gaugetype_factory_abs.hpp		\
			Action/action_fermiontype_factory_abs.hpp	\
			Action/action_fermiontype_factory.hpp		\
			Action/action_Factory.hpp			\
			Action/action_Nf2.hpp				\
			Action/action_Nf2_ratio.hpp			\
			Action/action_Nf2_DomainWall.hpp		\
			Measurements/GaugeM/staples.hpp			\
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
			Measurements/FermionicM/fermion_meas_factory_abs.hpp	\
			Measurements/FermionicM/quark_prop_meas_factory.hpp	\
                        Measurements/FermionicM/utils_DWF4d.hpp \
			Communicator/communicator.h		\
			Communicator/comm_io.hpp		\
			Communicator/fields_io.hpp		\
			Dirac_ops/dirac.hpp			\
			Dirac_ops/dirac_wilson.hpp		\
			Dirac_ops/dirac_wilson_EvenOdd.hpp	\
			Dirac_ops/dirac_clover.hpp		\
			Dirac_ops/dirac_DomainWall.hpp   	\
			Dirac_ops/dirac_DomainWall_EvenOdd.hpp	\
			Dirac_ops/dirac_DomainWall_4D_fullSolv.hpp \
			Dirac_ops/dirac_DomainWall_4D_eoSolv.hpp \
			Dirac_ops/dirac_overlap_Zolotarev.hpp	\
			Dirac_ops/dirac_Preconditioners.hpp	\
			Dirac_ops/eoUtils.hpp            	\
			Dirac_ops/dirac_Operator_Factory.hpp	\
			EigenModes/eigenModes_IRL.hpp		\
			EigenModes/sortEigen.h			\
			EigenModes/eigenProc_Zolotarev.h	\
			Main/Geometry/siteIndex.h		\
			Main/Geometry/siteIndex_eo.h		\
			Main/Geometry/mapper.hpp		\
			$(MAPPER_HEAD)		     	        \
			Main/gaugeConf.hpp			\
			Main/commandline.cpp			\
			Fields/field_expressions.hpp		\
			HMC/hmcGeneral.hpp                      \
			HMC/mdExec.hpp		                \
			HMC/mdExec_leapfrog.hpp			\
			HMC/mdExec_leapfrog_eigen.hpp		\
			HMC/mdExec_factory.hpp			\
			HMC/mdExec_factory_abs.hpp		\
			Solver/solver_Factory.hpp		\
			Solver/solver_BiCGStab.hpp		\
			Solver/solver.hpp			\
			Solver/solver_CG.hpp			\
			Solver/multiShiftSolver_CG.hpp		\
			Solver/multiShiftSolver.hpp		\
			Smearing/APEsmear.hpp			\
			Smearing/BaseSmear.hpp			\
			Smearing/StoutSmear.hpp			\
			Smearing/SmartConf.hpp			\
			Smearing/Smearing_Factories.hpp		\
			Tools/sunMat.hpp			\
			Tools/sunMatUtils.hpp			\
			Tools/randNum_MP.h			\
			Tools/randNum.h				\
			Tools/sunVec.hpp			\
			Tools/randNum_MT19937.h			\
			Tools/randNum_Factory.h			\
			Tools/RAIIFactory.hpp			\
			Tools/EnumToString.hpp			\
			Tools/byteswap.hpp			\
			Tools/gammaMatrices.hpp			\
			Tools/Remez/alg_remez.hpp		\
			Tools/Remez/bigfloat.hpp





