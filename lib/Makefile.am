SUBDIRS = . 

CXXFLAGS +=
LDFLAGS +=
OPTLIB =

DWF_BGQ_SRC = 

if USE_MPI
COMMUNICATOR_SRCS = Communicator/communicator_mpi.cpp

if HITACHI	
CXXFLAGS += -q64 -qarch=pwr7 -qtune=pwr7 -qmaxmem=-1 
LDFLAGS += -q64
endif

if IBM_BGQ
DRIVER=/bgsys/drivers/ppcfloor
CXXFLAGS += -q64 -qarch=qp -qtune=qp -I$(DRIVER) \
	-I$(DRIVER)/spi/include/kernel/cnk \
	-I$(DRIVER)/spi/include/mu/default \
	-I$(DRIVER)/comm/xl/include	

if IBM_BGQWILSON_NOSMP  
DWF_BGQ_SRC += Dirac_ops/dirac_DomainWall_BGQ.cpp
endif

if IBM_BGQWILSON_SMP
DWF_BGQ_SRC += Dirac_ops/dirac_DomainWall_BGQ.cpp
endif

if IBM_BGQWILSON_SMP
CXXFLAGS += -qsmp=omp:auto
endif

LDFLAGS += -L$(DRIVER)/comm/xl/lib -L$(DRIVER)/comm/sys/lib \
	 -L$(DRIVER)/spi/lib -L$(DRIVER)/bgpm/lib \
	 -L/bghome/scbadm/ibm-doi/lib
OPTLIB += -lkek -lbgnet -lpthread -lmpich -lopa -lmpl -lpami -lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++ 

if IBM_BGQWILSON_SMP
OPTLIB += -lbgqwilson3x4smp -lbgqthread
else
OPTLIB += -lbgqwilson3x4
endif

endif

else
COMMUNICATOR_SRCS = Communicator/communicator_dummy.cpp	
endif


if USE_EVENODD
SITE_INDEXING_SRCS = Main/Geometry/siteIndexEO.cpp
else
SITE_INDEXING_SRCS = Main/Geometry/siteIndex.cpp
endif


AM_CPPFLAGS=-I$(top_srcdir)
AM_CFLAGS  =-I$(top_srcdir)
AM_LDFLAGS =


lib_LIBRARIES = libcommcode.a

EXTRA_DIST =    HMC/mdExec_leapfrog_eigen.cpp	\
		Main/documentation_pages.h	\
                Dirac_ops/dirac_wilson_improved.code    \
                Dirac_ops/dirac_wilson_standard.code    \
		Fields/field.inl

libcommcode_a_SOURCES = Main/commonPrms.cpp		\
		$(COMMUNICATOR_SRCS)     	        \
		$(SITE_INDEXING_SRCS)			\
		Action/action.cpp			\
	 	Action/action_gauge_wilson.cpp 		\
		Action/action_gauge_rect.cpp            \
		Action/action_gaugetype_factory.cpp     \
		Action/action_Nf2.cpp			\
		Action/action_Nf2_ratio.cpp		\
		Action/action_Nf.cpp			\
		Action/action_Nf_ratio.cpp		\
		Action/action_staggered.cpp		\
		Action/action_staggered_ratio.cpp	\
		Action/action_Nf2_DomainWall.cpp	\
		Action/action_fermiontype_factory.cpp	\
		Communicator/comm_io.cpp		\
		Dirac_ops/dirac.cpp     		\
		Dirac_ops/dirac_wilson.cpp		\
		Dirac_ops/dirac_wilson_EvenOdd.cpp	\
		Dirac_ops/dirac_wilson_Brillouin.cpp	\
		Dirac_ops/dirac_staggered_EvenOdd.cpp	\
		Dirac_ops/dirac_clover.cpp		\
		Dirac_ops/fopr_signH_Zolotarev.cpp	\
		Dirac_ops/dirac_overlap_Zolotarev.cpp	\
		Dirac_ops/dirac_DomainWall.cpp   	\
		$(DWF_BGQ_SRC)				\
		Dirac_ops/dirac_DomainWall_EvenOdd.cpp  \
		Dirac_ops/dirac_DomainWall_4D_fullSolv.cpp \
		Dirac_ops/dirac_DomainWall_4D_eoSolv.cpp \
		Dirac_ops/eoUtils.cpp           	\
		Dirac_ops/dirac_Operator_Factory.cpp	\
		Dirac_ops/boundaryCond.cpp      	\
		EigenModes/eigenModesSolver_Factory.cpp	\
		EigenModes/eigenModesSolver_IRL.cpp	\
		EigenModes/eigenProc_Zolotarev.cpp	\
		EigenModes/eigenSorter.cpp		\
		EigenModes/fopr_chebyshev_DdagD.cpp	\
		EigenModes/foprHermFactory.cpp   	\
		EigenModes/findminmax.cpp		\
		Measurements/measGeneral.cpp		\
		Measurements/GaugeM/staples.cpp		\
		Measurements/GaugeM/wilsonFlow.cpp	\
		Measurements/GaugeM/polyakovLoop.cpp			\
		Measurements/GaugeM/gaugeFixingStep.cpp			\
		Measurements/GaugeM/gaugeFixing_Coulomb.cpp		\
		Measurements/GaugeM/gaugeFixing_Landau.cpp		\
		Measurements/GaugeM/gfixFactory.cpp			\
		Measurements/FermionicM/qprop_MultiShift.cpp 		\
		Measurements/FermionicM/qprop.cpp 			\
		Measurements/FermionicM/qprop_mom.cpp 			\
		Measurements/FermionicM/qprop_EvenOdd.cpp		\
		Measurements/FermionicM/qprop_DomainWall.cpp		\
		Measurements/FermionicM/utils_DWF4d.cpp          	\
		Measurements/FermionicM/quark_prop_meas_factory.cpp 	\
		Main/Geometry/siteIndex_EvenOdd.cpp	\
		Main/Geometry/siteMap.cpp		\
		Main/Geometry/shiftField.cpp		\
		Main/gaugeConf.cpp			\
		Main/gaugeGlobal.cpp			\
		Main/commandline.cpp			\
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
		Solver/rationalSolver.cpp		\
		Solver/solver.cpp			\
		Smearing/APEsmear.cpp			\
		Smearing/stoutSmear.cpp			\
		Smearing/smartConf.cpp			\
		Smearing/smearingFactories.cpp		\
		Tools/fieldUtils.cpp			\
		Tools/jobUtils.cpp			\
		Tools/sunMatUtils.cpp			\
		Tools/sunAdjMatUtils.cpp		\
		Tools/randNum_MT19937.cpp		\
		Tools/randNum_Factory.cpp		\
		Tools/gammaMatrices.cpp			\
		Tools/Remez/alg_remez.cpp		\
		Tools/RationalApprox/rationalapprox.cpp	


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
			$(top_srcdir)/include/foprHermFactory.hpp	\
			$(top_srcdir)/include/format_A.h		\
			$(top_srcdir)/include/format_F.h		\
			$(top_srcdir)/include/format_G.h		\
			$(top_srcdir)/include/format_S.h		\
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
			Measurements/FermionicM/fermion_meas_factory_abs.hpp \
			Measurements/FermionicM/quark_prop_meas_factory.hpp \
                        Measurements/FermionicM/utils_DWF4d.hpp \
			Communicator/communicator.h		\
			Communicator/comm_io.hpp		\
			Communicator/fields_io.hpp		\
			Main/Geometry/siteIndex.hpp     	\
			Main/Geometry/siteIndex_EvenOdd.hpp	\
			Main/Geometry/siteMap.hpp		\
			Main/Geometry/mapping.hpp		\
			Main/Geometry/shiftField.hpp		\
			Main/gaugeConf.hpp			\
			Main/gaugeGlobal.hpp			\
			Dirac_ops/dirac_Operator_Factory.hpp	\
			Dirac_ops/dirac_overlap_Zolotarev.hpp	\
			Dirac_ops/dirac.hpp			\
			Dirac_ops/dirac_wilson.hpp		\
			Dirac_ops/dirac_wilson_EvenOdd.hpp	\
			Dirac_ops/dirac_wilson_Brillouin.hpp	\
			Dirac_ops/dirac_staggered_EvenOdd.hpp	\
			Dirac_ops/dirac_clover.hpp		\
			Dirac_ops/dirac_DomainWall.hpp   	\
			Dirac_ops/dirac_DomainWall_EvenOdd.hpp	\
			Dirac_ops/dirac_DomainWall_4D_fullSolv.hpp \
			Dirac_ops/dirac_DomainWall_4D_eoSolv.hpp \
			Dirac_ops/dirac_Preconditioners.hpp	\
			Dirac_ops/eoUtils.hpp            	\
			Dirac_ops/boundaryCond.hpp      	\
			EigenModes/eigenModesSolver_Factory.hpp	\
			EigenModes/eigenModesSolver.hpp  	\
			EigenModes/eigenModesSolver_IRL.hpp	\
			EigenModes/eigenSorter.hpp		\
			EigenModes/eigenProc_Zolotarev.h	\
			EigenModes/findminmax.hpp		\
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
			Solver/rationalSolver.hpp		\
			Smearing/APEsmear.hpp			\
			Smearing/baseSmear.hpp			\
			Smearing/stoutSmear.hpp			\
			Smearing/smartConf.hpp			\
			Smearing/smearingFactories.hpp		\
			Tools/jobUtils.hpp			\
			Tools/sunMat.hpp			\
			Tools/sunMatUtils.hpp   		\
			Tools/sunAdjMat.hpp			\
			Tools/sunAdjMatUtils.hpp		\
			Tools/sunVec.hpp			\
			Tools/randNum_MP.h			\
			Tools/randNum.h				\
			Tools/randNum_MT19937.h			\
			Tools/randNum_Factory.h			\
			Tools/fieldUtils.hpp			\
			Tools/RAIIFactory.hpp			\
			Tools/EnumToString.hpp			\
			Tools/byteswap.hpp			\
			Tools/gammaMatrices.hpp			\
			Tools/Remez/alg_remez.hpp		\
			Tools/Remez/bigfloat.hpp		\
			Tools/RationalApprox/rationalapprox.hpp	





