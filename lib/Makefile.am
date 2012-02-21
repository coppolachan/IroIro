SUBDIRS = . 

if USE_MPI
COMMUNICATOR_SRCS = Communicator/communicator_mpi.cpp
MAPPER_SRCS	  = Main/Geometry/map_parallelscalar.cpp
MAPPER_HEAD	  = Main/Geometry/map_parallelscalar.hpp
if HITACHI	
CXXFLAGS += -q64 -qarch=pwr7 -qtune=pwr7 -qmaxmem=-1 
LDFLAGS += -q64
else
CXXFLAGS += 
LDFLAGS += 
endif
else
COMMUNICATOR_SRCS = Communicator/communicator_dummy.cpp	
MAPPER_SRCS	  = Main/Geometry/map_scalar.cpp
MAPPER_HEAD	  = Main/Geometry/map_scalar.hpp
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
		Main/common_fields.cpp			\
		$(COMMUNICATOR_SRCS)     	        \
		Action/action.cpp			\
	 	Action/action_gauge_wilson.cpp 		\
		Communicator/comm_io.cpp		\
		Communicator/fields_io.cpp		\
		Measurements/GaugeM/staples.cpp		\
		Main/Geometry/siteIndex.cpp		\
		Main/Geometry/siteIndex_eo.cpp		\
		Main/Geometry/mapper.cpp		\
		$(MAPPER_SRCS)		     	        \
		Main/gaugeConf.cpp			\
		Main/commandline.cpp			\
		PugiXML/pugixml.cpp			\
		PugiXML/pugi_interface.cpp		\
		Tools/sunMatUtils.cpp			\
		Tools/randNum_MT19937.cpp		\
		Tools/randNum_Factory.cpp		\
		Tools/gammaMatrices.cpp			


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
			Action/action.hpp				\
			Action/action_gaugetype_factory.hpp		\
			Action/action_Factory.hpp			\
			Measurements/GaugeM/staples.hpp				\
			Communicator/communicator.h		\
			Communicator/comm_io.hpp		\
			Communicator/fields_io.hpp		\
			Main/Geometry/siteIndex.h		\
			Main/Geometry/siteIndex_eo.h		\
			Main/Geometry/mapper.hpp		\
			$(MAPPER_HEAD)		     	        \
			Main/gaugeConf.hpp			\
			Main/commandline.cpp			\
			Fields/field_expressions.hpp		\
			Tools/sunMat.hpp			\
			Tools/sunMatUtils.hpp			\
			Tools/randNum_MP.h			\
			Tools/randNum.h				\
			Tools/sunVec.h				\
			Tools/randNum_MT19937.h			\
			Tools/randNum_Factory.h			\
			Tools/RAIIFactory.hpp			\
			Tools/EnumToString.hpp			\
			Tools/byteswap.hpp			\
			Tools/gammaMatrices.hpp			





