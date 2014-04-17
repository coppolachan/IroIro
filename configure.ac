#                         -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
# Project IROIRO JLQCD package (KEK) 
# Lattice QCD MonteCarlo simulations
#
# Time-stamp: <2014-04-17 12:56:03 cossu>

AC_PREREQ([2.59])
AC_INIT([IroIro], [1.0.7], [cossu@post.kek.jp])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE([1.10 foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])]) ##backward compatibility

AC_CONFIG_HEADERS([include/config.h])
AC_PROG_RANLIB


AC_MSG_NOTICE([

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Configuring $PACKAGE_STRING for $host
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
])

# Declares code language
AC_LANG(C++)

#####################################################################
# Save some environment variables
cpp=$CXX
cpp_flags=$CXXFLAGS
ld_flags=$LDFLAGS
ar_flags=$ARFLAGS
hitachi=false  # default value 
bluegeneq=false   # default value
SR_INCLUDE=/srhome/.tools/include
SR_LIBS=/srhome/.tools/lib

BLUEGENE_INCLUDE=/bghome/scbadm/ibm-doi/lib
BLUEGENE_LIBS=/bghome/scbadm/ibm-doi/lib
#BLUEGENE_LIBS=/usr/local/bg/lib
#BLUEGENE_INCLUDE=/usr/local/bg/include
#####################################################################
# identify architecture
echo           Checking architecture
echo :::::::::::::::::::::::::::::::::::::::::::
if test "$arch" = ""
then
        case $host in
        ia64-*-linux-gnu )      arch=ia64     ;;
        x86_64-*-linux-gnu )    arch=amd64    ;;
        *-pc-linux-gnu )        arch=ia32     ;;
        *-ibm-aix* )            arch=aix      ;;
	powerpc64-* )           arch=bgq      ;;
        powerpc-apple-darwin* ) arch=macppc   ;;
        i*86-apple-darwin* )    arch=macx86   ;;
	x86_64-apple-darwin* )  arch=mac86_64 ;;
	* )
	AC_MSG_WARN(unsupported architecture. Setting defaults...)
	use_cpp="g++"
	use_cpp_flags="-g"
	;;
esac
fi
echo checking architecture result... $arch
###################################################################
# MPI communicator flags
AC_ARG_ENABLE(mpi,
   AC_HELP_STRING([--enable-mpi],[Enable MPI (parallelization) support routines.]),
   [case "${enableval}" in
     yes) mpi=true  ;;
     no)  mpi=false ;;
     *) AC_MSG_ERROR([bad value ${enableval} for --enable-mpi . MPI not used.]) ;;
   esac],[mpi=false])

###################################################################
# Enabling Improved version routines
AC_ARG_ENABLE(improved,
	AC_HELP_STRING([--enable-improved], [Enable compilation of improved version of kernels.]),
	[case "${enableval}" in
     	yes) improved_kernels=true  ;;
     	no)  improved_kernels=false ;;
     	*) AC_MSG_ERROR([bad value ${enableval} for --enable-improved . Compiling with standard kernels.]) ;;
   esac],[improved_kernels=false])

###################################################################
# Enabling BlueGene Wilson improved version routines single thread
AC_ARG_ENABLE(bgq-improved,
	AC_HELP_STRING([--enable-bgq-improved], [Enable compilation of BlueGene/Q improved version of kernels. Single threaded version.]),
	[case "${enableval}" in
     	yes) bgq_kernels_nosmp=true  ;;
     	no)  bgq_kernels_nosmp=false ;;
     	*) AC_MSG_ERROR([bad value ${enableval} for --enable-bgq-improved . Compiling with standard kernels.]) ;;
   esac],[bgq_kernels_nosmp=false])

###################################################################
# Enabling Hitachi SR16K Wilson improved version routines single thread
AC_ARG_ENABLE(sr16k-improved,
	AC_HELP_STRING([--enable-sr16k-improved], [Enable compilation of Hitachi SR16K improved version of kernels.]),
	[case "${enableval}" in
     	yes) sr16k_kernels=true  ;;
     	no)  sr16k_kernels=false ;;
     	*) AC_MSG_ERROR([bad value ${enableval} for --enable-sr16k-improved . Compiling with standard kernels.]) ;;
   esac],[sr16k_kernels=false])

###################################################################
# Enabling BlueGene Wilson improved version routines multi thread
AC_ARG_ENABLE(bgq-improved-smp,
	AC_HELP_STRING([--enable-bgq-improved-smp], [Enable compilation of BlueGene/Q improved version of kernels. Multithreaded SMP version.]),
	[case "${enableval}" in
     	yes) bgq_kernels_smp=true  ;;
     	no)  bgq_kernels_smp=false ;;
     	*) AC_MSG_ERROR([bad value ${enableval} for --enable-bgq-improved-smp . Compiling with standard kernels.]) ;;
   esac],[bgq_kernels_smp=false])

###################################################################
# Enabling Even Odd indexing 
AC_ARG_ENABLE(evenodd-idx,
        AC_HELP_STRING([--enable-evenodd-idx], [Enable even odd checkerboard indexing for the gauge field.]),
        [case "${enableval}" in
        yes) even_odd_idx=true  ;;
        no)  even_odd_idx=false ;;
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-evenodd-idx . Using lexicographic indexing.]) ;;
   esac],[even_odd_idx=false])


###################################################################
# Enabling DCMT random number generator linking 
AC_ARG_WITH(dcmt,
        AC_HELP_STRING([--with-dcmt=DCMT_BASEDIR], [Enable Dynamic Creation Mersenne Twister RNG linkage. If DCMT_BASEDIR is not specified it assumes that libraries and header are somewhere in the path]),
        [ dcmt_link=true
	  dcmt_dir=${with_dcmt}
	  dcmt_include=${with_dcmt}/include
	  dcmt_lib=${with_dcmt}/lib	
	],[dcmt_link=false])

###################################################################
# Enabling Even Odd indexing BGQ fast optimizations
AC_ARG_ENABLE(bgq-fast,
        AC_HELP_STRING([--enable-bgq-fast], [Enable automatically all optimizations for BGQ. It is equivalent to switching on --enable-evenodd-idx --enable-bgq-improved-smp --enable-improved -with-dcmt=<your home>/dcmt]),
        [case "${enableval}" in
        yes) even_odd_idx=true
	     bgq_kernels_smp=true
	     improved_kernels=true	
	     dcmt_link=true
	     dcmt_dir=${HOME}/dcmt
	     dcmt_include=${dcmt_dir}/include
	     dcmt_lib=${dcmt_dir}/lib
	;;  
        no)  
	;;   
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-bgq-fast. Using basic settings.]) ;;
   	esac],[])

###################################################################
# Enabling BFM (P.Boyle) Bagel Fermion Matrix classes linking 
AC_ARG_WITH(bfm,
        AC_HELP_STRING([--with-bfm=BFM_BASEDIR], [Enable BFM Bagel Fermion Matrix libraries (P.Boyle) linkage (Only on BGQ). If BFM_BASEDIR is not specified it assumes that libraries and header are somewhere in the path]),
        [ bfm_link=true
	  bfm_dir=${with_bfm}
	  bfm_include=${with_bfm}/include
	  bfm_lib=${with_bfm}/lib	
	],[bfm_link=false])


######################################################################
# many HPC systems are configured so that running parallel programs
# interactively is disabled: on those systems, AC_PROG_F77 and AC_PROG_CC
# would fail because they can't run the compiled executables.
# to work around that, let's pretend we are cross-compiling even if we aren't
# !!! this relies on undocumented Autoconf behavior !!!
#
# default flags are also defined

echo
echo           Setting up flags
echo :::::::::::::::::::::::::::::::::::::::::::

case "$arch" in
ia64 | amd64 | ia32 )
   if [ $mpi ]
    then
      # list of tested MPI compilers	
      use_cpp="mpicxx.openmpi mpicxx mpic++" 
    else
     # list of tested compilers
      use_cpp="g++ icpc"; 
    fi
    use_cpp_flags="-O3 -fopenmp -I${HOME}/lime/include" 
    ld_flags="${LDFLAGS} -L${HOME}/lime/lib"       
    use_ar_flags="cru"
    ;;
macppc | macx86 | mac86_64 )
   if [ $mpi ]
    then
      # list of tested MPI compilers	
      use_cpp="mpicxx.openmpi mpicxx mpic++" 
    else
     # list of tested compilers
      use_cpp="icpc g++"; 
    fi
    use_cpp_flags="-O3 -I/opt/local/include -I${HOME}/lime/include"
    ld_flags="${LDFLAGS} -L/opt/local/lib -L${HOME}/lime/lib"
    use_ar_flags="cru"
    ;;
aix )
     host=$build 
     cross_compiling=yes  # this will do the trick
     mpi=true
     hitachi=true
     AC_MSG_NOTICE([forcing cross compilation on AIX])
     use_cpp="mpCC";
     use_cpp_flags="-O3 -q64 -qlanglvl=stdc99 -qrtti=type -I${SR_INCLUDE} -I${HOME}/gmp/include -I${HOME}/gsl/include -I${HOME}/mpfr/include -I${HOME}/lime/include -qflag=i:s"
     ld_flags="${LDFLAGS} -L${HOME}/gsl/lib -L${HOME}/gmp/lib -L${HOME}/mpfr/lib -L${HOME}/lime/lib -L${SR_LIBS}"
     use_ar_flags="-X64 -c -ru"
     AC_DEFINE([BIG_ENDIAN_TYPE],[],[Using Big Endian machine type])
     ;;
bgq ) 
     host=$build
     cross_compiling=yes  # this will do the trick     
     mpi=true
     bluegeneq=true
     AC_MSG_NOTICE([forcing cross compilation on BG/Q])
     if [ $bgq_kernels_smp ]
     then 
     	use_cpp="bgxlC_r";
     else
	use_cpp="bgxlC";
     fi
     #use_cpp_flags="-O3 -qlanglvl=stdc99 -qrtti=type -qhot=level=1 -qsimd=auto -q64 -qarch=qp -qtune=qp -I${BLUEGENE_INCLUDE} -I${HOME}/gmp/include -I${HOME}/gsl/include -I${HOME}/mpfr/include -I${HOME}/lime/include -qflag=i:s"
     use_cpp_flags="-O3 -qlanglvl=stdc99 -qrtti=type -qhot=level=1 -qsimd=auto -q64 -qarch=qp -qtune=qp -I${BLUEGENE_INCLUDE} -I${HOME}/gmp/include -I${HOME}/gsl/include -I${HOME}/mpfr/include -I${HOME}/lime/include -qflag=i:s" 
     ld_flags="${LDFLAGS} -L${HOME}/gsl/lib -L${HOME}/gmp/lib -L${HOME}/mpfr/lib -L${HOME}/lime/lib"
     use_ar_flags="cru"
     AC_DEFINE([BIG_ENDIAN_TYPE],[],[Using Big Endian machine type])
     ;;

esac

###################################################################
# Automake flags
AM_CONDITIONAL([USE_MPI], [(test x$mpi = xtrue)])
AM_CONDITIONAL([HITACHI], [(test x$hitachi = xtrue)])
AM_CONDITIONAL([IBM_BGQ], [(test x$bluegeneq = xtrue)])
AM_CONDITIONAL([IBM_BGQWILSON_NOSMP], [(test x$bgq_kernels_nosmp = xtrue)])
AM_CONDITIONAL([IBM_BGQWILSON_SMP], [(test x$bgq_kernels_smp = xtrue)])
AM_CONDITIONAL([SR_WILSON], [(test x$sr16k_kernels = xtrue)])
AM_CONDITIONAL([IBM_BGQKERNELS], [[ ($bgq_kernels_nosmp == true) || ($bgq_kernels_smp == true)]])
AM_CONDITIONAL([USE_EVENODD], [(test x$even_odd_idx = xtrue)])	 
AM_CONDITIONAL([USE_BFM], [(test x$bfm_link = xtrue)&&(test x$bluegeneq = xtrue)])	 
AM_CONDITIONAL([USE_DCMT], [(test x$dcmt_link = xtrue)])	 

###################################################################
# Checks for programs.
echo
echo        Checking compiler features
echo :::::::::::::::::::::::::::::::::::::::::::

if test "$cpp" = "" ; then cpp="$use_cpp" ; fi
AC_PROG_CXX($cpp)
cpp=$CXX

# Setting up global compiler and linker flags
if test "$cpp_flags" = "" ; then cpp_flags="$use_cpp_flags" ; fi
CXXFLAGS=$cpp_flags
CPPFLAGS=$cpp_flags # preprocessof flags
AM_COND_IF([USE_BFM], [CXXFLAGS+=" -I$bfm_include"],,)
AM_COND_IF([USE_DCMT], [CXXFLAGS+=" -I$dcmt_include"],,)
AM_COND_IF([IBM_BGQWILSON_SMP], [CXXFLAGS+=" -qsmp=omp:auto"])
AM_COND_IF([IBM_BGQ], [DRIVER="/bgsys/drivers/ppcfloor"
		      CXXFLAGS+=" -I$DRIVER -I$DRIVER/spi/include/kernel/cnk -I$DRIVER/spi/include/mu/default -I$DRIVER/comm/xl/include"],,)
echo configure: setting c++ compiler flags... done

LDFLAGS=$ld_flags
AM_COND_IF([USE_BFM], [LDFLAGS+=" -L$bfm_lib"],,)
AM_COND_IF([USE_DCMT], [LDFLAGS+=" -L$dcmt_lib"],,)
#AM_COND_IF([IBM_BGQ], [LDFLAGS+=" -L$DRIVER/comm/lib -L$DRIVER/comm/xl/lib -L$DRIVER/comm/sys/lib -L$DRIVER/spi/lib -L$DRIVER/bgpm/lib -L$BLUEGENE_LIBS"],,)
AM_COND_IF([IBM_BGQ], [LDFLAGS+=" -L$DRIVER/comm/lib -L$DRIVER/spi/lib -L$DRIVER/bgpm/lib -L$BLUEGENE_LIBS"],,)
AM_COND_IF([SR_WILSON], [LDFLAGS+=" -L$SR_LIBS"],,)
echo configure: setting library flags... done

if test "$ar_flags" = "" ; then ar_flags="$use_ar_flags" ; fi
ARFLAGS=$ar_flags
echo configure: setting linker flags... done



###################################################################
# Checks for libraries and headers
echo
echo Checking libraries 
echo :::::::::::::::::::::::::::::::::::::::::::

AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([gslcblas],
	[cblas_dgemm],,
	[AC_MSG_ERROR(GNU Scientific Library BLAS not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.gnu.org/s/gsl/)])

AC_CHECK_LIB([gsl],[gsl_sf_elljac_e],,
	[AC_MSG_ERROR(GNU Scientific Library not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.gnu.org/s/gsl/)])
AC_CHECK_FUNC([gsl_sf_ellint_Kcomp],,
	[AC_MSG_ERROR(Missing function in GNU Scientific Library)])

AC_CHECK_HEADER(gsl/gsl_sf_ellint.h,,
	[AC_MSG_ERROR(GNU Scientific Library headers not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.gnu.org/s/gsl/)])

AC_CHECK_HEADER(gsl/gsl_sf_elljac.h,,
	[AC_MSG_ERROR(GNU Scientific Library headers not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.gnu.org/s/gsl/)])

AC_CHECK_LIB([gmp],[__gmpf_init],,
	[AC_MSG_ERROR(GNU Multiple Precision GMP library was not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.gmplib.org)])

AC_CHECK_LIB([mpfr],[mpfr_init],,
	[AC_MSG_ERROR(GNU Multiple Precision MPFR library was not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.mpfr.org/)])

AC_CHECK_LIB([lime],[limeCreateReader],,
	[AC_MSG_ERROR(C-LIME library was not found in your system.
Please install or provide the correct path to your installation [default search path ~/lime/]
Info at: http://usqcd.jlab.org/usqcd-docs/c-lime/)])

if [ ($dcmt_link) ]
then
AC_CHECK_LIB([dcmt], [global_mt19937],,  
	[AC_MSG_ERROR(Dynamic Creation Mersenne Twister library was not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/DC/dc.html)]) 
fi

# Hitashi SR16K specific stuff
if [ $sr16k_kernels ]
then
echo :::: Hitachi SR16K related stuff
AC_CHECK_LIB([srmwilson],[SRWilsonLA_SetConst],,
	[AC_MSG_ERROR(SRWilson library not found in your system or missing dependences.
	Please install or provide the correct path to your installation
	or switch off SRWilson support.)],[])
echo enabled SRWilson support on Hitachi SR16K...
fi

# BGQ specific stuff
if [ $bluegeneq ]
then
echo :::: BGQ related stuff
AC_CHECK_LIB([bgnet],[BGNET_Init],,
	[AC_MSG_ERROR(BGNET library not found in your system or missing dependences.
	Please install or provide the correct path to your installation
	or switch off BGNET support.)],[-lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++])
echo enabled BGNET support on BGQ...


if [ $bgq_kernels_smp ]
then
#AC_CHECK_LIB([bgqwilsonsmp],[BGWilson_Init],,
#	[AC_MSG_ERROR(BGWilsonSMP library not found in your system or missing dependences.
#	Please install or provide the correct path to your installation
#	or switch off BGNET support.)],[-lbgnet -lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++])
AC_CHECK_LIB([bgqwilson3x4smp],[BGWilson_Init],,
	[AC_MSG_ERROR(BGWilson3x4SMP library not found in your system or missing dependences.
	Please install or provide the correct path to your installation
	or switch off BGNET support.)],[-lbgnet -lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++])
echo enabled BGWilson Multiprocessor SMP version support on BGQ...
fi
if [ $bgq_kernels_nosmp ]
then
#AC_CHECK_LIB([bgqwilson],[BGWilson_Init],,
#	[AC_MSG_ERROR(BGWilson library not found in your system or missing dependences.
#	Please install or provide the correct path to your installation
#	or switch off BGNET support.)],[-lbgnet -lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++])
AC_CHECK_LIB([bgqwilson3x4],[BGWilson_Init],,
	[AC_MSG_ERROR(BGWilson3x4 library not found in your system or missing dependences.
	Please install or provide the correct path to your installation
	or switch off BGNET support.)],[-lbgnet -lSPI_cnk -lSPI -lSPI_l1p -lrt -lbgpm -lstdc++])
echo enabled BGWilson single threaded version support on BGQ...
fi

fi #bluegeneq

# BFM Bagel Fermion Class Peter Boyle's library
if [ ($bfm_link) ]
then
if [ $bluegeneq ]
then
AC_CHECK_LIB([bfm],[vmx_vaxpy],,
	[AC_MSG_ERROR(BFM library not found in your system.
	Please install or provide the correct path to your installation
	or switch off BFM support.)])
echo enabled BFM Bagel Fermion Matrix support on BGQ...
else
echo WARN: BFM support enabled only on BGQ, switching off...
bfm_link=false
fi
fi	
###################################################################
# Checks for header files.
echo
echo Checking headers 
echo :::::::::::::::::::::::::::::::::::::::::::

AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt])


# Checks for SRWilson specific headers 
if [ $sr16k_kernels ]
then
AC_CHECK_HEADERS([srmwilson.h],,
	[AC_MSG_ERROR(SRWilson header not found in your system.
	Please install or provide the correct path to your installation
	or switch to non optimized libraries.)])
fi

# Checks for BGQ specific headers 
if [ $bgq_kernels_nosmp || $bgq_kernels_smp ]
then
# BGQ IBM Wilson operator 
AC_CHECK_HEADERS([bgqwilson.h],,
	[AC_MSG_ERROR(BGQ Wilson header not found in your system.
	Please install or provide the correct path to your installation
	or switch to non optimized libraries.)])
# BGQ IBM Communication library
AC_CHECK_HEADERS([bgnet.h],,
	[AC_MSG_ERROR(BGNet header not found in your system.
	Please install or provide the correct path to your installation
	or switch to non optimized libraries.)])
fi

if [ $bgq_kernels_smp ]
then
AC_CHECK_HEADERS([bgqthread.h],,
	[AC_MSG_ERROR(BGQThread header not found in your system.
	Please install or provide the correct path to your installation
	or switch to non optimized libraries.)])
fi

##################################################################
# Define macros for compilation
if [ $mpi ]
then
AC_DEFINE([HAVE_MPI],[],[Processing MPI routines])
fi

if [ $hitachi ]
then
AC_DEFINE([HITACHISR16K],[],[Specific tweaks for Hitachi SR16K])

if [ $sr16k_kernels ]
then
AC_DEFINE([SR16K_WILSON],[],[Enable Hitachi SR16K optimized Wilson kernel])
fi
fi #hitachi


if [ $bluegeneq ]
then
AC_DEFINE([IBM_BGQ],[],[Specific tweaks for IBM BG/Q])
fi


if [ $even_odd_idx ]
then
AC_DEFINE([USING_EVENODD_IDX],[],[Using even odd storage indexing for the Gauge field])
fi

# Turn on all improved kernels
if [ $improved_kernels ]
then 
AC_DEFINE([IMPROVED_WILSON],[], [Enable improved Wilson kernel])
fi

# BGQ specifics - improved kernels
if [ $bgq_kernels_nosmp || $bgq_kernels_smp ]
then 
AC_DEFINE([IBM_BGQ_WILSON],[], [Enable BlueGene/Q optimized Wilson kernel])
fi
###################################################################
# Checks for doxygen support
# if present enables the "make doxyfile" command
echo
echo Checking doxygen support 
echo :::::::::::::::::::::::::::::::::::::::::::
AC_PROG_DOXYGEN

if test -n "$DOXYGEN"
then
AC_CONFIG_FILES([docs/doxy.cfg])
fi

# Verbosity level - default level is 1
AC_ARG_VAR(code_verbosity, Sets code verbosity level 0..5)
AC_DEFINE_UNQUOTED([VERBOSITY],[${code_verbosity=1}],[Code Verbosity level])
###################################################################
AC_SUBST(ARFLAGS)
AC_SUBST(DRIVER)

AC_CONFIG_FILES([Makefile
 		lib/Makefile
		tests/Makefile
		docs/Makefile])

echo
echo Creating makefiles
echo :::::::::::::::::::::::::::::::::::::::::::		 
AC_OUTPUT

echo "
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary of configuration for $PACKAGE $VERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following features are enabled:

- architecture                  : $arch
- cross compilation             : `if test "x$cross_compiling" = xyes; then echo yes; else echo no; fi`
- gauge field indexing          : `if test "x$even_odd_idx" = xtrue; then echo even-odd; else echo lexicographic; fi`
- multicore MPI version enabled : $mpi"
if [ $bgq_kernels_nosmp || $bgq_kernels_smp ]
then 
echo -ne "- IBM kernels compilation       : "; if test "x$bgq_kernels_nosmp" = xtrue; then echo single_threaded; else echo smp; fi
else
echo "- improved kernels compilation  : $improved_kernels"
fi
if [ $hitachi ]
then 
echo -ne "- SR kernels compilation        : "; if test "x$sr16k_kernels" = xtrue; then echo enabled; else echo not enabled; fi
fi

echo "- build DOXYGEN documentation   : `if test "x$enable_doc" = xyes; then echo yes; else echo no; fi`
- graphs and diagrams           : `if test "x$enable_dot" = xyes; then echo yes; else echo no; fi`
- optional BFM library linking  : `if test "x$bfm_link" = xtrue; then echo yes; else echo no; fi`
- optional DCMT Random numbers  : `if test "x$dcmt_link" = xtrue; then echo yes; else echo no; fi`

Directories:

- Source code location          : $srcdir
- Install path                  : $prefix

Compiler flags setup:
-------------------------------
CXX       (Compiler)          : $CXX
CXXFLAGS  (Compiler flags)    : $CXXFLAGS
LDFLAGS	  (Library flags)     : $LDFLAGS
ARFLAGS	  (Linker flags)      : $ARFLAGS
-------------------------------

The configure script has generated makefiles for your platform.
You have several options now:

make               Depends on configure's flags given, builds the docs,
                   builds $PACKAGE version $VERSION
make doxyfile      Makes the docs using doxygen (if enabled)
make dist	   Creates a tar.gz file to be distributed
make install       Install code in $prefix
make uninstall     Uninstall IroIro.
"
