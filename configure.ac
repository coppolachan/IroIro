#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([IroIro], [0.9.8], [cossu@post.kek.jp])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE([1.8 -Wall -Werror foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])]) ##backward compatibility

AC_CONFIG_HEADERS([include/config.h])
AC_PROG_RANLIB

AC_MSG_NOTICE([

Configuring $PACKAGE_STRING for $host
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
#####################################################################
# identify architecture
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
echo checking architecture... $arch
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

######################################################################
# many HPC systems are configured so that running parallel programs
# interactively is disabled: on those systems, AC_PROG_F77 and AC_PROG_CC
# would fail because they can't run the compiled executables.
# to work around that, let's pretend we are cross-compiling even if we aren't
# !!! this relies on undocumented Autoconf behavior !!!
#
# default flags are also defined

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
    use_cpp_flags="-O3" 
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
    use_cpp_flags="-O3 -I/opt/local/include"
    ld_flags="${LDFLAGS} -L/opt/local/lib"
    use_ar_flags="cru"
    ;;
aix )
     host=$build 
     cross_compiling=yes  # this will do the trick
     mpi=true
     hitachi=true
     AC_MSG_NOTICE([forcing cross compilation on AIX])
     use_cpp="mpCC";
     use_cpp_flags="-O3 -q64 -qlanglvl=stdc99 -qrtti=type  -I${HOME}/gmp/include -I${HOME}/gsl/include -I${HOME}/mpfr/include -qflag=i:s"
     ld_flags="${LDFLAGS} -L${HOME}/gsl/lib -L${HOME}/gmp/lib -L${HOME}/mpfr/lib"
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
     use_cpp_flags="-O3 -qlanglvl=stdc99 -qrtti=type -qhot=level=1 -qsimd=auto -I/bghome/scbadm/ibm-doi/lib -I${HOME}/gmp/include -I${HOME}/gsl/include -I${HOME}/mpfr/include -qflag=i:s"
     ld_flags="${LDFLAGS} -L${HOME}/gsl/lib -L${HOME}/gmp/lib -L${HOME}/mpfr/lib"
     use_ar_flags="cru"
     AC_DEFINE([BIG_ENDIAN_TYPE],[],[Using Big Endian machine type])
     ;;

esac
###################################################################
# Checks for programs.
if test "$cpp" = "" ; then cpp="$use_cpp" ; fi
AC_PROG_CXX($cpp)
cpp=$CXX

# Setting up compiler and linker flags
if test "$cpp_flags" = "" ; then cpp_flags="$use_cpp_flags" ; fi
CXXFLAGS=$cpp_flags
echo configure: setting c++ compiler flags... $CXXFLAGS

LDFLAGS=$ld_flags
echo configure: setting library flags... $LDFLAGS

if test "$ar_flags" = "" ; then ar_flags="$use_ar_flags" ; fi
ARFLAGS=$ar_flags
echo configure: setting linker flags... $ARFLAGS


###################################################################
# Automake flags
AM_CONDITIONAL([USE_MPI], [(test x$mpi = xtrue)])
AM_CONDITIONAL([HITACHI], [(test x$hitachi = xtrue)])
AM_CONDITIONAL([IBM_BGQ], [(test x$bluegeneq = xtrue)])
AM_CONDITIONAL([IBM_BGQWILSON_NOSMP], [(test x$bgq_kernels_nosmp = xtrue)])
AM_CONDITIONAL([IBM_BGQWILSON_SMP], [(test x$bgq_kernels_smp = xtrue)])
AM_CONDITIONAL([USE_EVENODD], [(test x$even_odd_idx = xtrue)])	 

###################################################################
# Checks for libraries and headers
if [ $bgq_kernels_nosmp || $bgq_kernels_smp ]
then
AC_CHECK_HEADERS([bgqwilson.h],,
	[AC_MSG_ERROR(BGQ Wilson library not found in your system.
	Please install or provide the correct path to your installation
	or switch to non optimized libraries.)])
fi

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


AC_CHECK_LIB([gmp],[__gmpf_init],,
	[AC_MSG_ERROR(GNU Multiple Precision GMP library was not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.gmplib.org)])
AC_CHECK_LIB([mpfr],[mpfr_init],,
	[AC_MSG_ERROR(GNU Multiple Precision MPFR library was not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.mpfr.org/)])
	
###################################################################
# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt])

##################################################################
# Define macros for compilation
if [ $mpi ]
then
AC_DEFINE([HAVE_MPI],[],[Processing MPI routines])
fi

if [ $hitachi ]
then
AC_DEFINE([HITACHISR16K],[],[Specific tweaks for Hitachi SR16K])
fi

if [ $bluegeneq ]
then
AC_DEFINE([IBM_BGQ],[],[Specific tweaks for IBM BG/Q])
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

AC_CONFIG_FILES([Makefile
 		lib/Makefile
		tests/Makefile
		docs/Makefile])
		 
AC_OUTPUT

echo "
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary of configuration for $PACKAGE $VERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
echo "- build DOXYGEN documentation   : `if test "x$enable_doc" = xyes; then echo yes; else echo no; fi`
- graphs and diagrams           : `if test "x$enable_dot" = xyes; then echo yes; else echo no; fi`

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
