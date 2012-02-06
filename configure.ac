#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([IroIro], [0.9], [cossu@post.kek.jp])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE([1.8 -Wall -Werror foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES]) ##backward compatibility
AM_SILENT_RULES([yes]) ##compile silently

AC_CONFIG_HEADERS([include/config.h])
AC_PROG_RANLIB

AC_MSG_NOTICE([

Configuring $PACKAGE_STRING for $host
---------------------------------------------
])

# Declares code language
AC_LANG(C++)

#####################################################################
# Save some environment variables
cpp=$CPP
cpp_flags=$CPPFLAGS
ld_flags=$LDFLAGS
ar_flags=$ARFLAGS
#####################################################################
# identify architecture
if test "$arch" = ""
then
        case $host in
        ia64-*-linux-gnu )      arch=ia64    ;;
        x86_64-*-linux-gnu )    arch=amd64   ;;
        *-pc-linux-gnu )        arch=ia32    ;;
        *-ibm-aix* )            arch=aix     ;;
        powerpc-apple-darwin* ) arch=macppc  ;;
        i686-apple-darwin* )    arch=mac686  ;;
	* )
	AC_MSG_WARN(unsupported architecture. Setting defaults...)
	cpp="g++"
	cpp_flags="-g"
	ld_flags=""
	;;
esac
fi
echo checking architecture... $arch
###################################################################
# MPI communicator flags
     AC_ARG_ENABLE([mpi],
     [  --enable-mpi    Turn on mpi communicator routines],
     [case "${enableval}" in
       yes) mpi=true ;;
       no)  mpi=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-mpi]) ;;
     esac],[mpi=false])


     AC_ARG_ENABLE([mpi-hitachi],
     [  --enable-mpi-hitachi   Turn on mpi communicator routines on Hitachi SR16K],
     [case "${enableval}" in
       yes) mpihitachi=true ;;
       no)  mpihitachi=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-mpi-hitachi]) ;;
     esac],[mpihitachi=false])
     AM_CONDITIONAL([USE_MPI], [(test x$mpihitachi = xtrue)||(test x$mpi = xtrue)])
     AM_CONDITIONAL([HITACHI], [test x$mpihitachi = xtrue])
######################################################################
# many HPC systems are configured so that running parallel programs
# interactively is disabled: on those systems, AC_PROG_F77 and AC_PROG_CC
# would fail because they can't run the compiled executables.
# to work around that, let's pretend we are cross-compiling even if we aren't
# !!! this relies on undocumented Autoconf behavior !!!
#
# default flags are also defined

case "$arch" in
ia64 | amd64 | ia32 | macppc | mac686 )
    if [ $mpi ]
    then
      use_cpp="mpicxx mpic++" # list of possible MPI compilers
    else
      use_cpp="icpc g++"; 
    fi
    use_cpp_flags="-O3"
    use_ar_flags="cru"
    ;;
aix )
     host=$build ;
     cross_compiling=yes ; # this will do the trick
     echo forcing cross compilation on AIX
     use_cpp="mpCC";
     use_cpp_flags="-O3 -q64 -qlanglvl=stdc99 -qrtti=type -I${HOME}/gsl/include -qflag=i:s"
     ld_flags="${LDFLAGS} -L${HOME}/gsl/lib"
     use_ar_flags="-X64 -c -ru"
     AC_DEFINE([BIG_ENDIAN],[],[Using Big Endian machine type])
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
# Checks for libraries.
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

###################################################################
# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt])

##################################################################
# Make decisions on MPI flags
if [ $mpi || $mpihitachi ]
then
AC_CHECK_HEADERS([/usr/include/mpich2/mpi.h])
AC_DEFINE([HAVE_MPI],[],[Processing MPI routines])
fi

if [ $mpihitachi ]
then
AC_DEFINE([HITACHISR16K],[],[Specific tweaks for Hitachi SR16K])
fi
###################################################################
# Checks for doxygen support
# if present enables the "make doxyfile" command
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

AM_CONDITIONAL([HAVE_DOXYGEN],[test -n "$DOXYGEN"])

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
