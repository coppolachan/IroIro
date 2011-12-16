#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([CommonCode++], [0.6.6], [cossu@post.kek.jp])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE([1.8 -Wall -Werror foreign])
AC_CONFIG_HEADERS([include/config.h])
AC_PROG_RANLIB


# Checks for programs.
AC_PROG_CXX

case $host in 
rs6000-*-aix*) AC_PROG_CXX(mpCC);;
esac

AC_PROG_CC
AC_LANG(C++)

# Checks for libraries.
AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([gslcblas],
	[cblas_dgemm],,
	[AC_MSG_ERROR(GNU Scientific Library BLAS not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.gnu.org/s/gsl/)])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,
	[AC_MSG_ERROR(GNU Scientific Library not found in your system.
Please install or provide the correct path to your installation
Info at: http://www.gnu.org/s/gsl/)])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt])

# MPI communicator
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

if [ $mpi || $mpihitachi ]
then

AC_CHECK_HEADERS([/usr/include/mpich2/mpi.h])
AC_DEFINE([HAVE_MPI],[],[Processing MPI routines])
fi

if [ $mpihitachi ]
then
AC_DEFINE([BIG_ENDIAN],[],[Using Big Endian machine type])
AC_DEFINE([HITACHISR16K],[],[Specific tweaks for Hitachi SR16K])
fi
# Checks for doxygen support
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

AM_CONDITIONAL([HAVE_DOXYGEN],
[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([docs/doxy.cfg])])

if [ $mpihitachi ] 
then
ARFLAGS="-X64 -c -ru"
else
ARFLAGS="cru"
fi

AC_SUBST(ARFLAGS)

AC_CONFIG_FILES([Makefile
 		lib/Makefile
		tests/Makefile
		docs/Makefile])

		 
AC_OUTPUT
