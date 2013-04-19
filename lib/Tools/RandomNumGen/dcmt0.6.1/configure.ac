#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([dcmt], [0.6.1], [BUG-REPORT-ADDRESS])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE([1.8 -Wall -Werror foreign])
AC_CONFIG_SRCDIR([include/dc.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([lib/Makefile])
AC_OUTPUT
