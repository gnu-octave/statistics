#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([statistics], [1.4.0])

AC_CANONICAL_TARGET

if test "x${build}" != "x${host}"; then
  AC_SUBST([crosscompiling], [true])
else
  AC_SUBST([crosscompiling], [false])
fi  

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
