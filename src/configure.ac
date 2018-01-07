#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([statistics], [1.3.0])

AC_CANONICAL_TARGET

if test "x${build}" != "x${host}"; then
  AC_SUBST([crosscompiling], [true])
else
  AC_SUBST([crosscompiling], [false])
fi  

# Checks for programs.
AC_CHECK_PROG(OCTAVE, octave-cli, octave-cli)
if test -z "$OCTAVE"; then
AC_MSG_ERROR([octave-cli not found], 1);
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
