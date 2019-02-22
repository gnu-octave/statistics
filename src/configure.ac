#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([statistics], [1.4.1])

AC_CANONICAL_TARGET

if test "x${build}" != "x${host}"; then
  AC_SUBST([crosscompiling], [true])
else
  AC_SUBST([crosscompiling], [false])
fi

AC_PROG_CC


# Check for Windows
AC_MSG_CHECKING([for windows])
have_windows=no
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
   #include <windows.h>
   #ifndef __WIN32__
     #error "Not windows!"
   #endif
  ]], [])],
 [AC_MSG_RESULT([yes])
  have_windows=yes],
 [AC_MSG_RESULT([no])])

# if windows, we can build some of the packages just by having
# windows.h
if test $have_windows = yes; then
  AC_SUBST([HAVE_WINDOWS], [1], [build for windows target])

  AC_CHECK_PROG(OCTAVE_CLI, octave-cli, octave-cli)
else
  AC_SUBST([HAVE_WINDOWS], [0], [build for not a windows target])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
