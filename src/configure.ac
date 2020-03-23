#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
### Copyright (C) 2019-2020 John Donoghue <john.donoghue@ieee.org>
### Copyright (C) 2018      Olaf Till <i7tiol@t-online.de>
###
### This program is free software; you can redistribute it and/or
### modify it under the terms of the GNU General Public License as
### published by the Free Software Foundation; either version 3 of the
### License, or (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, see
### <http://www.gnu.org/licenses/>.

AC_PREREQ([2.69])
AC_INIT([statistics], [1.4.2])

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
