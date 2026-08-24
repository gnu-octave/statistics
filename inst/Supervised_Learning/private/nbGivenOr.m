## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.
##

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{v} =} nbGivenOr (@var{given}, @var{dflt})
##
## The argument as it was given, or the default when it was not given at all.
##
## Used to record a naive Bayes model's arguments in @qcode{ModelParameters},
## which reports what was asked for rather than what it was resolved to.
##
## @seealso{ClassificationNaiveBayes}
## @end deftypefn

function v = nbGivenOr (given, dflt)

  if (isempty (given))
    v = dflt;
  else
    v = given;
  endif

endfunction
