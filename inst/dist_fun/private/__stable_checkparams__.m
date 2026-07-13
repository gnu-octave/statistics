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

## -*- texinfo -*-
## @deftypefn {Function File} {@var{msg} =} __stable_checkparams__ (@var{alpha}, @var{beta}, @var{gam}, @var{delta})
##
## Validate the parameters of a stable distribution.
##
## Returns an empty string if the tail index @var{alpha}, skewness @var{beta},
## scale @var{gam}, and location @var{delta} are valid, otherwise a message body
## describing the first violation.  The calling function prepends its own name.
## This is a private helper shared by the @code{stbl*} functions.
##
## @end deftypefn

function msg = __stable_checkparams__ (alpha, beta, gam, delta)

  msg = "";
  if (! (isscalar (alpha) && isreal (alpha) && alpha > 0 && alpha <= 2))
    msg = "ALPHA must be a scalar in the range (0, 2].";
  elseif (! (isscalar (beta) && isreal (beta) && beta >= -1 && beta <= 1))
    msg = "BETA must be a scalar in the range [-1, 1].";
  elseif (! (isscalar (gam) && isreal (gam) && gam > 0))
    msg = "GAM must be a positive scalar.";
  elseif (! (isscalar (delta) && isreal (delta)))
    msg = "DELTA must be a real scalar.";
  endif

endfunction
