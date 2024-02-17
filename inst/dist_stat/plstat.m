## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} plstat (@var{x}, @var{Fx})
##
## Compute statistics of the piecewise linear distribution.
##
## @code{[@var{m}, @var{v}] = plstat (@var{x}, @var{Fx})} returns the mean,
## @var{m}, and variance, @var{v}, of the piecewise linear distribution with a
## vector of @var{x} values at which the CDF changes slope and a vector of CDF
## values @var{Fx} that correspond to each value in @var{x}.  Both @var{x} and
## @var{Fx} must be vectors of the same size and at least 2-elements long.
##
## Further information about the piecewise linear distribution can be found at
## @url{https://en.wikipedia.org/wiki/Piecewise_linear_function}
##
## @seealso{plcdf, plinv, plpdf, plrnd}
## @end deftypefn

function [m, v] = plstat (x, Fx)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("plstat: function called with too few input arguments.");
  endif

  ## Check for common size of X and FX
  if (! isvector (x) || ! isvector (Fx) || ! isequal (size (x), size (Fx)))
    error ("plstat: X and FX must be vectors of equal size.");
  endif

  ## Check for X and FX being at least 2-elements long
  if (length (x) < 2 || length (Fx) < 2)
    error ("plstat: X and FX must be at least two-elements long.");
  endif

  ## Check for Fx being bounded in [0, 1]
  if (any (Fx < 0) || any (Fx > 1))
    error ("plstat: FX must be bounded in the range [0, 1].");
  endif

  ## Check for X and FX being reals
  if (iscomplex (x) || iscomplex (Fx))
    error ("plstat: X and FX must not be complex.");
  endif

  ## Compute the mean and variance
  x_m = (x(1:end-1) + x(2:end)) / 2;
  dFx = diff (Fx);
  m = dot (dFx, x_m);
  x_v = diff(x) .^ 2 / 12;
  v = dot (dFx, x_v + (x_m - m) .^ 2);

endfunction

## Test output
%!shared x, Fx
%! x = [0, 1, 3, 4, 7, 10];
%! Fx = [0, 0.2, 0.5, 0.6, 0.7, 1];
%!assert (plstat (x, Fx), 4.15)
%!test
%! [m, v] = plstat (x, Fx);
%! assert (v, 10.3775, 1e-14)

## Test input validation
%!error<plstat: function called with too few input arguments.> plstat ()
%!error<plstat: function called with too few input arguments.> plstat (1)
%!error<plstat: X and FX must be vectors of equal size.> ...
%! plstat ([0, 1, 2], [0, 1])
%!error<plstat: X and FX must be at least two-elements long.> ...
%! plstat ([0], [1])
%!error<plstat: FX must be bounded in the range> ...
%! plstat ([0, 1, 2], [0, 1, 1.5])
%!error<plstat: FX must be bounded in the range> ...
%! plstat ([0, 1, 2], [0, i, 1])
%!error<plstat: X and FX must not be complex.> ...
%! plstat ([0, i, 2], [0, 0.5, 1])
%!error<plstat: X and FX must not be complex.> ...
%! plstat ([0, i, 2], [0, 0.5i, 1])
