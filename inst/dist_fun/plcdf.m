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
## @deftypefn  {statistics} {@var{p} =} plcdf (@var{data}, @var{x}, @var{Fx})
## @deftypefnx {statistics} {@var{p} =} plcdf (@var{data}, @var{x}, @var{Fx}, @qcode{"upper"})
##
## Piecewise linear cumulative distribution function (CDF).
##
## For each element of @var{data}, compute the cumulative distribution function
## (CDF) of the piecewise linear distribution with a vector of @var{x} values at
## which the CDF changes slope and a vector of CDF values @var{Fx} that
## correspond to each value in @var{x}.  Both @var{x} and @var{Fx} must be
## vectors of the same size and at least 2-elements long.  The size of @var{p}
## is the same as @var{data}.
##
## @code{@var{p} = plcdf (@var{data}, @var{x}, @var{Fx}, "upper")} computes
## the upper tail probability of the piecewise linear distribution with
## parameters @var{x} and @var{Fx}, at the values in @var{data}.
##
## Further information about the piecewise linear distribution can be found at
## @url{https://en.wikipedia.org/wiki/Piecewise_linear_function}
##
## @seealso{plinv, plpdf, plrnd, plstat}
## @end deftypefn

function p = plcdf (data, x, Fx, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("plcdf: function called with too few input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin == 4 && strcmpi (uflag, "upper"))
    uflag = true;
  elseif (nargin == 4  && ! strcmpi (uflag, "upper"))
    error ("plcdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Check for common size of X and FX
  if (! isvector (x) || ! isvector (Fx) || ! isequal (size (x), size (Fx)))
    error ("plcdf: X and FX must be vectors of equal size.");
  endif

  ## Check for X and FX being at least 2-elements long
  if (length (x) < 2 || length (Fx) < 2)
    error ("plcdf: X and FX must be at least two-elements long.");
  endif

  ## Check for Fx being bounded in [0, 1]
  if (any (Fx < 0) || any (Fx > 1))
    error ("plcdf: FX must be bounded in the range [0, 1].");
  endif

  ## Check for DATA, X, and FX being reals
  if (iscomplex (data) || iscomplex (x) || iscomplex (Fx))
    error ("plcdf: DATA, X, and FX must not be complex.");
  endif

  ## Check for class type
  if (isa (data, "single") || isa (x, "single") || isa (Fx, "single"));
    p = zeros (size (data), "single");
  else
    p = zeros (size (data));
  endif

  ## Find data within supported range
  support = (data >= x(1) & data <= x(end));
  p(support) = interp1 (x, Fx, data(support), "linear");

  ## Force right side outside support to 1 and invalid data to NaN
  p(data > x(end)) = 1;
  p(isnan(data)) = NaN;

  ## Return upper tail (if requested)
  if (uflag)
    p = 1 - p;
  endif

endfunction

%!demo
%! ## Plot various CDFs from the Piecewise linear distribution
%! data = 0:0.01:10;
%! x1 = [0, 1, 3, 4, 7, 10];
%! Fx1 = [0, 0.2, 0.5, 0.6, 0.7, 1];
%! x2 = [0, 2, 5, 6, 7, 8];
%! Fx2 = [0, 0.1, 0.3, 0.6, 0.9, 1];
%! p1 = plcdf (data, x1, Fx1);
%! p2 = plcdf (data, x2, Fx2);
%! plot (data, p1, "-b", data, p2, "g")
%! grid on
%! ylim ([0, 1])
%! xlim ([0, 10])
%! legend ({"x1, Fx1", "x2, Fx2"}, "location", "southeast")
%! title ("Piecewise linear CDF")
%! xlabel ("values in data")
%! ylabel ("probability")

## Test output
%!test
%! data = 0:0.2:1;
%! p = plcdf (data, [0, 1], [0, 1]);
%! assert (p, data);
%!test
%! data = 0:0.2:1;
%! p = plcdf (data, [0, 2], [0, 1]);
%! assert (p, 0.5 * data);
%!test
%! data = 0:0.2:1;
%! p = plcdf (data, [0, 1], [0, 0.5]);
%! assert (p, 0.5 * data);
%!test
%! data = 0:0.2:1;
%! p = plcdf (data, [0, 0.5], [0, 1]);
%! assert (p, [0, 0.4, 0.8, 1, 1, 1]);
%!test
%! data = 0:0.2:1;
%! p = plcdf (data, [0, 1], [0, 1], "upper");
%! assert (p, 1 - data);

## Test input validation
%!error<plcdf: function called with too few input arguments.> plcdf ()
%!error<plcdf: function called with too few input arguments.> plcdf (1)
%!error<plcdf: function called with too few input arguments.> plcdf (1, 2)
%!error<plcdf: invalid argument for upper tail.> plcdf (1, 2, 3, "uper")
%!error<plcdf: invalid argument for upper tail.> plcdf (1, 2, 3, 4)
%!error<plcdf: X and FX must be vectors of equal size.> ...
%! plcdf (1, [0, 1, 2], [0, 1])
%!error<plcdf: X and FX must be at least two-elements long.> ...
%! plcdf (1, [0], [1])
%!error<plcdf: FX must be bounded in the range> ...
%! plcdf (1, [0, 1, 2], [0, 1, 1.5])
%!error<plcdf: FX must be bounded in the range> ...
%! plcdf (1, [0, 1, 2], [0, i, 1])
%!error<plcdf: DATA, X, and FX must not be complex.> ...
%! plcdf (i, [0, 1, 2], [0, 0.5, 1])
%!error<plcdf: DATA, X, and FX must not be complex.> ...
%! plcdf (1, [0, i, 2], [0, 0.5, 1])
%!error<plcdf: DATA, X, and FX must not be complex.> ...
%! plcdf (1, [0, 1, 2], [0, 0.5i, 1])
