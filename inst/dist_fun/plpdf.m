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
## @deftypefn  {statistics} {@var{y} =} plpdf (@var{data}, @var{x}, @var{Fx})
##
## Piecewise linear probability density function (PDF).
##
## For each element of @var{data}, compute the probability density function
## (PDF) of the piecewise linear distribution with a vector of @var{x} values at
## which the CDF changes slope and a vector of CDF values @var{Fx} that
## correspond to each value in @var{x}.  Both @var{x} and @var{Fx} must be
## vectors of the same size and at least 2-elements long.  The size of @var{p}
## is the same as @var{data}.
##
## Further information about the piecewise linear distribution can be found at
## @url{https://en.wikipedia.org/wiki/Piecewise_linear_function}
##
## @seealso{plcdf, plinv, plrnd, plstat}
## @end deftypefn

function y = plpdf (data, x, Fx)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("plpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X and FX
  if (! isvector (x) || ! isvector (Fx) || ! isequal (size (x), size (Fx)))
    error ("plpdf: X and FX must be vectors of equal size.");
  endif

  ## Check for X and FX being at least 2-elements long
  if (length (x) < 2 || length (Fx) < 2)
    error ("plpdf: X and FX must be at least two-elements long.");
  endif

  ## Check for Fx being bounded in [0, 1]
  if (any (Fx < 0) || any (Fx > 1))
    error ("plpdf: FX must be bounded in the range [0, 1].");
  endif

  ## Check for DATA, X, and FX being reals
  if (iscomplex (data) || iscomplex (x) || iscomplex (Fx))
    error ("plpdf: DATA, X, and FX must not be complex.");
  endif

  ## Check for class type
  if (isa (data, "single") || isa (x, "single") || isa (Fx, "single"));
    y = zeros (size (data), "single");
  else
    y = zeros (size (data));
  endif

  ## Bin data according to X
  [~, bin] = histc (data, [-Inf, x, Inf]);

  ## Compute piecewise densities
  dense = diff(Fx) ./ diff(x);
  bin_d = [0, dense, 0];

  ## Fix densities
  xlen = length (x);
  bin(bin > xlen) = xlen + 1;
  y(bin>0) = bin_d(bin(bin>0));

  ## Force invalid data to NaN
  y(isnan(data)) = NaN;

endfunction

%!demo
%! ## Plot various PDFs from the Piecewise linear distribution
%! data = 0:0.01:10;
%! x1 = [0, 1, 3, 4, 7, 10];
%! Fx1 = [0, 0.2, 0.5, 0.6, 0.7, 1];
%! x2 = [0, 2, 5, 6, 7, 8];
%! Fx2 = [0, 0.1, 0.3, 0.6, 0.9, 1];
%! y1 = plpdf (data, x1, Fx1);
%! y2 = plpdf (data, x2, Fx2);
%! plot (data, y1, "-b", data, y2, "g")
%! grid on
%! ylim ([0, 0.6])
%! xlim ([0, 10])
%! legend ({"x1, Fx1", "x2, Fx2"}, "location", "northeast")
%! title ("Piecewise linear CDF")
%! xlabel ("values in data")
%! ylabel ("density")

## Test output
%!shared x, Fx
%! x = [0, 1, 3, 4, 7, 10];
%! Fx = [0, 0.2, 0.5, 0.6, 0.7, 1];
%!assert (plpdf (0.5, x, Fx), 0.2, eps);
%!assert (plpdf (1.5, x, Fx), 0.15, eps);
%!assert (plpdf (3.5, x, Fx), 0.1, eps);
%!assert (plpdf (5, x, Fx), 0.1/3, eps);
%!assert (plpdf (8, x, Fx), 0.1, eps);

## Test input validation
%!error<plpdf: function called with too few input arguments.> plpdf ()
%!error<plpdf: function called with too few input arguments.> plpdf (1)
%!error<plpdf: function called with too few input arguments.> plpdf (1, 2)
%!error<plpdf: X and FX must be vectors of equal size.> ...
%! plpdf (1, [0, 1, 2], [0, 1])
%!error<plpdf: X and FX must be at least two-elements long.> ...
%! plpdf (1, [0], [1])
%!error<plpdf: FX must be bounded in the range> ...
%! plpdf (1, [0, 1, 2], [0, 1, 1.5])
%!error<plpdf: FX must be bounded in the range> ...
%! plpdf (1, [0, 1, 2], [0, i, 1])
%!error<plpdf: DATA, X, and FX must not be complex.> ...
%! plpdf (i, [0, 1, 2], [0, 0.5, 1])
%!error<plpdf: DATA, X, and FX must not be complex.> ...
%! plpdf (1, [0, i, 2], [0, 0.5, 1])
%!error<plpdf: DATA, X, and FX must not be complex.> ...
%! plpdf (1, [0, 1, 2], [0, 0.5i, 1])
