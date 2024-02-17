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
## @deftypefn  {statistics} {@var{data} =} plinv (@var{p}, @var{x}, @var{Fx})
##
## Inverse of the piecewise linear distribution (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## of the piecewise linear distribution with a vector of @var{x} values at
## which the CDF changes slope and a vector of CDF values @var{Fx} that
## correspond to each value in @var{x}.  Both @var{x} and @var{Fx} must be
## vectors of the same_p size and at least 2-elements long..  The size of
## @var{data} is the same_p as @var{p}.
##
## Further information about the piecewise linear distribution can be found at
## @url{https://en.wikipedia.org/wiki/Piecewise_linear_function}
##
## @seealso{plcdf, plpdf, plrnd, plstat}
## @end deftypefn

function data = plinv (p, x, Fx)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("plinv: function called with too few input arguments.");
  endif

  ## Check for common size of X and FX
  if (! isvector (x) || ! isvector (Fx) || ! isequal (size (x), size (Fx)))
    error ("plinv: X and FX must be vectors of equal size.");
  endif

  ## Check for X and FX being at least 2-elements long
  if (length (x) < 2 || length (Fx) < 2)
    error ("plinv: X and FX must be at least two-elements long.");
  endif

  ## Check for Fx being bounded in [0, 1]
  if (any (Fx < 0) || any (Fx > 1))
    error ("plinv: FX must be bounded in the range [0, 1].");
  endif

  ## Check for P, X, and FX being reals
  if (iscomplex (p) || iscomplex (x) || iscomplex (Fx))
    error ("plinv: P, X, and FX must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (x, "single") || isa (Fx, "single"));
    data = zeros (size (p), "single");
  else
    data = zeros (size (p));
  endif

  ## Remove consecutive bins with almost zero probability
  pw_diff = diff (Fx);
  if any(pw_diff==0)
    zero_p = 2 * eps(Fx);
    same_p = pw_diff <= zero_p(1:end-1);
    remove = same_p(1:end-1) & same_p(2:end);
    while (any(remove))
      idx = find (remove);
      same_p(idx) = [];
      Fx(idx+1) = [];
      x(idx+1) = [];
      pw_diff = diff (Fx);
      remove = same_p(1:end-1) & same_p(2:end);
    endwhile
    idx = find (pw_diff==0);
    Fx(idx+1) = Fx(idx) + eps(Fx(idx));
  endif
  p(p < 0 | 1 < p) = NaN;
  data = interp1 (Fx, x, p, "linear");

endfunction

%!demo
%! ## Plot various iCDFs from the Piecewise linear distribution
%! p = 0.001:0.001:0.999;
%! x1 = [0, 1, 3, 4, 7, 10];
%! Fx1 = [0, 0.2, 0.5, 0.6, 0.7, 1];
%! x2 = [0, 2, 5, 6, 7, 8];
%! Fx2 = [0, 0.1, 0.3, 0.6, 0.9, 1];
%! data1 = plinv (p, x1, Fx1);
%! data2 = plinv (p, x2, Fx2);
%! plot (p, data1, "-b", p, data2, "-g")
%! grid on
%! legend ({"x1, Fx1", "x2, Fx2"}, "location", "northwest")
%! title ("Piecewise linear iCDF")
%! xlabel ("probability")
%! ylabel ("values in data")

## Test output
%!test
%! p = 0:0.2:1;
%! data = plinv (p, [0, 1], [0, 1]);
%! assert (data, p);
%!test
%! p = 0:0.2:1;
%! data = plinv (p, [0, 2], [0, 1]);
%! assert (data, 2 * p);
%!test
%! p = 0:0.2:1;
%! data_out = 1:6;
%! data = plinv (p, [0, 1], [0, 0.5]);
%! assert (data, [0, 0.4, 0.8, NA, NA, NA]);
%!test
%! p = 0:0.2:1;
%! data_out = 1:6;
%! data = plinv (p, [0, 0.5], [0, 1]);
%! assert (data, [0:0.1:0.5]);

## Test input validation
%!error<plinv: function called with too few input arguments.> plinv ()
%!error<plinv: function called with too few input arguments.> plinv (1)
%!error<plinv: function called with too few input arguments.> plinv (1, 2)
%!error<plinv: X and FX must be vectors of equal size.> ...
%! plinv (1, [0, 1, 2], [0, 1])
%!error<plinv: X and FX must be at least two-elements long.> ...
%! plinv (1, [0], [1])
%!error<plinv: FX must be bounded in the range> ...
%! plinv (1, [0, 1, 2], [0, 1, 1.5])
%!error<plinv: FX must be bounded in the range> ...
%! plinv (1, [0, 1, 2], [0, i, 1])
%!error<plinv: P, X, and FX must not be complex.> ...
%! plinv (i, [0, 1, 2], [0, 0.5, 1])
%!error<plinv: P, X, and FX must not be complex.> ...
%! plinv (1, [0, i, 2], [0, 0.5, 1])
%!error<plinv: P, X, and FX must not be complex.> ...
%! plinv (1, [0, 1, 2], [0, 0.5i, 1])
