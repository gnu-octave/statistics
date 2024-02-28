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
## @deftypefn  {statistics} {@var{x} =} riceinv (@var{p}, @var{s}, @var{sigma})
##
## Inverse of the Rician distribution (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## of the Rician distribution with non-centrality (distance) parameter @var{s}
## and scale parameter @var{sigma}.  The size of @var{x} is the common size of
## @var{x}, @var{s}, and @var{sigma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Further information about the Rician distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rice_distribution}
##
## @seealso{ricecdf, ricepdf, ricernd, ricefit, ricelike, ricestat}
## @end deftypefn

function x = riceinv (p, s, sigma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("riceinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, S, and B
  if (! isscalar (p) || ! isscalar (s) || ! isscalar (sigma))
    [retval, p, s, sigma] = common_size (p, s, sigma);
    if (retval > 0)
      error ("riceinv: P, S, and B must be of common size or scalars.");
    endif
  endif

  ## Check for P, S, and B being reals
  if (iscomplex (p) || iscomplex (s) || iscomplex (sigma))
    error ("riceinv: P, S, and B must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (s, "single") || isa (sigma, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  k = s < 0 | sigma <= 0 | p < 0 | p > 1 | ...
               isnan (p) | isnan (s) | isnan (sigma);
  x(k) = NaN;

  k = ! k;
  x(k) = sigma(k) .* sqrt (ncx2inv (p(k), 2, (s(k) ./ sigma(k)) .^ 2));

endfunction

%!demo
%! ## Plot various iCDFs from the Rician distribution
%! p = 0.001:0.001:0.999;
%! x1 = riceinv (p, 0, 1);
%! x2 = riceinv (p, 0.5, 1);
%! x3 = riceinv (p, 1, 1);
%! x4 = riceinv (p, 2, 1);
%! x5 = riceinv (p, 4, 1);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-m", p, x5, "-k")
%! grid on
%! legend ({"s = 0, σ = 1", "s = 0.5, σ = 1", "s = 1, σ = 1", ...
%!          "s = 2, σ = 1", "s = 4, σ = 1"}, "location", "northwest")
%! title ("Rician iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p
%! p = [-1 0 0.75 1 2];
%!assert (riceinv (p, ones (1,5), 2*ones (1,5)), [NaN 0 3.5354 Inf NaN], 1e-4)
%!assert (riceinv (p, 1, 2*ones (1,5)), [NaN 0 3.5354 Inf NaN], 1e-4)
%!assert (riceinv (p, ones (1,5), 2), [NaN 0 3.5354 Inf NaN], 1e-4)
%!assert (riceinv (p, [1 0 NaN 1 1], 2), [NaN 0 NaN Inf NaN])
%!assert (riceinv (p, 1, 2*[1 0 NaN 1 1]), [NaN NaN NaN Inf NaN])
%!assert (riceinv ([p(1:2) NaN p(4:5)], 1, 2), [NaN 0 NaN Inf NaN])

## Test class of input preserved
%!assert (riceinv ([p, NaN], 1, 2), [NaN 0 3.5354 Inf NaN NaN], 1e-4)
%!assert (riceinv (single ([p, NaN]), 1, 2), ...
%!        single ([NaN 0 3.5354 Inf NaN NaN]), 1e-4)
%!assert (riceinv ([p, NaN], single (1), 2), ...
%!        single ([NaN 0 3.5354 Inf NaN NaN]), 1e-4)
%!assert (riceinv ([p, NaN], 1, single (2)), ...
%!        single ([NaN 0 3.5354 Inf NaN NaN]), 1e-4)

## Test input validation
%!error<riceinv: function called with too few input arguments.> riceinv ()
%!error<riceinv: function called with too few input arguments.> riceinv (1)
%!error<riceinv: function called with too few input arguments.> riceinv (1,2)
%!error<riceinv: function called with too many inputs> riceinv (1,2,3,4)
%!error<riceinv: P, S, and B must be of common size or scalars.> ...
%! riceinv (ones (3), ones (2), ones (2))
%!error<riceinv: P, S, and B must be of common size or scalars.> ...
%! riceinv (ones (2), ones (3), ones (2))
%!error<riceinv: P, S, and B must be of common size or scalars.> ...
%! riceinv (ones (2), ones (2), ones (3))
%!error<riceinv: P, S, and B must not be complex.> riceinv (i, 2, 2)
%!error<riceinv: P, S, and B must not be complex.> riceinv (2, i, 2)
%!error<riceinv: P, S, and B must not be complex.> riceinv (2, 2, i)
