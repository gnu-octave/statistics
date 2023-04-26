## Copyright (C) 2018 John Donoghue
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{y} =} bisapdf (@var{x}, @var{a}, @var{b}, @var{mu})
##
## Birnbaum-Saunders probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Birnbaum-Saunders distribution with shape parameter
## @var{a}, scale parameter @var{b}, and location parameter @var{mu}.  The size
## of @var{y} is the common size of @var{x}, @var{a}, @var{b}, and @var{mu}.  A
## scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bbscdf, bbsinv, bbsrnd}
## @end deftypefn

function y = bisapdf (x, a, b, mu)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("bisapdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, A, B, and MU
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b) ...
                     || ! isscalar(mu))
    [retval, x, a, b, mu] = common_size (x, a, b, mu);
    if (retval > 0)
      error (strcat (["bisapdf: X, A, B, and MU must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, A, B, and MU being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b) ...
                    || iscomplex(mu))
    error ("bisapdf: X, A, B, and MU must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single") ...
                        || isa (mu, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters.
  k = isnan (x) | ! (mu > -Inf) | ! (mu < Inf) ...
                | ! (b > 0) | ! (b < Inf) | ! (a > 0) | ! (a < Inf);
  y(k) = NaN;

  ## Find valid values in parameters and data
  k = (x > mu) & (x < Inf) & (-Inf < mu) ...
               & (mu < Inf) & (b > 0) & (b < Inf) & (a > 0) & (a < Inf);

  if (isscalar (a) && isscalar (b) && isscalar (mu))
    x_m = x(k) - mu;
    z = (sqrt (x_m ./ b) - sqrt (b ./ x_m)) ./ a;
    w = (sqrt (x_m ./ b) + sqrt (b ./ x_m)) ./ a;
    y(k) = (exp (-0.5 .* z .^ 2) ./ sqrt (2 .* pi)) .* w ./ (2.*x);
  else
    x_m = x(k) - mu(k);
    z = (sqrt (x_m ./ b(k)) - sqrt (b(k) ./ x_m)) ./ a(k);
    w = (sqrt (x_m ./ b(k)) + sqrt (b(k) ./ x_m)) ./ a(k);
    y(k) = (exp (-0.5 .* z .^ 2) ./ sqrt (2 .* pi)) .* w ./ (2 .* x(k));
  endif

endfunction

%!demo
%! ## Plot various PDFs from the Birnbaum-Saunders distribution
%! x = 0.01:0.01:4;
%! y1 = bisapdf (x, 0.5, 1, 0);
%! y2 = bisapdf (x, 1, 1, 0);
%! y3 = bisapdf (x, 2, 1, 0);
%! y4 = bisapdf (x, 5, 1, 0);
%! y5 = bisapdf (x, 10, 1, 0);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
%! grid on
%! ylim ([0, 1.5])
%! legend ({"α = 0.5, β = 1, μ = 0", "α = 1,    β = 1, μ = 0", ...
%!          "α = 2,    β = 1, μ = 0", "α = 5,    β = 1, μ = 0", ...
%!          "α = 10,  β = 1, μ = 0"}, "location", "northeast")
%! title ("Birnbaum-Saunders CDF")
%! xlabel ("values in x")
%! ylabel ("density")

%!demo
%! ## Plot various PDFs from the Birnbaum-Saunders distribution
%! x = 0.01:0.01:6;
%! y1 = bisapdf (x, 0.3, 1, 0);
%! y2 = bisapdf (x, 0.3, 2, 0);
%! y3 = bisapdf (x, 0.5, 1, 0);
%! y4 = bisapdf (x, 0.5, 3, 0);
%! y5 = bisapdf (x, 0.5, 5, 0);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
%! grid on
%! ylim ([0, 1.5])
%! legend ({"α = 0.3, β = 1, μ = 0", "α = 0.3, β = 2, μ = 0", ...
%!          "α = 0.5, β = 1, μ = 0", "α = 0.5, β = 3, μ = 0", ...
%!          "α = 0.5, β = 5, μ = 0"}, "location", "northeast")
%! title ("Birnbaum-Saunders CDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test results
%!shared x, y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 0.3989422804014327, 0.1647717335503959, 0];
%!assert (bisapdf (x, ones (1,5), ones (1,5), zeros (1,5)), y, eps)
%!assert (bisapdf (x, 1, 1, zeros (1,5)), y, eps)
%!assert (bisapdf (x, 1, ones (1,5), 0), y, eps)
%!assert (bisapdf (x, ones (1,5), 1, 0), y, eps)
%!assert (bisapdf (x, 1, 1, 0), y, eps)
%!assert (bisapdf (x, 1, 1, [0, 0, NaN, 0, 0]), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisapdf (x, 1, [1, 1, NaN, 1, 1], 0), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisapdf (x, [1, 1, NaN, 1, 1], 1, 0), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisapdf ([x, NaN], 1, 1, 0), [y, NaN], eps)

## Test class of input preserved
%!assert (bisapdf (single ([x, NaN]), 1, 1, 0), single ([y, NaN]), eps('single'))
%!assert (bisapdf ([x, NaN], 1, 1, single (0)), single ([y, NaN]), eps('single'))
%!assert (bisapdf ([x, NaN], 1, single (1), 0), single ([y, NaN]), eps('single'))
%!assert (bisapdf ([x, NaN], single (1), 1, 0), single ([y, NaN]), eps('single'))

## Test input validation
%!error bisapdf ()
%!error bisapdf (1)
%!error bisapdf (1,2,3)
%!error bisapdf (1,2,3,4,5)
%!error bisapdf (ones (3), ones (2), ones(2), ones(2))
%!error bisapdf (ones (2), ones (3), ones(2), ones(2))
%!error bisapdf (ones (2), ones (2), ones(3), ones(2))
%!error bisapdf (ones (2), ones (2), ones(2), ones(3))
%!error bisapdf (i, 4, 3, 2)
%!error bisapdf (1, i, 3, 2)
%!error bisapdf (1, 4, i, 2)
%!error bisapdf (1, 4, 3, i)

