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
## @deftypefn  {statistics} {@var{y} =} bisapdf (@var{x}, @var{beta}, @var{gamma})
##
## Birnbaum-Saunders probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Birnbaum-Saunders distribution with scale parameter @var{beta} and
## shape parameter @var{gamma}.  The size of @var{y} is the common size of
## @var{x}, @var{beta}, and @var{gamma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bisacdf, bisapdf, bisarnd, bisafit, bisalike, bisastat}
## @end deftypefn

function y = bisapdf (x, beta, gamma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("bisapdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, BETA and GAMMA
  if (! isscalar (x) || ! isscalar (beta) || ! isscalar (gamma))
    [retval, x, beta, gamma] = common_size (x, beta, gamma);
    if (retval > 0)
      error (strcat (["bisapdf: X, BETA, and GAMMA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, BETA and GAMMA being reals
  if (iscomplex (x) || iscomplex (beta) || iscomplex (gamma))
    error ("bisapdf: X, BETA, and GAMMA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (beta, "single") || isa (gamma, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters.
  k = isnan (x) | ! (beta > 0) | ! (beta < Inf) ...
                | ! (gamma > 0) | ! (gamma < Inf);
  y(k) = NaN;

  ## Find valid values in parameters and data
  k = (x > 0) & (x < Inf) & (beta > 0) & (beta < Inf) ...
                          & (gamma > 0) & (gamma < Inf);
  xk = x(k);

  if (isscalar (beta) && isscalar (gamma))
    z = (sqrt (xk ./ beta) - sqrt (beta ./ xk)) ./ gamma;
    w = (sqrt (xk ./ beta) + sqrt (beta ./ xk)) ./ gamma;
    y(k) = (exp (-0.5 .* z .^ 2) ./ sqrt (2 .* pi)) .* w ./ (2.*xk);
  else
    z = (sqrt (xk ./ beta(k)) - sqrt (beta(k) ./ xk)) ./ gamma(k);
    w = (sqrt (xk ./ beta(k)) + sqrt (beta(k) ./ xk)) ./ gamma(k);
    y(k) = (exp (-0.5 .* z .^ 2) ./ sqrt (2 .* pi)) .* w ./ (2 .* xk);
  endif

endfunction

%!demo
%! ## Plot various PDFs from the Birnbaum-Saunders distribution
%! x = 0.01:0.01:4;
%! y1 = bisapdf (x, 1, 0.5);
%! y2 = bisapdf (x, 1, 1);
%! y3 = bisapdf (x, 1, 2);
%! y4 = bisapdf (x, 1, 5);
%! y5 = bisapdf (x, 1, 10);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
%! grid on
%! ylim ([0, 1.5])
%! legend ({"β = 1 ,γ = 0.5", "β = 1, γ = 1", "β = 1, γ = 2", ...
%!          "β = 1, γ = 5", "β = 1, γ = 10"}, "location", "northeast")
%! title ("Birnbaum-Saunders PDF")
%! xlabel ("values in x")
%! ylabel ("density")

%!demo
%! ## Plot various PDFs from the Birnbaum-Saunders distribution
%! x = 0.01:0.01:6;
%! y1 = bisapdf (x, 1, 0.3);
%! y2 = bisapdf (x, 2, 0.3);
%! y3 = bisapdf (x, 1, 0.5);
%! y4 = bisapdf (x, 3, 0.5);
%! y5 = bisapdf (x, 5, 0.5);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
%! grid on
%! ylim ([0, 1.5])
%! legend ({"β = 1, γ = 0.3", "β = 2, γ = 0.3", "β = 1, γ = 0.5", ...
%!          "β = 3, γ = 0.5", "β = 5, γ = 0.5"}, "location", "northeast")
%! title ("Birnbaum-Saunders CDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 0.3989422804014327, 0.1647717335503959, 0];
%!assert (bisapdf (x, ones (1,5), ones (1,5)), y, eps)
%!assert (bisapdf (x, 1, 1), y, eps)
%!assert (bisapdf (x, 1, ones (1,5)), y, eps)
%!assert (bisapdf (x, ones (1,5), 1), y, eps)
%!assert (bisapdf (x, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisapdf (x, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisapdf ([x, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (bisapdf (single ([x, NaN]), 1, 1), single ([y, NaN]), eps ("single"))
%!assert (bisapdf ([x, NaN], 1, single (1)), single ([y, NaN]), eps ("single"))
%!assert (bisapdf ([x, NaN], single (1), 1), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<bisapdf: function called with too few input arguments.> bisapdf ()
%!error<bisapdf: function called with too few input arguments.> bisapdf (1)
%!error<bisapdf: function called with too few input arguments.> bisapdf (1, 2)
%!error<bisapdf: function called with too many inputs> bisapdf (1, 2, 3, 4)
%!error<bisapdf: X, BETA, and GAMMA must be of common size or scalars.> ...
%! bisapdf (ones (3), ones (2), ones(2))
%!error<bisapdf: X, BETA, and GAMMA must be of common size or scalars.> ...
%! bisapdf (ones (2), ones (3), ones(2))
%!error<bisapdf: X, BETA, and GAMMA must be of common size or scalars.> ...
%! bisapdf (ones (2), ones (2), ones(3))
%!error<bisapdf: X, BETA, and GAMMA must not be complex.> bisapdf (i, 4, 3)
%!error<bisapdf: X, BETA, and GAMMA must not be complex.> bisapdf (1, i, 3)
%!error<bisapdf: X, BETA, and GAMMA must not be complex.> bisapdf (1, 4, i)
