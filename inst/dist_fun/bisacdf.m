## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2018 John Donoghue
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
## @deftypefn  {statistics} {@var{p} =} bisacdf (@var{x}, @var{beta}, @var{gamma})
## @deftypefnx {statistics} {@var{p} =} bisacdf (@var{x}, @var{beta}, @var{gamma}, @qcode{"upper"})
##
## Birnbaum-Saunders cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the Birnbaum-Saunders distribution with scale parameter @var{beta}
## and shape parameter @var{gamma}.  The size of @var{p} is the common size of
## @var{x}, @var{beta} and @var{gamma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## @code{@var{p} = bisacdf (@var{x}, @var{beta}, @var{gamma}, "upper")}
## computes the upper tail probability of the Birnbaum-Saunders distribution
## with parameters @var{beta} and @var{gamma} at the values in @var{x}.
##
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bisainv, bisapdf, bisarnd, bisafit, bisalike, bisastat}
## @end deftypefn

function p = bisacdf (x, beta, gamma, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("bisacdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("bisacdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, BETA, and GAMMA
  if (! isscalar (x) || ! isscalar (beta) || ! isscalar (gamma))
    [retval, x, beta, gamma] = common_size (x, beta, gamma);
    if (retval > 0)
      error (strcat (["bisacdf: X, BETA, and GAMMA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, BETA, and GAMMA being reals
  if (iscomplex (x) || iscomplex (beta) || iscomplex (gamma))
    error ("bisacdf: X, BETA, and GAMMA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (beta, "single") || isa (gamma, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters.
  k = isnan (x) | ! (beta > 0) | ! (beta < Inf) ...
                | ! (gamma > 0) | ! (gamma < Inf);
  p(k) = NaN;

  ## Find valid values in parameters and data
  k = (x > 0) & (x <= Inf) & (beta > 0) & (beta < Inf) ...
                           & (gamma > 0) & (gamma < Inf);
  xk = x(k);

  ## Compute Birnbaum-Saunders CDF
  if (isscalar (beta) && isscalar (gamma))
    if (uflag)
      z = (-sqrt (xk ./ beta) + sqrt (beta ./ xk)) ./ gamma;
    else
      z = (sqrt (xk ./ beta) - sqrt (beta ./ xk)) ./ gamma;
    endif
    p(k) = 0.5 * erfc (-z ./ sqrt (2));
  else
    if (uflag)
      z = (-sqrt (xk ./ beta(k)) + sqrt (beta(k) ./ xk)) ./ gamma(k);
    else
      z = (sqrt (xk ./ beta(k)) - sqrt (beta(k) ./ xk)) ./ gamma(k);
    endif
    p(k) = 0.5 * erfc (-z ./ sqrt (2));
  endif

endfunction

%!demo
%! ## Plot various CDFs from the Birnbaum-Saunders distribution
%! x = 0.01:0.01:4;
%! p1 = bisacdf (x, 1, 0.5);
%! p2 = bisacdf (x, 1, 1);
%! p3 = bisacdf (x, 1, 2);
%! p4 = bisacdf (x, 1, 5);
%! p5 = bisacdf (x, 1, 10);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-m")
%! grid on
%! legend ({"β = 1, γ = 0.5", "β = 1, γ = 1", "β = 1, γ = 2", ...
%!          "β = 1, γ = 5", "β = 1, γ = 10"}, "location", "southeast")
%! title ("Birnbaum-Saunders CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

%!demo
%! ## Plot various CDFs from the Birnbaum-Saunders distribution
%! x = 0.01:0.01:6;
%! p1 = bisacdf (x, 1, 0.3);
%! p2 = bisacdf (x, 2, 0.3);
%! p3 = bisacdf (x, 1, 0.5);
%! p4 = bisacdf (x, 3, 0.5);
%! p5 = bisacdf (x, 5, 0.5);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-m")
%! grid on
%! legend ({"β = 1, γ = 0.3", "β = 2, γ = 0.3", "β = 1, γ = 0.5", ...
%!          "β = 3, γ = 0.5", "β = 5, γ = 0.5"}, "location", "southeast")
%! title ("Birnbaum-Saunders CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test results
%!shared x, y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 1/2, 0.76024993890652337, 1];
%!assert (bisacdf (x, ones (1,5), ones (1,5)), y, eps)
%!assert (bisacdf (x, 1, 1), y, eps)
%!assert (bisacdf (x, 1, ones (1,5)), y, eps)
%!assert (bisacdf (x, ones (1,5), 1), y, eps)
%!assert (bisacdf (x, 1, 1), y, eps)
%!assert (bisacdf (x, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisacdf (x, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisacdf ([x, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (bisacdf (single ([x, NaN]), 1, 1), single ([y, NaN]), eps ("single"))
%!assert (bisacdf ([x, NaN], 1, single (1)), single ([y, NaN]), eps ("single"))
%!assert (bisacdf ([x, NaN], single (1), 1), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<bisacdf: function called with too few input arguments.> bisacdf ()
%!error<bisacdf: function called with too few input arguments.> bisacdf (1)
%!error<bisacdf: function called with too few input arguments.> bisacdf (1, 2)
%!error<bisacdf: function called with too many inputs> ...
%! bisacdf (1, 2, 3, 4, 5)
%!error<bisacdf: invalid argument for upper tail.> bisacdf (1, 2, 3, "tail")
%!error<bisacdf: invalid argument for upper tail.> bisacdf (1, 2, 3, 4)
%!error<bisacdf: X, BETA, and GAMMA must be of common size or scalars.> ...
%! bisacdf (ones (3), ones (2), ones(2))
%!error<bisacdf: X, BETA, and GAMMA must be of common size or scalars.> ...
%! bisacdf (ones (2), ones (3), ones(2))
%!error<bisacdf: X, BETA, and GAMMA must be of common size or scalars.> ...
%! bisacdf (ones (2), ones (2), ones(3))
%!error<bisacdf: X, BETA, and GAMMA must not be complex.> bisacdf (i, 4, 3)
%!error<bisacdf: X, BETA, and GAMMA must not be complex.> bisacdf (1, i, 3)
%!error<bisacdf: X, BETA, and GAMMA must not be complex.> bisacdf (1, 4, i)

