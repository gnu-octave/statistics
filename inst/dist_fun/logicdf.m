## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} logicdf (@var{x}, @var{mu}, @var{sigma})
## @deftypefnx {statistics} {@var{p} =} logicdf (@var{x}, @var{mu}, @var{sigma}, @qcode{"upper"})
##
## Logistic cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the logistic distribution with location parameter @var{mu} and scale
## parameter @var{sigma}.  The size of @var{p} is the common size of @var{x},
## @var{mu}, and @var{sigma}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Both parameters must be reals and @qcode{@var{sigma} > 0}.
## For @qcode{@var{sigma} <= 0}, @qcode{NaN} is returned.
##
## @code{@var{p} = logicdf (@var{x}, @var{mu}, @var{sigma}, "upper")} computes
## the upper tail probability of the logistic distribution with parameters
## @var{mu} and @var{sigma}, at the values in @var{x}.
##
## Further information about the logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Logistic_distribution}
##
## @seealso{logiinv, logipdf, logirnd, logifit, logilike, logistat}
## @end deftypefn

function p = logicdf (x, mu, sigma, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("logicdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("logicdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, MU, and SIGMA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(sigma))
    [retval, x, mu, sigma] = common_size (x, mu, sigma);
    if (retval > 0)
      error ("logicdf: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("logicdf: X, MU, and SIGMA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"));
    p = NaN (size (x), "single");
  else
    p = NaN (size (x));
  endif

  ## Find normal and edge cases
  k1 = (x == -Inf) & (sigma > 0);
  k2 = (x == Inf) & (sigma > 0);
  k = ! k1 & ! k2 & (sigma > 0);

  ## Compute logistic CDF
  if (uflag)
    p(k1) = 1;
    p(k2) = 0;
    p(k) = 1 ./ (1 + exp ((x(k) - mu(k)) ./ sigma(k)));
  else
    p(k1) = 0;
    p(k2) = 1;
    p(k) = 1 ./ (1 + exp (- (x(k) - mu(k)) ./ sigma(k)));
  endif

endfunction

%!demo
%! ## Plot various CDFs from the logistic distribution
%! x = -5:0.01:20;
%! p1 = logicdf (x, 5, 2);
%! p2 = logicdf (x, 9, 3);
%! p3 = logicdf (x, 9, 4);
%! p4 = logicdf (x, 6, 2);
%! p5 = logicdf (x, 2, 1);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-m")
%! grid on
%! legend ({"μ = 5, σ = 2", "μ = 9, σ = 3", "μ = 9, σ = 4", ...
%!          "μ = 6, σ = 2", "μ = 2, σ = 1"}, "location", "southeast")
%! title ("Logistic CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-Inf -log(3) 0 log(3) Inf];
%! y = [0, 1/4, 1/2, 3/4, 1];
%!assert (logicdf ([x, NaN], 0, 1), [y, NaN], eps)
%!assert (logicdf (x, 0, [-2, -1, 0, 1, 2]), [nan(1, 3), 0.75, 1], eps)

## Test class of input preserved
%!assert (logicdf (single ([x, NaN]), 0, 1), single ([y, NaN]), eps ("single"))
%!assert (logicdf ([x, NaN], single (0), 1), single ([y, NaN]), eps ("single"))
%!assert (logicdf ([x, NaN], 0, single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<logicdf: function called with too few input arguments.> logicdf ()
%!error<logicdf: function called with too few input arguments.> logicdf (1)
%!error<logicdf: function called with too few input arguments.> ...
%! logicdf (1, 2)
%!error<logicdf: invalid argument for upper tail.> logicdf (1, 2, 3, "tail")
%!error<logicdf: invalid argument for upper tail.> logicdf (1, 2, 3, 4)
%!error<logicdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! logicdf (1, ones (2), ones (3))
%!error<logicdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! logicdf (ones (2), 1, ones (3))
%!error<logicdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! logicdf (ones (2), ones (3), 1)
%!error<logicdf: X, MU, and SIGMA must not be complex.> logicdf (i, 2, 3)
%!error<logicdf: X, MU, and SIGMA must not be complex.> logicdf (1, i, 3)
%!error<logicdf: X, MU, and SIGMA must not be complex.> logicdf (1, 2, i)
