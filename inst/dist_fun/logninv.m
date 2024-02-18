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
## @deftypefn  {statistics} {@var{x} =} logninv (@var{p})
## @deftypefnx {statistics} {@var{x} =} logninv (@var{p}, @var{mu})
## @deftypefnx {statistics} {@var{x} =} logninv (@var{p}, @var{mu}, @var{sigma})
##
## Inverse of the log-normal cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the log-normal distribution with mean parameter @var{mu} and standard
## deviation parameter @var{sigma}, each corresponding to the associated normal
## distribution.  The size of @var{x} is the common size of @var{p}, @var{mu},
## and @var{sigma}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## If a random variable follows this distribution, its logarithm is normally
## distributed with mean @var{mu} and standard deviation @var{sigma}.
##
## Default parameter values are @qcode{@var{mu} = 0} and
## @qcode{@var{sigma} = 1}.  Both parameters must be reals and
## @qcode{@var{sigma} > 0}.  For @qcode{@var{sigma} <= 0}, @qcode{NaN} is
## returned.
##
## Further information about the log-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-normal_distribution}
##
## @seealso{logncdf, lognpdf, lognrnd, lognfit, lognlike, lognstat}
## @end deftypefn

function x = logninv (p, mu = 0, sigma = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of P, MU, and SIGMA
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar (sigma))
    [retval, p, mu, sigma] = common_size (p, mu, sigma);
    if (retval > 0)
      error ("logninv: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (sigma))
    error ("logninv: X, MU, and SIGMA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") || isa (sigma, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Compute lognormal iCDF
  k = !(p >= 0) | !(p <= 1) | !(sigma > 0) | !(sigma < Inf);
  x(k) = NaN;

  k = (p == 1) & (sigma > 0) & (sigma < Inf);
  x(k) = Inf;

  k = (p >= 0) & (p < 1) & (sigma > 0) & (sigma < Inf);
  if (isscalar (mu) && isscalar (sigma))
    x(k) = exp (mu) .* exp (sigma .* (-sqrt (2) * erfcinv (2 * p(k))));
  else
    x(k) = exp (mu(k)) .* exp (sigma(k) .* (-sqrt (2) * erfcinv (2 * p(k))));
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the log-normal distribution
%! p = 0.001:0.001:0.999;
%! x1 = logninv (p, 0, 1);
%! x2 = logninv (p, 0, 0.5);
%! x3 = logninv (p, 0, 0.25);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r")
%! grid on
%! ylim ([0, 3])
%! legend ({"μ = 0, σ = 1", "μ = 0, σ = 0.5", "μ = 0, σ = 0.25"}, ...
%!         "location", "northwest")
%! title ("Log-normal iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (logninv (p, ones (1,5), ones (1,5)), [NaN 0 e Inf NaN])
%!assert (logninv (p, 1, ones (1,5)), [NaN 0 e Inf NaN])
%!assert (logninv (p, ones (1,5), 1), [NaN 0 e Inf NaN])
%!assert (logninv (p, [1 1 NaN 0 1], 1), [NaN 0 NaN Inf NaN])
%!assert (logninv (p, 1, [1 0 NaN Inf 1]), [NaN NaN NaN NaN NaN])
%!assert (logninv ([p(1:2) NaN p(4:5)], 1, 2), [NaN 0 NaN Inf NaN])

## Test class of input preserved
%!assert (logninv ([p, NaN], 1, 1), [NaN 0 e Inf NaN NaN])
%!assert (logninv (single ([p, NaN]), 1, 1), single ([NaN 0 e Inf NaN NaN]))
%!assert (logninv ([p, NaN], single (1), 1), single ([NaN 0 e Inf NaN NaN]))
%!assert (logninv ([p, NaN], 1, single (1)), single ([NaN 0 e Inf NaN NaN]))

## Test input validation
%!error logninv ()
%!error logninv (1,2,3,4)
%!error logninv (ones (3), ones (2), ones (2))
%!error logninv (ones (2), ones (3), ones (2))
%!error logninv (ones (2), ones (2), ones (3))
%!error logninv (i, 2, 2)
%!error logninv (2, i, 2)
%!error logninv (2, 2, i)
