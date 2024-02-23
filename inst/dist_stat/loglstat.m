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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} loglstat (@var{mu}, @var{sigma})
##
## Compute statistics of the loglogistic distribution.
##
## @code{[@var{m}, @var{v}] = loglstat (@var{mu}, @var{sigma})} returns the mean
## and variance of the loglogistic distribution with mean parameter @var{mu} and
## scale parameter @var{sigma}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the loglogistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-logistic_distribution}
##
## OCTAVE/MATLAB use an alternative parameterization given by the pair
## @math{μ, σ}, i.e. @var{mu} and @var{sigma}, in analogy with the logistic
## distribution.  Their relation to the @math{α} and @math{b} parameters used
## in Wikipedia are given below:
##
## @itemize
## @item @qcode{@var{mu} = log (@var{a})}
## @item @qcode{@var{sigma} = 1 / @var{a}}
## @end itemize
##
## @seealso{logncdf, logninv, lognpdf, lognrnd, lognfit, lognlike}
## @end deftypefn

function [m, v] = loglstat (mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("loglstat: function called with too few input arguments.");
  endif

  ## Check for MU and SIGMA being numeric
  if (! (isnumeric (mu) && isnumeric (sigma)))
    error ("loglstat: MU and SIGMA must be numeric.");
  endif

  ## Check for MU and SIGMA being real
  if (iscomplex (mu) || iscomplex (sigma))
    error ("loglstat: MU and SIGMA must not be complex.");
  endif

  ## Check for common size of MU and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("loglstat: MU and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  pib = pi .* sigma;
  m = (exp (mu) .* pib) ./ sin (pib);
  pib2 = 2 * pib;
  v = exp (mu) .^ 2 .* (pib2 ./ sin (pib2) - pib .^ 2 ./ sin (pib) .^ 2);

  ## Continue argument check
  m(sigma >= 1) = Inf;
  v(sigma >= 0.5) = Inf;
  m(sigma <= 0) = NaN;
  v(sigma <= 0) = NaN;

endfunction

## Input validation tests
%!error<loglstat: function called with too few input arguments.> loglstat ()
%!error<loglstat: function called with too few input arguments.> loglstat (1)
%!error<loglstat: MU and SIGMA must be numeric.> loglstat ({}, 2)
%!error<loglstat: MU and SIGMA must be numeric.> loglstat (1, "")
%!error<loglstat: MU and SIGMA must not be complex.> loglstat (i, 2)
%!error<loglstat: MU and SIGMA must not be complex.> loglstat (1, i)
%!error<loglstat: MU and SIGMA must be of common size or scalars.> ...
%! loglstat (ones (3), ones (2))
%!error<loglstat: MU and SIGMA must be of common size or scalars.> ...
%! loglstat (ones (2), ones (3))

## Output validation tests
%!test
%! [m, v] = loglstat (0, 1);
%! assert (m, Inf, 0.001);
%! assert (v, Inf, 0.001);
%!test
%! [m, v] = loglstat (0, 0.8);
%! assert (m, 4.2758, 0.001);
%! assert (v, Inf, 0.001);
%!test
%! [m, v] = loglstat (0, 0.6);
%! assert (m, 1.9820, 0.001);
%! assert (v, Inf, 0.001);
%!test
%! [m, v] = loglstat (0, 0.4);
%! assert (m, 1.3213, 0.001);
%! assert (v, 2.5300, 0.001);
%!test
%! [m, v] = loglstat (0, 0.2);
%! assert (m, 1.0690, 0.001);
%! assert (v, 0.1786, 0.001);
