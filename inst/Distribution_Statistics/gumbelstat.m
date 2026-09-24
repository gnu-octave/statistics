## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} gumbelstat (@var{mu}, @var{beta})
##
## Compute statistics of the Gumbel distribution.
##
## @code{[@var{m}, @var{v}] = gumbelstat (@var{mu}, @var{beta})} returns the
## mean and variance of the Gumbel distribution (also known as the extreme value
## or the type I generalized extreme value distribution) with location
## parameter @var{mu} and scale parameter @var{beta}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## The mean is @code{@var{mu} + @var{gamma} * @var{beta}}, where @var{gamma}
## is the Euler-Mascheroni constant, and the variance is
## @code{pi^2 * @var{beta}^2 / 6}.  A non-positive @var{beta} returns
## @code{NaN} for both.
##
## The Gumbel distribution is used to model the distribution of the maximum (or
## the minimum) of a number of samples of various distributions.  This version
## is suitable for modeling maxima.  For modeling minima, use the alternative
## extreme value statistics, @code{evstat}.
##
## Further information about the Gumbel distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gumbel_distribution}
##
## Input arguments must be @qcode{double} or @qcode{single}; integer, logical,
## and character arrays are rejected.
##
## @seealso{gumbelcdf, gumbelinv, gumbelpdf, gumbelrnd, gumbelfit, gumbellike,
## evstat}
## @end deftypefn

function [m, v] = gumbelstat (mu, beta)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("gumbelstat: too few input arguments.");
  endif

  ## Check for common size of MU and BETA
  if (! isscalar (mu) || ! isscalar (beta))
    [err, mu, beta] = common_size (mu, beta);
    if (err > 0)
      error ("gumbelstat: MU and BETA must be of common size or scalars.");
    endif
  endif

  ## Check for MU and BETA being double or single
  if (! (isfloat (mu) && isfloat (beta)))
    error ("gumbelstat: MU and BETA must be double or single.");
  endif

  ## Check for MU and BETA being reals
  if (iscomplex (mu) || iscomplex (beta))
    error ("gumbelstat: MU and BETA must not be complex.");
  endif

  ## Return NaNs for out of range parameters
  beta(beta <= 0) = NaN;

  ## Calculate mean and variance; psi (1) is minus the Euler-Mascheroni
  ## constant
  m = mu - psi (1) .* beta;
  v = (pi .* beta) .^ 2 ./ 6;

endfunction

## Test output
%!test
%! [m, v] = gumbelstat (0, 1);
%! assert_equal (m, 0.577215664901533, 1e-15);
%! assert_equal (v, pi ^ 2 / 6, eps);
%!test
%! [m, v] = gumbelstat ([-5, 0, 1, 2, 3], [0, 1, 2, -1, 3]);
%! assert_equal (m, [NaN, 0.577215664901533, 2.154431329803066, NaN, ...
%!                   4.731646994704599], 1e-14);
%! assert_equal (v, [NaN, pi^2/6, 2*pi^2/3, NaN, 3*pi^2/2], 1e-14);
%!test
%! [m, v] = gumbelstat (2, 3);
%! f = @(x) gumbelpdf (x, 2, 3);
%! assert_equal (m, integral (@(x) x .* f(x), -Inf, Inf), 1e-8);
%! assert_equal (v, integral (@(x) (x - m) .^ 2 .* f(x), -Inf, Inf), 1e-8);
%!test
%! [m, v] = gumbelstat (1, 2);
%! [me, ve] = evstat (-1, 2);
%! assert_equal (m, -me, 1e-14);
%! assert_equal (v, ve, 1e-14);
%!test
%! [m, v] = gumbelstat (0, single (1));
%! assert_equal (class (m), 'single');
%! assert_equal (class (v), 'single');

## Test input validation
%!error<gumbelstat: too few input arguments.> gumbelstat ()
%!error<gumbelstat: too few input arguments.> gumbelstat (1)
%!error<gumbelstat: MU and BETA must be of common size or scalars.> ...
%! gumbelstat (ones (3), ones (2))
%!error<gumbelstat: MU and BETA must be double or single.> ...
%! gumbelstat (int32 (1), 2)
%!error<gumbelstat: MU and BETA must be double or single.> ...
%! gumbelstat (1, true)
%!error<gumbelstat: MU and BETA must be double or single.> ...
%! gumbelstat ('a', 2)
%!error<gumbelstat: MU and BETA must not be complex.> gumbelstat (i, 2)
%!error<gumbelstat: MU and BETA must not be complex.> gumbelstat (1, i)
