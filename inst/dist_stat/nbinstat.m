## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} nbinstat (@var{r}, @var{ps})
##
## Compute statistics of the negative binomial distribution.
##
## @code{[@var{m}, @var{v}] = nbinstat (@var{r}, @var{ps})} returns the mean
## and variance of the negative binomial distribution with parameters @var{r}
## and @var{ps}, where @var{r} is the number of successes until the experiment
## is stopped and @var{ps} is the probability of success in each experiment,
## given the number of failures in @var{x}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the negative binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Negative_binomial_distribution}
##
## @seealso{nbincdf, nbininv, nbininv, nbinrnd, nbinfit, nbinlike}
## @end deftypefn

function [m, v] = nbinstat (r, ps)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("nbinstat: function called with too few input arguments.");
  endif

  ## Check for R and PS being numeric
  if (! (isnumeric (r) && isnumeric (ps)))
    error ("nbinstat: R and PS must be numeric.");
  endif

  ## Check for R and PS being real
  if (iscomplex (r) || iscomplex (ps))
    error ("nbinstat: R and PS must not be complex.");
  endif

  ## Check for common size of R and PS
  if (! isscalar (r) || ! isscalar (ps))
    [retval, r, ps] = common_size (r, ps);
    if (retval > 0)
      error ("nbinstat: R and PS must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  q = 1 - ps;
  m = r .* q ./ ps;
  v = r .* q ./ (ps .^ 2);

  ## Continue argument check
  k = find (! (r > 0) | ! (r < Inf) | ! (ps > 0) | ! (ps < 1));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<nbinstat: function called with too few input arguments.> nbinstat ()
%!error<nbinstat: function called with too few input arguments.> nbinstat (1)
%!error<nbinstat: R and PS must be numeric.> nbinstat ({}, 2)
%!error<nbinstat: R and PS must be numeric.> nbinstat (1, "")
%!error<nbinstat: R and PS must not be complex.> nbinstat (i, 2)
%!error<nbinstat: R and PS must not be complex.> nbinstat (1, i)
%!error<nbinstat: R and PS must be of common size or scalars.> ...
%! nbinstat (ones (3), ones (2))
%!error<nbinstat: R and PS must be of common size or scalars.> ...
%! nbinstat (ones (2), ones (3))

## Output validation tests
%!test
%! r = 1:4;
%! ps = 0.2:0.2:0.8;
%! [m, v] = nbinstat (r, ps);
%! expected_m = [ 4.0000, 3.0000, 2.0000, 1.0000];
%! expected_v = [20.0000, 7.5000, 3.3333, 1.2500];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
%!test
%! r = 1:4;
%! [m, v] = nbinstat (r, 0.5);
%! expected_m = [1, 2, 3, 4];
%! expected_v = [2, 4, 6, 8];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
