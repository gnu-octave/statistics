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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} unidstat (@var{df})
##
## Compute statistics of the discrete uniform cumulative distribution.
##
## @code{[@var{m}, @var{v}] = unidstat (@var{df})} returns the mean and variance
## of the discrete uniform cumulative distribution with parameter @var{N}, which
## corresponds to the maximum observable value and must be a positive natural
## number.
##
## The size of @var{m} (mean) and @var{v} (variance) is the same size of the
## input argument.
##
## Further information about the discrete uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Discrete_uniform_distribution}
##
## @seealso{unidcdf, unidinv, unidpdf, unidrnd, unidfit}
## @end deftypefn

function [m, v] = unidstat (N)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("unidstat: function called with too few input arguments.");
  endif

  ## Check for N being numeric
  if (! isnumeric (N))
    error ("unidstat: N must be numeric.");
  endif

  ## Check for N being real
  if (iscomplex (N))
    error ("unidstat: N must not be complex.");
  endif

  ## Calculate moments
  m = (N + 1) ./ 2;
  v = ((N .^ 2) - 1) ./ 12;

  ## Continue argument check
  k = find (! (N > 0) | ! (N < Inf) | ! (N == round (N)));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<unidstat: function called with too few input arguments.> unidstat ()
%!error<unidstat: N must be numeric.> unidstat ({})
%!error<unidstat: N must be numeric.> unidstat ("")
%!error<unidstat: N must not be complex.> unidstat (i)

## Output validation tests
%!test
%! N = 1:6;
%! [m, v] = unidstat (N);
%! expected_m = [1.0000, 1.5000, 2.0000, 2.5000, 3.0000, 3.5000];
%! expected_v = [0.0000, 0.2500, 0.6667, 1.2500, 2.0000, 2.9167];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
