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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} tstat (@var{df})
##
## Compute statistics of the Student's T distribution.
##
## @code{[@var{m}, @var{v}] = tstat (@var{df})} returns the mean and variance of
## the Student's T distribution with @var{df} degrees of freedom.
##
## The size of @var{m} (mean) and @var{v} (variance) is the same size of the
## input argument.
##
## Further information about the Student's T distribution can be found at
## @url{https://en.wikipedia.org/wiki/Student%27s_t-distribution}
##
## @seealso{tcdf, tinv, tpdf, trnd}
## @end deftypefn

function [m, v] = tstat (df)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("tstat: function called with too few input arguments.");
  endif

  ## Check for DF being numeric
  if (! isnumeric (df))
    error ("tstat: DF must be numeric.");
  endif

  ## Check for DF being real
  if (iscomplex (df))
    error ("tstat: DF must not be complex.");
  endif

  ## Calculate moments
  m = zeros (size (df));
  v = df ./ (df - 2);

  ## Continue argument check
  k = find (! (df > 1) | ! (df < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif
  k = find (! (df > 2) & (df < Inf));
  if (any (k))
    v(k) = Inf;
  endif

endfunction

## Input validation tests
%!error<tstat: function called with too few input arguments.> tstat ()
%!error<tstat: DF must be numeric.> tstat ({})
%!error<tstat: DF must be numeric.> tstat ("")
%!error<tstat: DF must not be complex.> tstat (i)

## Output validation tests
%!test
%! df = 3:8;
%! [m, v] = tstat (df);
%! expected_m = [0, 0, 0, 0, 0, 0];
%! expected_v = [3.0000, 2.0000, 1.6667, 1.5000, 1.4000, 1.3333];
%! assert (m, expected_m);
%! assert (v, expected_v, 0.001);
