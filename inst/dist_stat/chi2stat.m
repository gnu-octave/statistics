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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} chi2stat (@var{df})
##
## Compute statistics of the chi-squared distribution.
##
## @code{[@var{m}, @var{v}] = chi2stat (@var{df})} returns the mean and
## variance of the chi-squared distribution with @var{df} degrees of freedom.
##
## The size of @var{m} (mean) and @var{v} (variance) is the same size of the
## input argument.
##
## Further information about the chi-squared distribution can be found at
## @url{https://en.wikipedia.org/wiki/Chi-squared_distribution}
##
## @seealso{chi2cdf, chi2inv, chi2pdf, chi2rnd}
## @end deftypefn

function [m, v] = chi2stat (df)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("chi2stat: function called with too few input arguments.");
  endif

  ## Check for DF being numeric
  if (! isnumeric (df))
    error ("chi2stat: DF must be numeric.");
  endif

  ## Check for DF being real
  if (iscomplex (df))
    error ("chi2stat: DF must not be complex.");
  endif

  ## Calculate moments
  m = df;
  v = 2 .* df;

  ## Continue argument check
  k = find (! (df > 0) | ! (df < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<chi2stat: function called with too few input arguments.> chi2stat ()
%!error<chi2stat: DF must be numeric.> chi2stat ({})
%!error<chi2stat: DF must be numeric.> chi2stat ("")
%!error<chi2stat: DF must not be complex.> chi2stat (i)

## Output validation tests
%!test
%! df = 1:6;
%! [m, v] = chi2stat (df);
%! assert (m, df);
%! assert (v, [2, 4, 6, 8, 10, 12], 0.001);
