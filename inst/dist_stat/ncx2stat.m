## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} ncx2stat (@var{df}, @var{lambda})
##
## Compute statistics for the noncentral chi-squared distribution.
##
## @code{[@var{m}, @var{v}] = ncx2stat (@var{df}, @var{lambda})} returns the
## mean and variance of the noncentral chi-squared distribution with @var{df}
## degrees of freedom and noncentrality parameter @var{lambda}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the noncentral chi-squared distribution can be
## found at @url{https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution}
##
## @seealso{ncx2cdf, ncx2inv, ncx2pdf, ncx2rnd}
## @end deftypefn

function [m, v] = ncx2stat (df, lambda)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("ncx2stat: function called with too few input arguments.");
  endif

  ## Check for DF and LAMBDA being numeric
  if (! (isnumeric (df) && isnumeric (lambda)))
    error ("ncx2stat: DF and LAMBDA must be numeric.");
  endif

  ## Check for DF and LAMBDA being real
  if (iscomplex (df) || iscomplex (lambda))
    error ("ncx2stat: DF and LAMBDA must not be complex.");
  endif

  ## Check for common size of DF and LAMBDA
  if (! isscalar (df) || ! isscalar (lambda))
    [retval, df, lambda] = common_size (df, lambda);
    if (retval > 0)
      error ("ncx2stat: DF and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Initialize mean and variance
  if (isa (df, "single") || isa (lambda, "single"))
    m = NaN (size (df), "single");
    v = m;
  else
    m = NaN (size (df));
    v = m;
  endif

  ## Compute mean and variance for valid parameter values.
  k = (df > 0 & lambda >= 0);
  if (any (k(:)))
    m(k) = lambda(k) + df(k);
    v(k) = 2 * (df(k) + 2 * (lambda(k)));
  endif

endfunction

## Input validation tests
%!error<ncx2stat: function called with too few input arguments.> ncx2stat ()
%!error<ncx2stat: function called with too few input arguments.> ncx2stat (1)
%!error<ncx2stat: DF and LAMBDA must be numeric.> ncx2stat ({}, 2)
%!error<ncx2stat: DF and LAMBDA must be numeric.> ncx2stat (1, "")
%!error<ncx2stat: DF and LAMBDA must not be complex.> ncx2stat (i, 2)
%!error<ncx2stat: DF and LAMBDA must not be complex.> ncx2stat (1, i)
%!error<ncx2stat: DF and LAMBDA must be of common size or scalars.> ...
%! ncx2stat (ones (3), ones (2))
%!error<ncx2stat: DF and LAMBDA must be of common size or scalars.> ...
%! ncx2stat (ones (2), ones (3))

## Output validation tests
%!shared df, d1
%! df = [2, 0, -1, 1, 4];
%! d1 = [1, NaN, 3, -1, 2];
%!assert (ncx2stat (df, d1), [3, NaN, NaN, NaN, 6]);
%!assert (ncx2stat ([df(1:2), df(4:5)], 1), [3, NaN, 2, 5]);
%!assert (ncx2stat ([df(1:2), df(4:5)], 3), [5, NaN, 4, 7]);
%!assert (ncx2stat ([df(1:2), df(4:5)], 2), [4, NaN, 3, 6]);
%!assert (ncx2stat (2, [d1(1), d1(3:5)]), [3, 5, NaN, 4]);
%!assert (ncx2stat (0, [d1(1), d1(3:5)]), [NaN, NaN, NaN, NaN]);
%!assert (ncx2stat (1, [d1(1), d1(3:5)]), [2, 4, NaN, 3]);
%!assert (ncx2stat (4, [d1(1), d1(3:5)]), [5, 7, NaN, 6]);
