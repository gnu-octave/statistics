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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} ncfstat (@var{df1}, @var{df1}, @var{lambda})
##
## Compute statistics for the noncentral @math{F}-distribution.
##
## @code{[@var{m}, @var{v}] = ncfstat (@var{df1}, @var{df1}, @var{lambda})}
## returns the mean and variance of the noncentral @math{F}-distribution with
## @var{df1} and @var{df2} degrees of freedom and noncentrality parameter
## @var{lambda}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## Further information about the noncentral @math{F}-distribution can be found
## at @url{https://en.wikipedia.org/wiki/Noncentral_F-distribution}
##
## @seealso{ncfcdf, ncfinv, ncfpdf, ncfrnd, fstat}
## @end deftypefn

function [m, v] = ncfstat (df1, df2, lambda)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("ncfstat: function called with too few input arguments.");
  endif

  ## Check for DF1, DF2, and LAMBDA being numeric
  if (! (isnumeric (df1) && isnumeric (df2) && isnumeric (lambda)))
    error ("ncfstat: DF1, DF2, and LAMBDA must be numeric.");
  endif

  ## Check for DF1, DF2, and LAMBDA being reals
  if (iscomplex (df1) || iscomplex (df2) || iscomplex (lambda))
    error ("ncfstat: DF1, DF2, and LAMBDA must not be complex.");
  endif

  ## Check for common size of DF1, DF2, and LAMBDA
  if (! isscalar (df1) || ! isscalar (df2) || ! isscalar (lambda))
    [retval, df1, df2, lambda] = common_size (df1, df2, lambda);
    if (retval > 0)
      error ("ncfstat: DF1, DF2, and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Initialize mean and variance
  if (isa (df1, "single") || isa (df2, "single") || isa (lambda, "single"))
    m = zeros (size (df1), "single");
    v = m;
  else
    m = zeros (size (df1));
    v = m;
  endif

  ## Return NaNs for invalid df2 parameters
  m(df2 <= 2) = NaN;
  v(df2 <= 4) = NaN;

  ## Compute mean and variance for valid parameter values.
  k = (df2 > 2);
  if (any (k(:)))
    m(k) = df2(k) .* (df1(k) + lambda(k)) ./ (df1(k) .* (df2(k) - 2));
  endif
  k = (df2 > 4);
  if (any (k(:)))
    df1_idx = df1(k) + lambda(k);
    df2_idx = df2(k) - 2;
    df1_df2 = (df2(k) ./ df1(k)) .^ 2;
    v(k) = 2 * df1_df2 .* (df1_idx .^ 2 + (df1_idx + lambda(k)) .* ...
                           df2_idx) ./ ((df2(k) - 4) .* df2_idx .^ 2);
  endif

endfunction

## Input validation tests
%!error<ncfstat: function called with too few input arguments.> ncfstat ()
%!error<ncfstat: function called with too few input arguments.> ncfstat (1)
%!error<ncfstat: function called with too few input arguments.> ncfstat (1, 2)
%!error<ncfstat: DF1, DF2, and LAMBDA must be numeric.> ncfstat ({}, 2, 3)
%!error<ncfstat: DF1, DF2, and LAMBDA must be numeric.> ncfstat (1, "", 3)
%!error<ncfstat: DF1, DF2, and LAMBDA must be numeric.> ncfstat (1, 2, "")
%!error<ncfstat: DF1, DF2, and LAMBDA must not be complex.> ncfstat (i, 2, 3)
%!error<ncfstat: DF1, DF2, and LAMBDA must not be complex.> ncfstat (1, i, 3)
%!error<ncfstat: DF1, DF2, and LAMBDA must not be complex.> ncfstat (1, 2, i)
%!error<ncfstat: DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfstat (ones (3), ones (2), 3)
%!error<ncfstat: DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfstat (ones (2), 2, ones (3))
%!error<ncfstat: DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfstat (1, ones (2), ones (3))

## Output validation tests
%!shared df1, df2, lambda
%! df1 = [2, 0, -1, 1, 4, 5];
%! df2 = [2, 4, -1, 5, 6, 7];
%! lambda = [1, NaN, 3, 0, 2, -1];
%!assert (ncfstat (df1, df2, lambda), [NaN, NaN, NaN, 1.6667, 2.25, 1.12], 1e-4);
%!assert (ncfstat (df1(4:6), df2(4:6), 1), [3.3333, 1.8750, 1.6800], 1e-4);
%!assert (ncfstat (df1(4:6), df2(4:6), 2), [5.0000, 2.2500, 1.9600], 1e-4);
%!assert (ncfstat (df1(4:6), df2(4:6), 3), [6.6667, 2.6250, 2.2400], 1e-4);
%!assert (ncfstat (2, [df2(1), df2(4:6)], 5), [NaN,5.8333,5.2500,4.9000], 1e-4);
%!assert (ncfstat (0, [df2(1), df2(4:6)], 5), [NaN, Inf, Inf, Inf]);
%!assert (ncfstat (1, [df2(1), df2(4:6)], 5), [NaN, 10, 9, 8.4], 1e-14);
%!assert (ncfstat (4, [df2(1), df2(4:6)], 5), [NaN, 3.75, 3.375, 3.15], 1e-14);
