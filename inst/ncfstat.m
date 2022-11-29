## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} [@var{m}, @var{v}] = ncfstat (@var{df1}, @var{df1}, @var{delta})
##
## Mean and variance for the noncentral F distribution.
##
## @code{[@var{m}, @var{v}] = ncfstat (@var{df1}, @var{df1}, @var{delta})}
## returns the mean and variance of the noncentral F distribution with @var{df1}
## and @var{df2} degrees of freedom and noncentrality parameter @var{delta}.
##
## The size of @var{m} and @var{v} is the common size of the input arguments.
## Scalar input arguments @var{df1}, @var{df2}, and @var{delta} are regarded as
## constant matrices of the same size as the other input.
##
## @seealso{ncfcdf, ncfinv, ncfpdf, ncfrnd}
## @end deftypefn

function [m, v] = ncfstat (df1, df2, delta)

  ## Check for valid input arguments
  if (nargin <  3)
    error ("ncfstat: too few input arguments.");
  endif

  ## Check and fix size of input arguments
  [err, df1, df2, delta] = common_size (df1, df2, delta);
  if (err > 0)
    error ("ncfstat: input size mismatch.");
  endif

  ## Initialize mean and variance
  if (isa (df1, "single") || isa (df2, "single") || isa (delta, "single"))
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
    m(k) = df2(k) .* (df1(k) + delta(k)) ./ (df1(k) .* (df2(k) - 2));
  endif
  k = (df2 > 4);
  if (any (k(:)))
    df1_idx = df1(k) + delta(k);
    df2_idx = df2(k) - 2;
    df1_df2 = (df2(k) ./ df1(k)) .^ 2;
    v(k) = 2 * df1_df2 .* (df1_idx .^ 2 + (df1_idx + delta(k)) .* ...
                           df2_idx) ./ ((df2(k) - 4) .* df2_idx .^ 2);
  endif

endfunction

## Input validation tests
%!error<ncfstat: too few input arguments.> p = ncfstat ();
%!error<ncfstat: too few input arguments.> p = ncfstat (1);
%!error<ncfstat: too few input arguments.> p = ncfstat (1, 2);
%!error<ncfstat: input size mismatch.> p = ncfstat (5, [4, 3], [3, 4, 5]);

## Output validation tests
%!shared df1, df2, delta
%! df1 = [2, 0, -1, 1, 4, 5];
%! df2 = [2, 4, -1, 5, 6, 7];
%! delta = [1, NaN, 3, 0, 2, -1];
%!assert (ncfstat (df1, df2, delta), [NaN, NaN, NaN, 1.6667, 2.25, 1.12], 1e-4);
%!assert (ncfstat (df1(4:6), df2(4:6), 1), [3.3333, 1.8750, 1.6800], 1e-4);
%!assert (ncfstat (df1(4:6), df2(4:6), 2), [5.0000, 2.2500, 1.9600], 1e-4);
%!assert (ncfstat (df1(4:6), df2(4:6), 3), [6.6667, 2.6250, 2.2400], 1e-4);
%!assert (ncfstat (2, [df2(1), df2(4:6)], 5), [NaN,5.8333,5.2500,4.9000], 1e-4);
%!assert (ncfstat (0, [df2(1), df2(4:6)], 5), [NaN, Inf, Inf, Inf]);
%!assert (ncfstat (1, [df2(1), df2(4:6)], 5), [NaN, 10, 9, 8.4], 1e-14);
%!assert (ncfstat (4, [df2(1), df2(4:6)], 5), [NaN, 3.75, 3.375, 3.15], 1e-14);
