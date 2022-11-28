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
## @deftypefn {statistics} @var{y} = nctpdf (@var{x}, @var{df}, @var{delta})
##
## Noncentral Τ probability density function (pdf).
##
## @code{@var{y} = nctpdf (@var{x}, @var{df}, @var{delta})} returns
## the noncentral T distribution with @var{df} degrees of freedom and
## noncentrality parameter @var{delta}, at the values of @var{x}.
##
## The size of @var{y} is the common size of @var{df} and @var{delta}.  A scalar
## input functions as a constant matrix of the same size as the K_all inputs.
##
## @seealso{nctcdf, nctinv, nctrnd, nctstat, tpdf, pdf}
## @end deftypefn

function y = nctpdf (x, df, delta)

  ## Check for valid input arguments
  if (nargin <  3)
    error ("nctpdf: too few input arguments.");
  endif

  ## Check and fix size of input arguments
  [err, x, df, delta] = common_size (x, df, delta);
  if (err > 0)
    error ("nctpdf: input size mismatch.");
  endif

  ## Initialize Y
  if (isa (x, "single") || isa (df, "single") || isa (delta, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Find NaNs in input arguments (if any) and propagate them to p
  is_nan = isnan (x) | isnan (df) | isnan (delta);
  y(is_nan) = NaN;

  ## Force invalid parameter cases to NaN
  invalid = df <= 0 | ! isfinite (delta);
  y(invalid) = NaN;

  ## Use normal approximation for df > 1e6
  bigDF = df > 1e6 & ! is_nan & ! invalid;
  if (any (bigDF(:)))
    s = 1 - 1 ./ (4 * df);
    d = sqrt (1 + x .^ 2 ./ (2 * df));
    y(bigDF) = normpdf (x(bigDF) .* s(bigDF), delta(bigDF), d(bigDF));
  endif

  ## For negative x use left tail cdf
  x_neg = find ((x < 0) & isfinite (x) & df <= 1e6 & ! is_nan & ! invalid);
  if (any (x_neg))
    y(x_neg) = (df(x_neg) ./ x(x_neg)) .* ...
               (nctcdf (x(x_neg) .* sqrt ((df(x_neg) + 2) ./ df(x_neg)), ...
                        df(x_neg) + 2, delta(x_neg)) - ...
                nctcdf (x(x_neg), df(x_neg), delta(x_neg)));
  endif

  ## For positive x reflect about zero and use left tail cdf
  x_pos = find ((x > 0) & isfinite (x) & df <= 1e6 & ! is_nan & ! invalid);
  if (any (x_pos))
    y(x_pos) = (-df(x_pos) ./ x(x_pos)) .* ...
               (nctcdf (-x(x_pos) .* sqrt ((df(x_pos) + 2) ./ df(x_pos)), ...
                        df(x_pos) + 2, -delta(x_pos)) - ...
                nctcdf (-x(x_pos), df(x_pos), -delta(x_pos)));
  endif

  ## For x == 0 use power series
  xzero = find ((x == 0) & df <= 1e6 & ! is_nan & ! invalid);
  if (any (xzero))
    y(xzero) = exp (-0.5 * delta(xzero) .^ 2 - 0.5 * log (pi * df(xzero)) + ...
               gammaln (0.5 * (df(xzero) + 1)) - gammaln (0.5 * df(xzero)));
  endif

endfunction

## Input validation tests
%!error<nctpdf: too few input arguments.> y = nctpdf ();
%!error<nctpdf: too few input arguments.> y = nctpdf (2);
%!error<nctpdf: too few input arguments.> y = nctpdf (2, 4);
%!error<nctpdf: input size mismatch.> y = nctpdf (5, [4, 3], [3, 4, 5]);

## Output validation tests
%!shared x1, df, delta
%! x1 = [-Inf, 2, NaN, 4, Inf];
%! df = [2, 0, -1, 1, 4];
%! delta = [1, NaN, 3, -1, 2];
%!assert (nctpdf (x1, df, delta), [0, NaN, NaN, 0.00401787561306999, 0], 1e-14);
%!assert (nctpdf (x1, df, 1), [0, NaN, NaN, 0.0482312135423008, 0], 1e-14);
%!assert (nctpdf (x1, df, 3), [0, NaN, NaN, 0.1048493126401585, 0], 1e-14);
%!assert (nctpdf (x1, df, 2), [0, NaN, NaN, 0.08137377919890307, 0], 1e-14);
%!assert (nctpdf (x1, 3, delta), [0, NaN, NaN, 0.001185305171654381, 0], 1e-14);
%!assert (nctpdf (2, df, delta), [0.1791097459405861, NaN, NaN, ...
%!                             0.0146500727180389, 0.3082302682110299], 1e-14);
%!assert (nctpdf (4, df, delta), [0.04467929612254971, NaN, NaN, ...
%!                             0.00401787561306999, 0.0972086534042828], 1e-14);
