## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} nctpdf (@var{x}, @var{df}, @var{mu})
##
## Noncentral Τ probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the noncentral T distribution with @var{df} degrees of freedom and
## noncentrality parameter @var{mu}.  The size of @var{y} is the common size of
## @var{x}, @var{df}, and @var{mu}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Further information about the noncentral T distribution can be found at
## @url{https://en.wikipedia.org/wiki/Noncentral_t-distribution}
##
## @seealso{nctcdf, nctinv, nctrnd, nctstat, tpdf}
## @end deftypefn

function y = nctpdf (x, df, mu)

  ## Check for valid number of input arguments
  if (nargin <  3)
    error ("nctpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, DF, and MU
  [err, x, df, mu] = common_size (x, df, mu);
  if (err > 0)
    error ("nctpdf: X, DF, and MU must be of common size or scalars.");
  endif

  ## Check for X, DF, and MU being reals
  if (iscomplex (x) || iscomplex (df) || iscomplex (mu))
    error ("nctpdf: X, DF, and MU must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (df, "single") || isa (mu, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Find NaNs in input arguments (if any) and propagate them to p
  is_nan = isnan (x) | isnan (df) | isnan (mu);
  y(is_nan) = NaN;

  ## Force invalid parameter cases to NaN
  invalid = df <= 0 | ! isfinite (mu);
  y(invalid) = NaN;

  ## Use normal approximation for df > 1e6
  bigDF = df > 1e6 & ! is_nan & ! invalid;
  if (any (bigDF(:)))
    s = 1 - 1 ./ (4 * df);
    d = sqrt (1 + x .^ 2 ./ (2 * df));
    y(bigDF) = normpdf (x(bigDF) .* s(bigDF), mu(bigDF), d(bigDF));
  endif

  ## For negative x use left tail cdf
  x_neg = find ((x < 0) & isfinite (x) & df <= 1e6 & ! is_nan & ! invalid);
  if (any (x_neg))
    y(x_neg) = (df(x_neg) ./ x(x_neg)) .* ...
               (nctcdf (x(x_neg) .* sqrt ((df(x_neg) + 2) ./ df(x_neg)), ...
                        df(x_neg) + 2, mu(x_neg)) - ...
                nctcdf (x(x_neg), df(x_neg), mu(x_neg)));
  endif

  ## For positive x reflect about zero and use left tail cdf
  x_pos = find ((x > 0) & isfinite (x) & df <= 1e6 & ! is_nan & ! invalid);
  if (any (x_pos))
    y(x_pos) = (-df(x_pos) ./ x(x_pos)) .* ...
               (nctcdf (-x(x_pos) .* sqrt ((df(x_pos) + 2) ./ df(x_pos)), ...
                        df(x_pos) + 2, -mu(x_pos)) - ...
                nctcdf (-x(x_pos), df(x_pos), -mu(x_pos)));
  endif

  ## For x == 0 use power series
  xzero = find ((x == 0) & df <= 1e6 & ! is_nan & ! invalid);
  if (any (xzero))
    y(xzero) = exp (-0.5 * mu(xzero) .^ 2 - 0.5 * log (pi * df(xzero)) + ...
               gammaln (0.5 * (df(xzero) + 1)) - gammaln (0.5 * df(xzero)));
  endif

endfunction

%!demo
%! ## Plot various PDFs from the noncentral T distribution
%! x = -5:0.01:10;
%! y1 = nctpdf (x, 1, 0);
%! y2 = nctpdf (x, 4, 0);
%! y3 = nctpdf (x, 1, 2);
%! y4 = nctpdf (x, 4, 2);
%! plot (x, y1, "-r", x, y2, "-g", x, y3, "-k", x, y4, "-m")
%! grid on
%! xlim ([-5, 10])
%! ylim ([0, 0.4])
%! legend ({"df = 1, μ = 0", "df = 4, μ = 0", ...
%!          "df = 1, μ = 2", "df = 4, μ = 2"}, "location", "northeast")
%! title ("Noncentral T PDF")
%! xlabel ("values in x")
%! ylabel ("density")

%!demo
%! ## Compare the noncentral T PDF with MU = 1 to the T PDF
%! ## with the same number of degrees of freedom (10).
%!
%! x = -5:0.1:5;
%! y1 = nctpdf (x, 10, 1);
%! y2 = tpdf (x, 10);
%! plot (x, y1, "-", x, y2, "-");
%! grid on
%! xlim ([-5, 5])
%! ylim ([0, 0.4])
%! legend ({"Noncentral χ^2(4,2)", "χ^2(4)"}, "location", "northwest")
%! title ("Noncentral T vs T PDFs")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x1, df, mu
%! x1 = [-Inf, 2, NaN, 4, Inf];
%! df = [2, 0, -1, 1, 4];
%! mu = [1, NaN, 3, -1, 2];
%!assert (nctpdf (x1, df, mu), [0, NaN, NaN, 0.00401787561306999, 0], 1e-14);
%!assert (nctpdf (x1, df, 1), [0, NaN, NaN, 0.0482312135423008, 0], 1e-14);
%!assert (nctpdf (x1, df, 3), [0, NaN, NaN, 0.1048493126401585, 0], 1e-14);
%!assert (nctpdf (x1, df, 2), [0, NaN, NaN, 0.08137377919890307, 0], 1e-14);
%!assert (nctpdf (x1, 3, mu), [0, NaN, NaN, 0.001185305171654381, 0], 1e-14);
%!assert (nctpdf (2, df, mu), [0.1791097459405861, NaN, NaN, ...
%!                             0.0146500727180389, 0.3082302682110299], 1e-14);
%!assert (nctpdf (4, df, mu), [0.04467929612254971, NaN, NaN, ...
%!                             0.00401787561306999, 0.0972086534042828], 1e-14);

## Test input validation
%!error<nctpdf: function called with too few input arguments.> nctpdf ()
%!error<nctpdf: function called with too few input arguments.> nctpdf (1)
%!error<nctpdf: function called with too few input arguments.> nctpdf (1, 2)
%!error<nctpdf: X, DF, and MU must be of common size or scalars.> ...
%! nctpdf (ones (3), ones (2), ones (2))
%!error<nctpdf: X, DF, and MU must be of common size or scalars.> ...
%! nctpdf (ones (2), ones (3), ones (2))
%!error<nctpdf: X, DF, and MU must be of common size or scalars.> ...
%! nctpdf (ones (2), ones (2), ones (3))
%!error<nctpdf: X, DF, and MU must not be complex.> nctpdf (i, 2, 2)
%!error<nctpdf: X, DF, and MU must not be complex.> nctpdf (2, i, 2)
%!error<nctpdf: X, DF, and MU must not be complex.> nctpdf (2, 2, i)
