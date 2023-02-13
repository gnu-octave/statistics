## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{p} =} logistic_cdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} logistic_cdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} logistic_cdf (@var{x}, @var{mu}, @var{scale})
##
## Logistic cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the logistic distribution with mean @var{mu} and scale
## parameter @var{scale}.  The size of @var{p} is the common size of @var{x},
## @var{mu}, and @var{scale}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{beta} = 1.  Both parameters must be
## reals and @var{beta} > 0.  For @var{beta} <= 0, NaN is returned.
##
## @seealso{logistic_inv, logistic_pdf, logistic_rnd}
## @end deftypefn

function p = logistic_cdf (x, mu = 0, scale = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of X, MU, and SCALE
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(scale))
    [retval, x, mu, scale] = ...
        common_size (x, mu, scale);
    if (retval > 0)
      error (strcat (["logistic_cdf: X, MU, and SCALE must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, MU, and SCALE being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (scale))
    error ("logistic_cdf: X, MU, and SCALE must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (mu, "single") || isa (scale, "single"));
    p = NaN (size (x), "single");
  else
    p = NaN (size (x));
  endif

  ## Compute logistic CDF
  k1 = (x == -Inf) & (scale > 0);
  p(k1) = 0;

  k2 = (x == Inf) & (scale > 0);
  p(k2) = 1;

  k = ! k1 & ! k2 & (scale > 0);
  p(k) = 1 ./ (1 + exp (- (x(k) - mu(k)) ./ scale(k)));

endfunction


%!shared x,y
%! x = [-Inf -log(3) 0 log(3) Inf];
%! y = [0, 1/4, 1/2, 3/4, 1];
%!assert (logistic_cdf ([x, NaN]), [y, NaN], eps)
%!assert (logistic_cdf (x, 0, [-2, -1, 0, 1, 2]), [nan(1, 3), 0.75, 1])

## Test class of input preserved
%!assert (logistic_cdf (single ([x, NaN])), single ([y, NaN]), eps ("single"))

## Test input validation
%!error logistic_cdf ()
%!error logistic_cdf (1, 2, 3, 4)
%!error<logistic_cdf: X, MU, and SCALE must be of common size or scalars.> ...
%! logistic_cdf (1, ones (2), ones (3))
%!error<logistic_cdf: X, MU, and SCALE must be of common size or scalars.> ...
%! logistic_cdf (ones (2), 1, ones (3))
%!error<logistic_cdf: X, MU, and SCALE must be of common size or scalars.> ...
%! logistic_cdf (ones (2), ones (3), 1)
%!error<logistic_cdf: X, MU, and SCALE must not be complex.> ...
%! logistic_cdf (i, 2, 3)
%!error<logistic_cdf: X, MU, and SCALE must not be complex.> ...
%! logistic_cdf (1, i, 3)
%!error<logistic_cdf: X, MU, and SCALE must not be complex.> ...
%! logistic_cdf (1, 2, i)
