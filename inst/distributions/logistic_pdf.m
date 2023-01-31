## Copyright (C) 1995-2017 Kurt Hornik
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
## @deftypefn  {statistics} @var{y} = logistic_pdf (@var{x})
## @deftypefnx {statistics} @var{y} = logistic_pdf (@var{x}, @var{mu})
## @deftypefnx {statistics} @var{y} = logistic_pdf (@var{x}, @var{mu}, @var{scale})
##
## Logistic probability density function (PDF).
##
## For each element of @var{x}, compute the PDF at @var{x} of the logistic
## distribution with mean @var{mu} and scale parameter @var{scale}.  The size of
## @var{y} is the common size of @var{x}, @var{mu}, and @var{scale}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{scale} = 1.  Both parameters must be
## reals and @var{beta} > 0.  For @var{beta} <= 0, NaN is returned.
##
## @seealso{logistic_cdf, logistic_inv, logistic_rnd}
## @end deftypefn

function y = logistic_pdf (x, mu = 0, scale = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of P, MU, and SCALE
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(scale))
    [retval, x, mu, scale] = common_size (x, mu, scale);
    if (retval > 0)
      error (strcat (["logistic_pdf: X, MU, and SCALE must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, MU, and SCALE being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (scale))
    error ("logistic_pdf: X, MU, and SCALE must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (mu, "single") || isa (scale, "single"));
    y = NaN (size (x), "single");
  else
    y = NaN (size (x));
  endif

  ## Compute logistic PDF
  k1 = ((x == -Inf) & (scale > 0)) | ((x == Inf) & (scale > 0));
  y(k1) = 0;

  k = ! k1 & (scale > 0);
  y(k) = (1 ./ (4 .* scale(k))) .* ...
         (sech ((x(k) - mu(k)) ./ (2 .* scale(k))) .^ 2);

endfunction


%!shared x, y
%! x = [-Inf -log(4) 0 log(4) Inf];
%! y = [0, 0.16, 1/4, 0.16, 0];
%!assert (logistic_pdf ([x, NaN]), [y, NaN], eps)
%!assert (logistic_pdf (x, 0, [-2, -1, 0, 1, 2]), [nan(1, 3), y([4:5])], eps)

## Test class of input preserved
%!assert (logistic_pdf (single ([x, NaN])), single ([y, NaN]), eps ("single"))

## Test input validation
%!error logistic_pdf ()
%!error logistic_pdf (1, 2, 3, 4)
%!error<logistic_pdf: X, MU, and SCALE must be of common size or scalars.> ...
%! logistic_pdf (1, ones (2), ones (3))
%!error<logistic_pdf: X, MU, and SCALE must be of common size or scalars.> ...
%! logistic_pdf (ones (2), 1, ones (3))
%!error<logistic_pdf: X, MU, and SCALE must be of common size or scalars.> ...
%! logistic_pdf (ones (2), ones (3), 1)
%!error<logistic_pdf: X, MU, and SCALE must not be complex.> ...
%! logistic_pdf (i, 2, 3)
%!error<logistic_pdf: X, MU, and SCALE must not be complex.> ...
%! logistic_pdf (1, i, 3)
%!error<logistic_pdf: X, MU, and SCALE must not be complex.> ...
%! logistic_pdf (1, 2, i)
