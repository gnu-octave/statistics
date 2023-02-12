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
## @deftypefn  {statistics} {@var{y} =} laplace_pdf (@var{x})
## @deftypefnx {statistics} {@var{y} =} laplace_pdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{y} =} laplace_pdf (@var{x}, @var{mu}, @var{sigma})
##
## Lognormal probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the lognormal distribution with parameters @var{mu} and
## @var{sigma}.  The size of @var{p} is the common size of @var{p}, @var{mu},
## and @var{sigma}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## If a random variable follows this distribution, its logarithm is normally
## distributed with mean @var{mu} and standard deviation @var{sigma}.
##
## Default values are @var{mu} = 0, @var{sigma} = 1. Both parameters must be
## reals and @var{sigma} > 0.  For @var{sigma} <= 0, NaN is returned.
##
## @seealso{logncdf, logninv, lognrnd, lognstat}
## @end deftypefn

function y = lognpdf (x, mu = 0, sigma = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of P, MU, and SIGMA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (sigma))
    [retval, x, mu, sigma] = common_size (x, mu, sigma);
    if (retval > 0)
      error ("lognpdf: X, MU, and SIGMA must be of common size or scalars");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("lognpdf: X, MU, and SIGMA must not be complex");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Compute lognormal PDF
  k = isnan (x) | !(sigma > 0) | !(sigma < Inf);
  y(k) = NaN;

  k = (x > 0) & (x < Inf) & (sigma > 0) & (sigma < Inf);
  if (isscalar (mu) && isscalar (sigma))
    y(k) = normpdf (log (x(k)), mu, sigma) ./ x(k);
  else
    y(k) = normpdf (log (x(k)), mu(k), sigma(k)) ./ x(k);
  endif

endfunction


%!shared x,y
%! x = [-1 0 e Inf];
%! y = [0, 0, 1/(e*sqrt(2*pi)) * exp(-1/2), 0];
%!assert (lognpdf (x, zeros (1,4), ones (1,4)), y, eps)
%!assert (lognpdf (x, 0, ones (1,4)), y, eps)
%!assert (lognpdf (x, zeros (1,4), 1), y, eps)
%!assert (lognpdf (x, [0 1 NaN 0], 1), [0 0 NaN y(4)], eps)
%!assert (lognpdf (x, 0, [0 NaN Inf 1]), [NaN NaN NaN y(4)], eps)
%!assert (lognpdf ([x, NaN], 0, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (lognpdf (single ([x, NaN]), 0, 1), single ([y, NaN]), eps ("single"))
%!assert (lognpdf ([x, NaN], single (0), 1), single ([y, NaN]), eps ("single"))
%!assert (lognpdf ([x, NaN], 0, single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error lognpdf ()
%!error lognpdf (1,2,3,4)
%!error lognpdf (ones (3), ones (2), ones (2))
%!error lognpdf (ones (2), ones (3), ones (2))
%!error lognpdf (ones (2), ones (2), ones (3))
%!error lognpdf (i, 2, 2)
%!error lognpdf (2, i, 2)
%!error lognpdf (2, 2, i)
