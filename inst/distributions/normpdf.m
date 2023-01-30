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
## @deftypefn  {statistics} @var{y} = normpdf (@var{x})
## @deftypefnx {statistics} @var{y} = normpdf (@var{x}, @var{mu})
## @deftypefnx {statistics} @var{y} = (@var{x}, @var{mu}, @var{sigma})
##
## Normal probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the normal distribution with mean @var{mu} and standard
## deviation @var{sigma}.  The size of @var{p} is the common size of @var{p},
## @var{mu} and @var{sigma}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{sigma} = 1.
##
## @seealso{norminv, norminv, normrnd, normfit, normlike, normstat}
## @end deftypefn

function y = normpdf (x, mu = 0, sigma = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    error ("normpdf: invalid number of input arguments.");
  endif

  ## Check for common size of X, MU, and SIGMA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (sigma))
    [retval, x, mu, sigma] = common_size (x, mu, sigma);
    if (retval > 0)
      error ("normpdf: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("normpdf: X, MU, and SIGMA must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Compute normal PDF
  if (isscalar (mu) && isscalar (sigma))
    if (isfinite (mu) && (sigma > 0) && (sigma < Inf))
      y = stdnormal_pdf ((x - mu) / sigma) / sigma;
    else
      y = NaN (size (x), class (y));
    endif
  else
    k = isinf (mu) | !(sigma > 0) | !(sigma < Inf);
    y(k) = NaN;

    k = ! isinf (mu) & (sigma > 0) & (sigma < Inf);
    y(k) = stdnormal_pdf ((x(k) - mu(k)) ./ sigma(k)) ./ sigma(k);
  endif

endfunction


%!shared x,y
%! x = [-Inf 1 2 Inf];
%! y = 1/sqrt(2*pi)*exp (-(x-1).^2/2);
%!assert (normpdf (x, ones (1,4), ones (1,4)), y, eps)
%!assert (normpdf (x, 1, ones (1,4)), y, eps)
%!assert (normpdf (x, ones (1,4), 1), y, eps)
%!assert (normpdf (x, [0 -Inf NaN Inf], 1), [y(1) NaN NaN NaN], eps)
%!assert (normpdf (x, 1, [Inf NaN -1 0]), [NaN NaN NaN NaN], eps)
%!assert (normpdf ([x, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (normpdf (single ([x, NaN]), 1, 1), single ([y, NaN]), eps ("single"))
%!assert (normpdf ([x, NaN], single (1), 1), single ([y, NaN]), eps ("single"))
%!assert (normpdf ([x, NaN], 1, single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error normpdf ()
%!error normpdf (1,2)
%!error normpdf (1,2,3,4)
%!error normpdf (ones (3), ones (2), ones (2))
%!error normpdf (ones (2), ones (3), ones (2))
%!error normpdf (ones (2), ones (2), ones (3))
%!error normpdf (i, 2, 2)
%!error normpdf (2, i, 2)
%!error normpdf (2, 2, i)
