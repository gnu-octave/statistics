## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{x} =} nakacdf (@var{x}, @var{m}, @var{w})
##
## Inverse of the Nakagami cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Nakagami distribution with shape parameter @var{m} and
## scale parameter @var{w}.  The size of @var{p} is the common size of @var{x},
## @var{m}, and @var{w}.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## @seealso{nakacdf, nakapdf, nakarnd}
## @end deftypefn

function x = nakainv (p, m, w)

  ## Check for valid number of input arguments
  if (nargin != 3)
    print_usage ();
  endif

  ## Check for common size of P, M, and W
  if (! isscalar (p) || ! isscalar (m) || ! isscalar (w))
    [retval, p, m, w] = common_size (p, m, w);
    if (retval > 0)
      error ("nakainv: P, M and W must be of common size or scalars.");
    endif
  endif

  ## Check for P, M, and W being reals
  if (iscomplex (p) || iscomplex (m) || iscomplex (w))
    error ("nakainv: P, M, and W must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (p, "single") || isa (m, "single") || isa (w, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  ## Compute Nakagami iCDF
  k = isnan (p) | ! (0 <= p) | ! (p <= 1) | ! (-Inf < m) | ! (m < Inf) ...
    | ! (0 < w) | ! (w < Inf);
  x(k) = NaN;

  k = (p == 1) & (-Inf < m) & (m < Inf) & (0 < w) & (w < Inf);
  x(k) = Inf;

  k = (0 < p) & (p < 1) & (0 < m) & (m < Inf) & (0 < w) & (w < Inf);
  if (isscalar (m) && isscalar(w))
    m_gamma = m;
    w_gamma = w/m;
    x(k) = gaminv(p(k), m_gamma, w_gamma);
    x(k) = sqrt(x(k));
  else
    m_gamma = m;
    w_gamma = w./m;
    x(k) = gaminv(p(k), m_gamma(k), w_gamma(k));
    x(k) = sqrt(x(k));
  endif

endfunction


%!shared p,y
%! p = [-Inf, -1, 0, 1/2, 1, 2, Inf];
%! y = [NaN, NaN, 0, 0.83255461115769769, Inf, NaN, NaN];
%!assert (nakainv (p, ones (1,7), ones (1,7)), y, eps)
%!assert (nakainv (p, 1, 1), y, eps)
%!assert (nakainv (p, [1, 1, 1, NaN, 1, 1, 1], 1), [y(1:3), NaN, y(5:7)], eps)
%!assert (nakainv (p, 1, [1, 1, 1, NaN, 1, 1, 1]), [y(1:3), NaN, y(5:7)], eps)
%!assert (nakainv ([p, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (nakainv (single ([p, NaN]), 1, 1), single ([y, NaN]))
%!assert (nakainv ([p, NaN], single (1), 1), single ([y, NaN]))
%!assert (nakainv ([p, NaN], 1, single (1)), single ([y, NaN]))

## Test input validation
%!error nakainv ()
%!error nakainv (1)
%!error nakainv (1,2)
%!error nakainv (1,2,3,4)
%!error nakainv (ones (3), ones (2), ones(2))
%!error nakainv (ones (2), ones (3), ones(2))
%!error nakainv (ones (2), ones (2), ones(3))
%!error nakainv (i, 2, 2)
%!error nakainv (2, i, 2)
%!error nakainv (2, 2, i)

