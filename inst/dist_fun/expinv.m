## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} expinv (@var{p}, @var{mu})
##
## Inverse of the exponential cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the exponential distribution with mean @var{mu}.  The size of
## @var{x} is the common size of @var{p} and @var{mu}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## @seealso{expcdf, exppdf, exprnd, expfit, explike, expstat}
## @end deftypefn

function x = expinv (p, mu)

  if (nargin != 2)
    print_usage ();
  endif

  if (! isscalar (p) || ! isscalar (mu))
    [retval, p, mu] = common_size (p, mu);
    if (retval > 0)
      error ("expinv: P and MU must be of common size or scalars.");
    endif
  endif

  if (iscomplex (p) || iscomplex (mu))
    error ("expinv: P and MU must not be complex.");
  endif

  if (! isscalar (p))
    sz = size (p);
  else
    sz = size (mu);
  endif

  if (iscomplex (p) || iscomplex (mu))
    error ("expinv: P and MU must not be complex.");
  endif

  if (isa (p, "single") || isa (mu, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  k = (p == 1) & (mu > 0);
  x(k) = Inf;

  k = (p >= 0) & (p < 1) & (mu > 0);
  if (isscalar (mu))
    x(k) = - mu * log (1 - p(k));
  else
    x(k) = - mu(k) .* log (1 - p(k));
  endif

endfunction


%!shared p
%! p = [-1 0 0.3934693402873666 1 2];
%!assert (expinv (p, 2*ones (1,5)), [NaN 0 1 Inf NaN], eps)
%!assert (expinv (p, 2), [NaN 0 1 Inf NaN], eps)
%!assert (expinv (p, 2*[1 0 NaN 1 1]), [NaN NaN NaN Inf NaN], eps)
%!assert (expinv ([p(1:2) NaN p(4:5)], 2), [NaN 0 NaN Inf NaN], eps)

## Test class of input preserved
%!assert (expinv ([p, NaN], 2), [NaN 0 1 Inf NaN NaN], eps)
%!assert (expinv (single ([p, NaN]), 2), single ([NaN 0 1 Inf NaN NaN]), eps)
%!assert (expinv ([p, NaN], single (2)), single ([NaN 0 1 Inf NaN NaN]), eps)

## Test input validation
%!error expinv ()
%!error expinv (1)
%!error expinv (1,2,3)
%!error expinv (ones (3), ones (2))
%!error expinv (ones (2), ones (3))
%!error expinv (i, 2)
%!error expinv (2, i)
