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
## @deftypefn  {statistics} @var{x} = wblinv (@var{p})
## @deftypefnx {statistics} @var{x} = wblinv (@var{p}, @var{lambda})
## @deftypefnx {statistics} @var{x} = wblinv (@var{p}, @var{lambda}, @var{k})
##
## Inverse of the Weibull cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Weibull distribution with parameters @var{lambda} and
## @var{k}.  The size of @var{x} is the common size of @var{p}, @var{lambda},
## and @var{k}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## Default values are @var{lambda} = 1, @var{k} = 1.
##
## @seealso{wblcdf, wblpdf, wblrnd, wblstat, wblplot}
## @end deftypefn

function x = wblinv (p, lambda = 1, k = 1)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  if (! isscalar (p) || ! isscalar (lambda) || ! isscalar (k))
    [retval, p, lambda, k] = common_size (p, lambda, k);
    if (retval > 0)
      error ("wblinv: X, LAMBDA, and K must be of common size or scalars.");
    endif
  endif

  if (iscomplex (p) || iscomplex (lambda) || iscomplex (k))
    error ("wblinv: X, LAMBDA, and K must not be complex.");
  endif

  if (isa (p, "single") || isa (lambda, "single") || isa (k, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ok = (lambda > 0) & (lambda < Inf) & (k > 0) & (k < Inf);

  pk = (p == 0) & ok;
  x(pk) = 0;

  pk = (p == 1) & ok;
  x(pk) = Inf;

  pk = (p > 0) & (p < 1) & ok;
  if (isscalar (lambda) && isscalar (k))
    x(pk) = lambda * (- log (1 - p(pk))) .^ (1 / k);
  else
    x(pk) = lambda(pk) .* (- log (1 - p(pk))) .^ (1 ./ k(pk));
  endif

endfunction


%!shared p
%! p = [-1 0 0.63212055882855778 1 2];
%!assert (wblinv (p, ones (1,5), ones (1,5)), [NaN 0 1 Inf NaN], eps)
%!assert (wblinv (p, 1, ones (1,5)), [NaN 0 1 Inf NaN], eps)
%!assert (wblinv (p, ones (1,5), 1), [NaN 0 1 Inf NaN], eps)
%!assert (wblinv (p, [1 -1 NaN Inf 1], 1), [NaN NaN NaN NaN NaN])
%!assert (wblinv (p, 1, [1 -1 NaN Inf 1]), [NaN NaN NaN NaN NaN])
%!assert (wblinv ([p(1:2) NaN p(4:5)], 1, 1), [NaN 0 NaN Inf NaN])

## Test class of input preserved
%!assert (wblinv ([p, NaN], 1, 1), [NaN 0 1 Inf NaN NaN], eps)
%!assert (wblinv (single ([p, NaN]), 1, 1), single ([NaN 0 1 Inf NaN NaN]), eps ("single"))
%!assert (wblinv ([p, NaN], single (1), 1), single ([NaN 0 1 Inf NaN NaN]), eps ("single"))
%!assert (wblinv ([p, NaN], 1, single (1)), single ([NaN 0 1 Inf NaN NaN]), eps ("single"))

## Test input validation
%!error wblinv ()
%!error wblinv (1,2,3,4)
%!error wblinv (ones (3), ones (2), ones (2))
%!error wblinv (ones (2), ones (3), ones (2))
%!error wblinv (ones (2), ones (2), ones (3))
%!error wblinv (i, 2, 2)
%!error wblinv (2, i, 2)
%!error wblinv (2, 2, i)
