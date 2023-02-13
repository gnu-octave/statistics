## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1997-2016 Kurt Hornik
## Copyright (C) 2022 Nicholas R. Jankowski
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
## @deftypefn  {statistics} {@var{p} =} hygecdf (@var{x}, @var{t}, @var{m}, @var{n})
## @deftypefnx {statistics} {@var{p} =} hygecdf (@var{x}, @var{t}, @var{m}, @var{n}, "upper")
##
## Hypergeometric cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the hypergeometric distribution with parameters
## @var{t}, @var{k}, and @var{n}.  The size of @var{p} is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## This is the probability of obtaining not more than @var{x} marked items
## when randomly drawing a sample of size @var{n} without replacement from a
## population of total size @var{t} containing @var{m} marked items.
##
## The parameters @var{t}, @var{m}, and @var{n} must be positive integers
## with @var{m} and @var{n} not greater than @var{t}.
##
## @seealso{hygeinv, hygepdf, hygernd, hygestat}
## @end deftypefn

function p = hygecdf (x, t, m, n, uflag)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("hygecdf: too few input arguments.");
  endif

  ## Check for common size of X, T, M, and N
  if (! isscalar (x) || ! isscalar (t) || ! isscalar (m) || ! isscalar (n))
    [retval, x, t, m, n] = common_size (x, t, m, n);
    if (retval > 0)
      error ("hygecdf: X, T, M, and N must be of common size or scalars.");
    endif
  endif

  ## Check for X, T, M, and N being reals
  if (iscomplex (x) || iscomplex (t) || iscomplex (m) || iscomplex (n))
    error ("hygecdf: X, T, M, and N must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (t, "single")
      || isa (m, "single") || isa (n, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Check for "upper" flag
  if (nargin > 4 && strcmpi (uflag, "upper"))
    x = n - floor(x) - 1;
    m = t - m;
  elseif (nargin > 4 && ! strcmpi (uflag, "upper"))
    error ("hygecdf: invalid argument for upper tail.");
  endif

  ## Force 1 where required
  k1 = (x >= n | x >= m);
  p(k1) = 1;

  ## Force NaNs where required
  k2 = (isnan (x) | isnan (t) | isnan (m) | isnan (n) | ...
        t < 0 | m < 0 | n < 0 | round (t) != t | round (m) != m | ...
        round (n) != n | n > t | m > t);
  p(k2) = NaN;

  ## Get values for which P = 0
  k3 = (t - m - n + x + 1 <= 0 | x < 0);

  ok = ! (k1 | k2 | k3);

  ## Compute P
  if any(ok(:))
    ## For improved accuracy, compute the upper tail 1-p instead
    ## of the lower tail pfor x values that are larger than the mean
    lo = (x <= m .* n ./ t);
    ok_lo = ok & lo;
    if (any (ok_lo(:)))
      p(ok_lo) = localPDF (floor (x(ok_lo)), t(ok_lo), m(ok_lo), n(ok_lo));
    endif
    ok_hi = ok & ! lo;
    if (any (ok_hi(:)))
      p(ok_hi) = 1 - localPDF (n(ok_hi) - floor (x(ok_hi)) - 1, ...
                               t(ok_hi), t(ok_hi) - m(ok_hi), n(ok_hi));
    endif
  endif

endfunction

function p = localPDF (x, t, m, n)
  HPDF = hygepdf (x, t, m, n);

  ## Compute hygecdf(x,t,m,n)/hygepdf(x,t,m,n) with a series
  ## whose terms can be computed recursively, backwards.
  xmax = max (x(:));
  ybig = repmat((0:xmax)', 1, length (x));
  xbig = repmat(x(:)', xmax + 1, 1);
  mbig = repmat(t(:)', xmax + 1, 1);
  kbig = repmat(m(:)', xmax + 1, 1);
  nbig = repmat(n(:)', xmax + 1, 1);

  terms = ((ybig+1) .* (mbig-kbig-nbig+ybig+1)) ./ ((nbig-ybig) .* (kbig-ybig));
  terms(ybig >= xbig) = 1;
  terms = flip (cumprod (flip (terms)));

  terms(ybig > xbig) = 0;
  ratio = sum(terms,1);
  ratio = reshape(ratio,size(x));
  p = ratio.*HPDF;
  ## Correct round-off errors
  p(p > 1) = 1;
endfunction

%!shared x,y
%! x = [-1 0 1 2 3];
%! y = [0 1/6 5/6 1 1];
%!assert (hygecdf (x, 4*ones (1,5), 2, 2), y, 5*eps)
%!assert (hygecdf (x, 4, 2*ones (1,5), 2), y, 5*eps)
%!assert (hygecdf (x, 4, 2, 2*ones (1,5)), y, 5*eps)
%!assert (hygecdf (x, 4*[1 -1 NaN 1.1 1], 2, 2), [y(1) NaN NaN NaN y(5)], 5*eps)
%!assert (hygecdf (x, 4*[1 -1 NaN 1.1 1], 2, 2, "upper"), ...
%! [y(5) NaN NaN NaN y(1)], 5*eps)
%!assert (hygecdf (x, 4, 2*[1 -1 NaN 1.1 1], 2), [y(1) NaN NaN NaN y(5)], 5*eps)
%!assert (hygecdf (x, 4, 2*[1 -1 NaN 1.1 1], 2, "upper"), ...
%! [y(5) NaN NaN NaN y(1)], 5*eps)
%!assert (hygecdf (x, 4, 5, 2), [NaN NaN NaN NaN NaN])
%!assert (hygecdf (x, 4, 2, 2*[1 -1 NaN 1.1 1]), [y(1) NaN NaN NaN y(5)], 5*eps)
%!assert (hygecdf (x, 4, 2, 2*[1 -1 NaN 1.1 1], "upper"), ...
%! [y(5) NaN NaN NaN y(1)], 5*eps)
%!assert (hygecdf (x, 4, 2, 5), [NaN NaN NaN NaN NaN])
%!assert (hygecdf ([x(1:2) NaN x(4:5)], 4, 2, 2), [y(1:2) NaN y(4:5)], 5*eps)
%!test
%! p = hygecdf (x, 10, [1 2 3 4 5], 2, "upper");
%! assert (p, [1, 34/90, 2/30, 0, 0], 10*eps);
%!test
%! p = hygecdf (2*x, 10, [1 2 3 4 5], 2, "upper");
%! assert (p, [1, 34/90, 0, 0, 0], 10*eps);

## Test class of input preserved
%!assert (hygecdf ([x, NaN], 4, 2, 2), [y, NaN], 5*eps)
%!assert (hygecdf (single ([x, NaN]), 4, 2, 2), single ([y, NaN]), ...
%! eps ("single"))
%!assert (hygecdf ([x, NaN], single (4), 2, 2), single ([y, NaN]), ...
%! eps ("single"))
%!assert (hygecdf ([x, NaN], 4, single (2), 2), single ([y, NaN]), ...
%! eps ("single"))
%!assert (hygecdf ([x, NaN], 4, 2, single (2)), single ([y, NaN]), ...
%! eps ("single"))

## Test input validation
%!error<hygecdf: too few input arguments.> hygecdf ()
%!error<hygecdf: too few input arguments.> hygecdf (1)
%!error<hygecdf: too few input arguments.> hygecdf (1,2)
%!error<hygecdf: too few input arguments.> hygecdf (1,2,3)
%!error<hygecdf: invalid argument for upper tail.> hygecdf (1,2,3,4,5)
%!error<hygecdf: invalid argument for upper tail.> hygecdf (1,2,3,4,"uper")
%!error<hygecdf: X, T, M, and N must be of common size or scalars.> ...
%! hygecdf (ones (2), ones (3), 1, 1)
%!error<hygecdf: X, T, M, and N must be of common size or scalars.> ...
%! hygecdf (1, ones (2), ones (3), 1)
%!error<hygecdf: X, T, M, and N must be of common size or scalars.> ...
%! hygecdf (1, 1, ones (2), ones (3))
%!error<hygecdf: X, T, M, and N must not be complex.> hygecdf (i, 2, 2, 2)
%!error<hygecdf: X, T, M, and N must not be complex.> hygecdf (2, i, 2, 2)
%!error<hygecdf: X, T, M, and N must not be complex.> hygecdf (2, 2, i, 2)
%!error<hygecdf: X, T, M, and N must not be complex.> hygecdf (2, 2, 2, i)
