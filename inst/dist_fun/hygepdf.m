## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1996-2016 Kurt Hornik
## Copyright (C) 2022 Nicholas R. Jankowski
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} hygepdf (@var{x}, @var{m}, @var{k}, @var{n})
## @deftypefnx {statistics} {@var{y} =} hygepdf (@dots{}, @qcode{"vectorexpand"})
##
## Hypergeometric probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the hypergeometric distribution with parameters @var{m}, @var{k}, and
## @var{n}.  The size of @var{y} is the common size of @var{x}, @var{m},
## @var{k}, and @var{n}.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## This is the probability of obtaining @var{x} marked items when randomly
## drawing a sample of size @var{n} without replacement from a population of
## total size @var{m} containing @var{k} marked items.  The parameters @var{m},
## @var{k}, and @var{n} must be positive integers with @var{k} and @var{n} not
## greater than @var{m}.
##
## If the optional parameter @qcode{vectorexpand} is provided, @var{x} may be an
## array with size different from parameters @var{m}, @var{k}, and @var{n}
## (which must still be of a common size or scalar).  Each element of @var{x}
## will be evaluated against each set of parameters @var{m}, @var{k}, and
## @var{n} in columnwise order. The output @var{y} will be an array of size
## @qcode{@var{r} x @var{s}}, where @qcode{@var{r} = numel (@var{m})}, and
## @qcode{@var{s} = numel (@var{x})}.
##
## Further information about the hypergeometric distribution can be found at
## @url{https://en.wikipedia.org/wiki/Hypergeometric_distribution}
##
## @seealso{hygecdf, hygeinv, hygernd, hygestat}
## @end deftypefn

function y = hygepdf (x, m, k, n, vect_expand)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("hygepdf: function called with too few input arguments.");
  endif

  ## Check for X, T, M, and N being reals
  if (iscomplex (x) || iscomplex (m) || iscomplex (k) || iscomplex (n))
    error ("hygepdf: X, T, M, and N must not be complex.");
  endif

  ## Check for 5th argument or add default
  if (nargin < 5)
    vect_expand = [];
  endif

  if strcmpi (vect_expand, "vectorexpand")
    ## Expansion to improve vectorization of hyge calling functions.
    ## Project inputs over a 2D array with x(:) as a row vector and m,k,n as
    ## a column vector. each y(i,j) is hygepdf(x(j), m(i), k(i), n(i))
    ## Following expansion, remainder of algorithm processes as normal.

    if (! isscalar (m) || ! isscalar (k) || ! isscalar (n))
      [retval, m, k, n] = common_size (m, k, n);
      if (retval > 0)
        error ("hygepdf: T, M, and N must be of common size or scalars.");
      endif
      ## Ensure col vectors before expansion
      m = m(:);
      k = k(:);
      n = n(:);
    endif

    ## Expand x,m,k,n to arrays of size numel(m) x numel(x)
    sz = [numel(m), numel(x)];
    x = x(:)'; # ensure row vector before expansion

    x = x(ones (sz(1), 1), :);
    m = m(:, ones (sz(2), 1));
    k = k(:, ones (sz(2), 1));
    n = n(:, ones (sz(2), 1));

  else

    ## Check for common size of X, T, M, and N
    if (! isscalar (m) || ! isscalar (k) || ! isscalar (n))
      [retval, x, m, k, n] = common_size (x, m, k, n);
      if (retval > 0)
        error ("hygepdf: X, T, M, and N must be of common size or scalars.");
      endif
    endif

    sz = size (x);

  endif

  ## Check for class type
  if (isa (x, "single") || isa (m, "single")
                        || isa (k, "single") || isa (n, "single"))
    y = zeros (sz, "single");
  else
    y = zeros (sz);
  endif

  ## Everything in nel gives NaN
  nel = (isnan (x) | (m < 0) | (k < 0) | (n <= 0) | (k > m) | (n > m) |
        (m != fix (m)) | (k != fix (k)) | (n != fix (n)));
  ## Everything in zel gives 0 unless in nel
  zel = ((x != fix (x)) | (x < 0) | (x > k) | (n < x) | (n-x > m-k));

  y(nel) = NaN;

  ok = ! nel & ! zel;
  if (any (ok(:)))
    if (isscalar (m))
      y(ok) = exp (gammaln (k+1) - gammaln (k-x(ok)+1) - gammaln (x(ok)+1) + ...
                   gammaln (m-k+1) - gammaln (m-k-n+x(ok)+1) - ...
                   gammaln (n-x(ok)+1) - gammaln (m+1) + gammaln (m-n+1) + ...
                   gammaln (n+1));
    else
      y(ok) = exp (gammaln (k(ok)+1) - gammaln (k(ok)-x(ok)+1) - ...
                   gammaln (x(ok)+1) + gammaln (m(ok)-k(ok)+1) - ...
                   gammaln (m(ok)-k(ok)-n(ok)+x(ok)+1) - ...
                   gammaln (n(ok)-x(ok)+1) - gammaln (m(ok)+1) + ...
                   gammaln (m(ok)-n(ok)+1) + gammaln (n(ok)+1));
    endif
  endif

endfunction

%!demo
%! ## Plot various PDFs from the hypergeometric distribution
%! x = 0:60;
%! y1 = hygepdf (x, 500, 50, 100);
%! y2 = hygepdf (x, 500, 60, 200);
%! y3 = hygepdf (x, 500, 70, 300);
%! plot (x, y1, "*b", x, y2, "*g", x, y3, "*r")
%! grid on
%! xlim ([0, 60])
%! ylim ([0, 0.18])
%! legend ({"m = 500, k = 50, μ = 100", "m = 500, k = 60, μ = 200", ...
%!          "m = 500, k = 70, μ = 300"}, "location", "northeast")
%! title ("Hypergeometric PDF")
%! xlabel ("values in x (number of successes)")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 1 2 3];
%! y = [0 1/6 4/6 1/6 0];
%!assert (hygepdf (x, 4 * ones (1, 5), 2, 2), y, 3 * eps)
%!assert (hygepdf (x, 4, 2 * ones (1, 5), 2), y, 3 * eps)
%!assert (hygepdf (x, 4, 2, 2 * ones (1, 5)), y, 3 * eps)
%!assert (hygepdf (x, 4 * [1, -1, NaN, 1.1, 1], 2, 2), [0, NaN, NaN, NaN, 0])
%!assert (hygepdf (x, 4, 2 * [1, -1, NaN, 1.1, 1], 2), [0, NaN, NaN, NaN, 0])
%!assert (hygepdf (x, 4, 5, 2), [NaN, NaN, NaN, NaN, NaN], 3 * eps)
%!assert (hygepdf (x, 4, 2, 2 * [1, -1, NaN, 1.1, 1]), [0, NaN, NaN, NaN, 0])
%!assert (hygepdf (x, 4, 2, 5), [NaN, NaN, NaN, NaN, NaN], 3 * eps)
%!assert (hygepdf ([x, NaN], 4, 2, 2), [y, NaN], 3 * eps)

## Test class of input preserved
%!assert (hygepdf (single ([x, NaN]), 4, 2, 2), single ([y, NaN]), eps ("single"))
%!assert (hygepdf ([x, NaN], single (4), 2, 2), single ([y, NaN]), eps ("single"))
%!assert (hygepdf ([x, NaN], 4, single (2), 2), single ([y, NaN]), eps ("single"))
%!assert (hygepdf ([x, NaN], 4, 2, single (2)), single ([y, NaN]), eps ("single"))

## Test vector expansion
%!test
%! z = zeros(3,5);
%! z([4,5,6,8,9,12]) = [1, 0.5, 1/6, 0.5, 2/3, 1/6];
%! assert (hygepdf (x, 4, [0, 1, 2], 2, "vectorexpand"), z, 3 * eps);
%! assert (hygepdf (x, 4, [0, 1, 2]', 2, "vectorexpand"), z, 3 * eps);
%! assert (hygepdf (x', 4, [0, 1, 2], 2, "vectorexpand"), z, 3 * eps);
%! assert (hygepdf (2, 4, [0 ,1, 2], 2, "vectorexpand"), z(:,4), 3 * eps);
%! assert (hygepdf (x, 4, 1, 2, "vectorexpand"), z(2,:), 3 *eps);
%! assert (hygepdf ([NaN, x], 4, [0 1 2]', 2, "vectorexpand"), [NaN(3, 1), z], 3 * eps);

## Test input validation
%!error<hygepdf: function called with too few input arguments.> hygepdf ()
%!error<hygepdf: function called with too few input arguments.> hygepdf (1)
%!error<hygepdf: function called with too few input arguments.> hygepdf (1,2)
%!error<hygepdf: function called with too few input arguments.> hygepdf (1,2,3)
%!error<hygepdf: X, T, M, and N must be of common size or scalars.> ...
%! hygepdf (1, ones (3), ones (2), ones (2))
%!error<hygepdf: X, T, M, and N must be of common size or scalars.> ...
%! hygepdf (1, ones (2), ones (3), ones (2))
%!error<hygepdf: X, T, M, and N must be of common size or scalars.> ...
%! hygepdf (1, ones (2), ones (2), ones (3))
%!error<hygepdf: X, T, M, and N must not be complex.> hygepdf (i, 2, 2, 2)
%!error<hygepdf: X, T, M, and N must not be complex.> hygepdf (2, i, 2, 2)
%!error<hygepdf: X, T, M, and N must not be complex.> hygepdf (2, 2, i, 2)
%!error<hygepdf: X, T, M, and N must not be complex.> hygepdf (2, 2, 2, i)
