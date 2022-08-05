## Copyright (C) 2022 Nicholas R. Jankowski
## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1996-2016 Kurt Hornik
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
## @deftypefn {} {@var{pdf} =} hygepdf (@var{x}, @var{t}, @var{m}, @var{n})
## @deftypefnx {} {@var{pdf} =} hygepdf (@dots{}, "vectorexpand")
## Compute the probability density function (PDF) at @var{x} of the
## hypergeometric distribution with parameters @var{t}, @var{m}, and @var{n}.
##
## This is the probability of obtaining @var{x} marked items when randomly
## drawing a sample of size @var{n} without replacement from a population of
## total size @var{t} containing @var{m} marked items.
##
## The parameters @var{t}, @var{m}, and @var{n} must be positive integers
## with @var{m} and @var{n} not greater than @var{t}.  They and @var{x} may be
## scalars or arrays, but all non-scalars must be the same size.
##
## The output @var{pdf} will be the same size as the input array.
##
## If the optional parameter @code{vectorexpand} is provided, @var{x} may be an
## array with size different from parameters @var{t}, @var{m}, and @var{n}
## (which must still be of a common size or scalar).  Each element of @var{x}
## will be evaluated against each set of parameters @var{t}, @var{m}, and
## @var{n} in columnwise order. The output @var{pdf} will be an array of size
## @code{@var{r} x @var{s}}, where @code{@var{r} = numel (@var{t})}, and
## @code{@var{s} = numel (@var{x})}.
##
## @end deftypefn

## Author: KH <Kurt.Hornik@wu-wien.ac.at>
## Description: PDF of the hypergeometric distribution

function pdf = hygepdf (x, t, m, n, vect_expand = "")

  if (!any (nargin == [4,5]))
    print_usage ();
  endif

  if (iscomplex (x) || iscomplex (t) || iscomplex (m) || iscomplex (n))
    error ("hygepdf: X, T, M, and N must not be complex");
  endif

  if strcmpi (vect_expand, "vectorexpand")
    ## Expansion to improve vectorization of hyge calling functions.
    ## Project inputs over a 2D array with x(:) as a row vector and t,m,n as
    ## a column vector. each pdf(i,j) is hygepdf(x(j), t(i), m(i), n(i))
    ## Following expansion, remainder of algorithm processes as normal.

    if (! isscalar (t) || ! isscalar (m) || ! isscalar (n))
      [retval, t, m, n] = common_size (t, m, n);
      if (retval > 0)
        error ("hygepdf: T, M, and N must be of common size or scalars");
      endif
      ## ensure col vectors before expansion
      t = t(:);
      m = m(:);
      n = n(:);
    endif

    ## expand x,t,m,n to arrays of size numel(t) x numel(x)
    sz = [numel(t), numel(x)];
    x = x(:)'; # ensure row vector before expansion

    x = x(ones (sz(1), 1), :);
    t = t(:, ones (sz(2), 1));
    m = m(:, ones (sz(2), 1));
    n = n(:, ones (sz(2), 1));

  else

    if (! isscalar (t) || ! isscalar (m) || ! isscalar (n))
      [retval, x, t, m, n] = common_size (x, t, m, n);
      if (retval > 0)
        error ("hygepdf: X, T, M, and N must be of common size or scalars");
      endif
    endif

    sz = size (x);

  endif

  if (isa (x, "single") || isa (t, "single")
      || isa (m, "single") || isa (n, "single"))
    pdf = zeros (sz, "single");
  else
    pdf = zeros (sz);
  endif

  ## everything in nel gives NaN
  nel = (isnan (x) | (t < 0) | (m < 0) | (n <= 0) | (m > t) | (n > t) |
        (t != fix (t)) | (m != fix (m)) | (n != fix (n)));
  ## everything in zel gives 0 unless in nel
  zel = ((x != fix (x)) | (x < 0) | (x > m) | (n < x) | (n-x > t-m));

  pdf(nel) = NaN;

  k = ! nel & ! zel;
  if (any (k(:)))
    if (isscalar (t))
      pdf(k) = exp (gammaln (m+1) - gammaln (m-x(k)+1) - gammaln (x(k)+1) +...
                 gammaln (t-m+1) - gammaln (t-m-n+x(k)+1) - ...
                 gammaln (n-x(k)+1) - gammaln (t+1) + gammaln (t-n+1) + ...
                 gammaln (n+1));
    else
      pdf(k) = exp (gammaln (m(k)+1) - gammaln (m(k)-x(k)+1) - ...
                 gammaln (x(k)+1) + gammaln (t(k)-m(k)+1) - ...
                 gammaln (t(k)-m(k)-n(k)+x(k)+1) - gammaln (n(k)-x(k)+1) - ...
                 gammaln (t(k)+1) + gammaln (t(k)-n(k)+1) + gammaln (n(k)+1));
    endif
  endif


endfunction


%!shared x,y
%! x = [-1 0 1 2 3];
%! y = [0 1/6 4/6 1/6 0];
%!assert (hygepdf (x, 4*ones (1,5), 2, 2), y, eps)
%!assert (hygepdf (x, 4, 2*ones (1,5), 2), y, eps)
%!assert (hygepdf (x, 4, 2, 2*ones (1,5)), y, eps)
%!assert (hygepdf (x, 4*[1 -1 NaN 1.1 1], 2, 2), [0 NaN NaN NaN 0], eps)
%!assert (hygepdf (x, 4, 2*[1 -1 NaN 1.1 1], 2), [0 NaN NaN NaN 0], eps)
%!assert (hygepdf (x, 4, 5, 2), [NaN NaN NaN NaN NaN], eps)
%!assert (hygepdf (x, 4, 2, 2*[1 -1 NaN 1.1 1]), [0 NaN NaN NaN 0], eps)
%!assert (hygepdf (x, 4, 2, 5), [NaN NaN NaN NaN NaN], eps)
%!assert (hygepdf ([x, NaN], 4, 2, 2), [y, NaN], eps)

## Test class of input preserved
%!assert (hygepdf (single ([x, NaN]), 4, 2, 2), single ([y, NaN]), eps("single"))
%!assert (hygepdf ([x, NaN], single (4), 2, 2), single ([y, NaN]), eps("single"))
%!assert (hygepdf ([x, NaN], 4, single (2), 2), single ([y, NaN]), eps("single"))
%!assert (hygepdf ([x, NaN], 4, 2, single (2)), single ([y, NaN]), eps("single"))

## Test vector expansion
%!test
%! z = zeros(3,5);
%! z([4,5,6,8,9,12]) = [1, 0.5, 1/6, 0.5, 2/3, 1/6];
%! assert (hygepdf (x, 4, [0 1 2], 2,"vectorexpand"), z, eps);
%! assert (hygepdf (x, 4, [0 1 2]', 2,"vectorexpand"), z, eps);
%! assert (hygepdf (x', 4, [0 1 2], 2,"vectorexpand"), z, eps);
%! assert (hygepdf (2, 4, [0 1 2], 2,"vectorexpand"), z(:,4), eps);
%! assert (hygepdf (x, 4, 1, 2,"vectorexpand"), z(2,:), eps);
%! assert (hygepdf ([NaN,x], 4, [0 1 2]', 2,"vectorexpand"),[NaN(3,1), z], eps);

## Test input validation
%!error hygepdf ()
%!error hygepdf (1)
%!error hygepdf (1,2)
%!error hygepdf (1,2,3)
%!error hygepdf (1,2,3,4,5,6)
%!error hygepdf (1, ones (3), ones (2), ones (2))
%!error hygepdf (1, ones (2), ones (3), ones (2))
%!error hygepdf (1, ones (2), ones (2), ones (3))
%!error hygepdf (i, 2, 2, 2)
%!error hygepdf (2, i, 2, 2)
%!error hygepdf (2, 2, i, 2)
%!error hygepdf (2, 2, 2, i)
