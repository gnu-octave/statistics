## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{v} =} nanvar (@var{x})
## @deftypefnx {statistics} {@var{v} =} nanvar (@var{x}, @var{w})
## @deftypefnx {statistics} {@var{v} =} nanvar (@var{x}, @var{w}, @qcode{'all'})
## @deftypefnx {statistics} {@var{v} =} nanvar (@var{x}, @var{w}, @var{dim})
## @deftypefnx {statistics} {@var{v} =} nanvar (@var{x}, @var{w}, @var{vecdim})
##
## Compute the variance while ignoring NaN values.
##
## @code{@var{v} = nanvar (@var{x})} returns the variance of @var{x}, after
## removing @qcode{NaN} values.  If @var{x} is a vector, a scalar value is
## returned.  If @var{x} is a matrix, a row vector of column variances is
## returned.  If @var{x} is a multidimensional array, @code{nanvar} operates
## along the first nonsingleton dimension.  If a dimension contains fewer than
## two non-@qcode{NaN} values, the variance is returned as 0 for a single value
## and as @qcode{NaN} when all values are @qcode{NaN}.
##
## @code{@var{v} = nanvar (@var{x}, @var{w})} specifies the normalization.  When
## @var{w} is 0 (default), the variance is normalized by @math{N-1}, where
## @math{N} is the number of non-@qcode{NaN} observations.  When @var{w} is 1,
## it is normalized by @math{N}.  @var{w} may also be a vector of nonnegative
## weights whose length matches the operating dimension, in which case the
## weighted variance normalized by the sum of the weights is returned.
##
## @code{@var{v} = nanvar (@var{x}, @var{w}, @qcode{'all'})} returns the
## variance of all elements of @var{x}, after removing @qcode{NaN} values.  Use
## an empty value, @code{@var{w} = []}, to pass the default normalization.
##
## @code{@var{v} = nanvar (@var{x}, @var{w}, @var{dim})} operates along the
## dimension @var{dim} of @var{x}.
##
## @code{@var{v} = nanvar (@var{x}, @var{w}, @var{vecdim})} returns the variance
## over the dimensions specified in the vector @var{vecdim}.  A weight vector is
## not supported together with @qcode{'all'} or @var{vecdim}.  Any dimension in
## @var{vecdim} greater than @code{ndims (@var{x})} is ignored.
##
## @seealso{var, nanstd, nanmean, nansum}
## @end deftypefn

function v = nanvar (x, w, dim)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif
  if (! isnumeric (x) && ! islogical (x))
    error ("nanvar: X must be numeric.");
  endif
  if (isempty (x))
    v = NaN;
    return;
  endif

  ## Normalization / weight argument
  if (nargin < 2 || isempty (w))
    w = 0;
  endif
  if (! isnumeric (w) || ! isreal (w))
    error ("nanvar: W must be 0, 1, or a vector of nonnegative weights.");
  endif
  wvec = false;
  if (isscalar (w))
    if (w != 0 && w != 1)
      error ("nanvar: W must be 0 or 1 when it is a scalar.");
    endif
  elseif (isvector (w))
    if (any (w < 0) || any (isnan (w)))
      error ("nanvar: weight vector W must contain nonnegative values.");
    endif
    wvec = true;
  else
    error ("nanvar: W must be 0, 1, or a vector of nonnegative weights.");
  endif

  ## Operating dimension
  dimall = false;
  if (nargin < 3)
    dim = find (size (x) != 1, 1);
    if (isempty (dim))
      dim = 1;
    endif
  elseif (ischar (dim) && strcmpi (dim, 'all'))
    dimall = true;
  elseif (isscalar (dim))
    if (! isnumeric (dim) || fix (dim) != dim || dim <= 0)
      error ("nanvar: DIM must be a positive integer.");
    endif
  elseif (isnumeric (dim) && isvector (dim))
    if (any (fix (dim) != dim) || any (dim <= 0))
      error ("nanvar: VECDIM must be a vector of positive integers.");
    endif
    dim = sort (dim);
    if (! all (diff (dim)))
      error ("nanvar: VECDIM must contain non-repeating positive integers.");
    endif
  else
    error ("nanvar: invalid DIM argument.");
  endif

  if (dimall)
    dimarg = 'all';
  else
    dimarg = dim;
  endif

  if (wvec)
    ## Weighted variance (biased, normalized by the sum of the weights)
    if (dimall || ! isscalar (dim))
      error ("nanvar: a weight vector is supported only with a scalar DIM.");
    endif
    if (numel (w) != size (x, dim))
      error ("nanvar: the length of W must match the operating dimension.");
    endif
    wsz = ones (1, max (ndims (x), dim));
    wsz(dim) = numel (w);
    wr = reshape (w(:), wsz);
    mask = ! isnan (x);
    we = wr .* mask;
    sumw = sum (we, dim);
    x0 = x;
    x0(! mask) = 0;
    mu = sum (we .* x0, dim) ./ sumw;
    t = (x0 - mu) .^ 2;
    v = sum (we .* t, dim) ./ sumw;
  else
    mask = ! isnan (x);
    n = sum (mask, dimarg);
    x0 = x;
    x0(! mask) = 0;
    xm = sum (x0, dimarg) ./ n;
    d2 = (x - xm) .^ 2;
    d2(isnan (d2)) = 0;
    ss = sum (d2, dimarg);
    if (w == 1)
      v = ss ./ n;
    else
      v = ss ./ (n - 1);
    endif
    v(n == 1) = 0;
    v(n == 0) = NaN;
  endif

endfunction

%!demo
%! ## Find the column variances for a matrix with missing values.
%!
%! x = magic (3);
%! x([1, 6:9]) = NaN
%! v = nanvar (x)

%!demo
%! ## Find the row variances, normalized by N instead of N-1.
%!
%! x = magic (3);
%! x([1, 6:9]) = NaN
%! v = nanvar (x, 1, 2)

## Test output
%!assert_equal (nanvar ([]), NaN)
%!assert_equal (nanvar (NaN), NaN)
%!assert_equal (nanvar (5), 0)
%!assert_equal (nanvar ([2, 4, NaN, 8]), 9.333333333333334, 1e-14)
%!assert (nanvar ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12]), [15, 21, 9])
%!assert (nanvar ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12], 1), [11.25, 14, 6])
%!assert (nanvar ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12], 0, 2), [0.5; 2; 1; 1])
%!assert (nanvar ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12], [1 2 3 4]'), ...
%!        [9, 8.4375, 50/9], 1e-13)
%!assert_equal (nanvar (NaN (2, 3)), [NaN, NaN, NaN])
%!test
%! x = reshape (1:24, [2, 4, 3]);
%! x([5:6, 20]) = NaN;
%! assert (nanvar (x, 0, 'all'), nanvar (x(! isnan (x))(:)), 1e-12)

## Test input validation
%!error <Invalid call to nanvar> nanvar ()
%!error <nanvar: X must be numeric.> nanvar ({3})
%!error <nanvar: W must be 0 or 1 when it is a scalar.> nanvar (ones (3), 2)
%!error <nanvar: weight vector W must contain nonnegative values.> ...
%! nanvar (ones (1, 3), [1, -1, 2])
%!error <nanvar: DIM must be a positive integer.> nanvar (ones (3), 0, 1.5)
%!error <nanvar: VECDIM must contain non-repeating positive integers.> ...
%! nanvar (ones (3, 3, 3), 0, [2, 2, 3])
%!error <nanvar: a weight vector is supported only with a scalar DIM.> ...
%! nanvar (ones (2, 3), [1, 2], 'all')
%!error <nanvar: the length of W must match the operating dimension.> ...
%! nanvar (ones (2, 3), [1, 2, 3], 1)
