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
## @deftypefn  {statistics} {@var{m} =} nanmedian (@var{x})
## @deftypefnx {statistics} {@var{m} =} nanmedian (@var{x}, @qcode{'all'})
## @deftypefnx {statistics} {@var{m} =} nanmedian (@var{x}, @var{dim})
## @deftypefnx {statistics} {@var{m} =} nanmedian (@var{x}, @var{vecdim})
##
## Compute the median while ignoring NaN values.
##
## @code{@var{m} = nanmedian (@var{x})} returns the median of @var{x}, after
## removing @qcode{NaN} values.  If @var{x} is a vector, a scalar value is
## returned.  If @var{x} is a matrix, a row vector of column medians is
## returned.  If @var{x} is a multidimensional array, @code{nanmedian} operates
## along the first nonsingleton dimension.  If all values along a dimension are
## @qcode{NaN}, the median is returned as @qcode{NaN}.
##
## @code{@var{m} = nanmedian (@var{x}, @qcode{'all'})} returns the median of all
## elements of @var{x}, after removing @qcode{NaN} values.  It is the equivalent
## of @code{nanmedian (@var{x}(:))}.
##
## @code{@var{m} = nanmedian (@var{x}, @var{dim})} operates along the dimension
## @var{dim} of @var{x}.
##
## @code{@var{m} = nanmedian (@var{x}, @var{vecdim})} returns the median over
## the dimensions specified in the vector @var{vecdim}.  Each element of
## @var{vecdim} represents a dimension of the input array @var{x} and the output
## @var{m} has length 1 in the specified operating dimensions.  If @var{vecdim}
## indexes all dimensions of @var{x}, then it is equivalent to
## @code{nanmedian (@var{x}, @qcode{'all'})}.  Any dimension in @var{vecdim}
## greater than @code{ndims (@var{x})} is ignored.
##
## @seealso{median, nanmean, nansum}
## @end deftypefn

function m = nanmedian (x, dim)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif
  if (! isnumeric (x) && ! islogical (x))
    error ("nanmedian: X must be numeric.");
  endif
  if (isempty (x))
    m = NaN;
    return;
  endif

  ## Operating dimension(s)
  if (nargin < 2)
    dim = find (size (x) != 1, 1);
    if (isempty (dim))
      dim = 1;
    endif
  elseif (ischar (dim) && strcmpi (dim, 'all'))
    m = __colmedian__ (x(:).');
    return;
  elseif (isscalar (dim))
    if (! isnumeric (dim) || fix (dim) != dim || dim <= 0)
      error ("nanmedian: DIM must be a positive integer.");
    endif
  elseif (isnumeric (dim) && isvector (dim))
    if (any (fix (dim) != dim) || any (dim <= 0))
      error ("nanmedian: VECDIM must be a vector of positive integers.");
    endif
    dim = sort (dim);
    if (! all (diff (dim)))
      error ("nanmedian: VECDIM must contain non-repeating positive integers.");
    endif
  else
    error ("nanmedian: invalid DIM argument.");
  endif

  ## Ignore dimensions larger than the array (they are singleton)
  vecdim = dim;
  vecdim(vecdim > ndims (x)) = [];
  if (isempty (vecdim))
    m = x;
    return;
  endif

  szx = size (x);
  remdims = 1:ndims (x);
  remdims(vecdim) = [];
  if (isempty (remdims))
    m = __colmedian__ (x(:).');
    return;
  endif

  ## Permute the operating dimensions to the back and squash them together
  perm = [remdims, vecdim];
  xp = permute (x, perm);
  R = prod (szx(remdims));
  K = prod (szx(vecdim));
  xr = reshape (xp, R, K);
  mr = __colmedian__ (xr);

  ## Restore the original layout with the operating dimensions collapsed
  m = reshape (mr, [szx(remdims), ones(1, numel (vecdim))]);
  m = ipermute (m, perm);

endfunction

## Median over the second dimension of a 2-D matrix, ignoring NaN values.
function m = __colmedian__ (M)

  R = rows (M);
  Ms = sort (M, 2);              # NaN values are sorted to the end
  n = sum (! isnan (M), 2);
  m = NaN (R, 1);
  v = find (n > 0);
  nn = n(v);
  iL = floor ((nn + 1) / 2);
  iU = ceil ((nn + 1) / 2);
  linL = v + (iL - 1) * R;
  linU = v + (iU - 1) * R;
  m(v) = (Ms(linL) + Ms(linU)) / 2;

endfunction

%!demo
%! ## Find the column medians for a matrix with missing values.
%!
%! x = magic (3);
%! x([1, 6:9]) = NaN
%! m = nanmedian (x)

%!demo
%! ## Find the median of all elements, ignoring missing values.
%!
%! x = reshape (1:30, [2, 5, 3]);
%! x([10:12, 25]) = NaN
%! m = nanmedian (x, 'all')

## Test output
%!assert_equal (nanmedian ([]), NaN)
%!assert_equal (nanmedian (NaN), NaN)
%!assert_equal (nanmedian (5), 5)
%!assert_equal (nanmedian ([2, 4, NaN, 8]), 4)
%!assert_equal (nanmedian ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12]), [5.5, 8, 9])
%!assert_equal (nanmedian ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12], 2), ...
%!              [1.5; 5; 8; 11])
%!assert_equal (nanmedian ([1 NaN; NaN NaN; 3 NaN]), [2, NaN])
%!assert_equal (nanmedian (reshape (1:12, [2, 3, 2])(:)), 6.5)
%!test
%! x = reshape (1:24, [2, 4, 3]);
%! x([5:6, 20]) = NaN;
%! assert_equal (nanmedian (x, 'all'), nanmedian (x(! isnan (x))(:)));
%! assert_equal (nanmedian (x, [1, 2, 3]), nanmedian (x, 'all'));
%!test
%! x = magic (4);
%! x([1, 6, 11, 16]) = NaN;
%! assert_equal (nanmedian (x, 2), [3; 8; 9; 14]);

## Test input validation
%!error <Invalid call to nanmedian> nanmedian ()
%!error <nanmedian: X must be numeric.> nanmedian ({3})
%!error <nanmedian: DIM must be a positive integer.> nanmedian (ones (3), 0)
%!error <nanmedian: DIM must be a positive integer.> nanmedian (ones (3), 1.5)
%!error <nanmedian: VECDIM must contain non-repeating positive integers.> ...
%! nanmedian (ones (3, 3, 3), [2, 2, 3])
