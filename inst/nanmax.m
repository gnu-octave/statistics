## Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{v} =} nanmax (@var{x})
## @deftypefnx {statistics} {@var{v} =} nanmax (@var{x}, [], @var{dim})
## @deftypefnx {statistics} {[@var{v}, @var{idx}] =} nanmax (@dots{})
## @deftypefnx {statistics} {@var{v} =} nanmax (@var{x}, [], @qcode{'all'})
## @deftypefnx {statistics} {@var{v} =} nanmax (@var{x}, [], @var{vecdim})
## @deftypefnx {statistics} {@var{v} =} nanmax (@var{x}, @var{y})
##
## Find the maximum while ignoring NaN values.
##
## @code{@var{v} = nanmax (@var{x})} returns the maximum of @var{x}, after
## removing @qcode{NaN} values.  If @var{x} is a vector, a scalar maximum value
## is returned.  If @var{x} is a matrix, a row vector of column maxima is
## returned.  If @var{x} is a multidimentional array, the @code{nanmax} operates
## along the first nonsigleton dimension.  If all values in a column are
## @qcode{NaN}, the maximum is returned as @qcode{NaN} rather than @qcode{[]}.
##
## @code{@var{v} = nanmax (@var{x}, [], @var{dim})} operates along the dimension
## @var{dim} of @var{x}.
##
## @code{[@var{v}, @var{idx}] = nanmax (@dots{})} also returns the row indices
## of the maximum values for each column in the vector @var{idx}.  When @var{x}
## is a vector, then @var{idx} is a scalar value as @var{v}.
##
## @code{@var{v} = nanmax (@var{x}, [], @qcode{'all'})} returns the maximum of
## all elements of @var{x}, after removing @qcode{NaN} values.  It is the
## equivalent of @code{nanmax (@var{x}(:))}.   The optional flag @qcode{'all'}
## cannot be used together with @var{dim} or @var{vecdim} input arguments.
##
## @code{@var{v} = nanmax (@var{x}, [], @var{vecdim})} returns the maximum over
## the dimensions specified in the vector @var{vecdim}.  Each element of
## @var{vecdim} represents a dimension of the input array @var{x} and the output
## @var{v} has length 1 in the specified operating dimensions.  The lengths of
## the other dimensions are the same for @var{x} and @var{y}.  For example, if
## @var{x} is a 2-by-3-by-4 array, then @code{nanmax (@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the maximum of the
## elements on the corresponding page of @var{x}.  If @var{vecdim} indexes all
## dimensions of @var{x}, then it is equivalent to
## @code{nanmax (@var{x}, @qcode{'all'})}.  Any dimension in @var{vecdim}
## greater than @code{ndims (@var{x})} is ignored.
##
## @seealso{max, nanmin, nansum}
## @end deftypefn

function [v, idx] = nanmax (x, y, dim)
  if (nargin < 1 || nargin > 3)
    print_usage;
  elseif (nargin == 1 || (nargin == 2 && isempty (y)))
    nanvals = isnan (x);
    x(nanvals) = -Inf;
    [v, idx] = max (x);
    v(all (nanvals)) = NaN;
  elseif (nargin == 3 && strcmpi (dim, "all") && isempty (y))
    x = x(:);
    nanvals = isnan (x);
    x(nanvals) = -Inf;
    [v, idx] = max (x);
    v(all (nanvals)) = NaN;
  elseif (nargin == 3 && isempty (y))
    if (isscalar (dim))
      nanvals = isnan (x);
      x(nanvals) = -Inf;
      [v, idx] = max (x, [], dim);
      v(all (nanvals, dim)) = NaN;
    else
      vecdim = sort (dim);
      if (! all (diff (vecdim)))
         error ("nanmax: VECDIM must contain non-repeating positive integers.");
      endif
      ## Ignore dimensions in VECDIM larger than actual array
      vecdim(find (vecdim > ndims (x))) = [];

      if (isempty (vecdim))
        v = x;
        if (nargout > 1)
          idx = reshape ([1:numel(x)], size (x));
        endif
      else

        ## Calculate permutation vector
        szx = size (x);
        remdims = 1:ndims (x);      # All dimensions
        remdims(vecdim) = [];       # Delete dimensions specified by vecdim
        nremd = numel (remdims);

        ## If all dimensions are given, it is equivalent to 'all' flag
        if (nremd == 0)
          x = x(:);
          nanvals = isnan (x);
          x(nanvals) = -Inf;
          [v, idx] = max (x);
          v(all (nanvals)) = NaN;

        else
          ## Permute to push vecdims to back
          perm = [remdims, vecdim];
          x = permute (x, perm);

          ## Reshape to squash all vecdims in final dimension
          sznew = [szx(remdims), prod(szx(vecdim))];
          x = reshape (x, sznew);

          ## Calculate nanmax on final dimension
          dim = nremd + 1;
          nanvals = isnan (x);
          x(nanvals) = -Inf;
          [v, idx] = max (x, [], dim);
          v(all (nanvals, dim)) = NaN;

          ## Inverse permute back to correct dimensions
          v = ipermute (v, perm);
          idx = ipermute (idx, perm);
        endif
      endif
    endif
  else
    if (nargout > 1)
      error ("nanmax: a second output is not supported with this syntax.");
    endif
    Xnan = isnan (x);
    Ynan = isnan (y);
    x(Xnan) = -Inf;
    y(Ynan) = -Inf;
    v = max (x, y);
    v(Xnan & Ynan) = NaN;
  endif
endfunction

%!demo
%! ## Find the column maximum values and their indices
%! ## for matrix data with missing values.
%!
%! x = magic (3);
%! x([1, 6:9]) = NaN
%! [y, ind] = nanmax (x)

%!demo
%! ## Find the maximum of all the values in an array, ignoring missing values.
%! ## Create a 2-by-5-by-3 array x with some missing values.
%!
%! x = reshape (1:30, [2, 5, 3]);
%! x([10:12, 25]) = NaN
%!
%! ## Find the maximum of the elements of x.
%!
%! y = nanmax (x, [], 'all')

## Test output
%!assert (nanmax ([2, 4, NaN, 7]), 7)
%!assert (nanmax ([2, 4, NaN, Inf]), Inf)
%!assert (nanmax ([1, NaN, 3; NaN, 5, 6; 7, 8, NaN]), [7, 8, 6])
%!assert (nanmax ([1, NaN, 3; NaN, 5, 6; 7, 8, NaN]'), [3, 6, 8])
%!assert (nanmax (single ([1, NaN, 3; NaN, 5, 6; 7, 8, NaN])), single ([7, 8, 6]))
%!shared x, y
%! x(:,:,1) = [1.77, -0.005, NaN, -2.95; NaN, 0.34, NaN, 0.19];
%! x(:,:,2) = [1.77, -0.005, NaN, -2.95; NaN, 0.34, NaN, 0.19] + 5;
%! y = x;
%! y(2,3,1) = 0.51;
%!assert (nanmax (x, [], [1, 2])(:), [1.77;6.77])
%!assert (nanmax (x, [], [1, 3])(:), [6.77;5.34;NaN;5.19])
%!assert (nanmax (x, [], [2, 3])(:), [6.77;5.34])
%!assert (nanmax (x, [], [1, 2, 3]), 6.77)
%!assert (nanmax (x, [], 'all'), 6.77)
%!assert (nanmax (y, [], [1, 3])(:), [6.77;5.34;0.51;5.19])
%!assert (nanmax (x(1,:,1), x(2,:,1)), [1.77, 0.34, NaN, 0.19])
%!assert (nanmax (x(1,:,2), x(2,:,2)), [6.77, 5.34, NaN, 5.19])
%!assert (nanmax (y(1,:,1), y(2,:,1)), [1.77, 0.34, 0.51, 0.19])
%!assert (nanmax (y(1,:,2), y(2,:,2)), [6.77, 5.34, NaN, 5.19])

## Test dimension indexing with vecdim in N-dimensional arrays
%!test
%! xx = repmat ([1:20;6:25], [5 2 6 3]);
%! assert (size (nanmax (xx, [], [3, 2])), [10, 1, 1, 3]);
%! assert (size (nanmax (xx, [], [1, 2])), [1, 1, 6, 3]);
%! assert (size (nanmax (xx, [], [1, 2, 4])), [1, 1, 6]);
%! assert (size (nanmax (xx, [], [1, 4, 3])), [1, 40]);
%! assert (size (nanmax (xx, [], [1, 2, 3, 4])), [1, 1]);

## Test exceeding dimensions
%!assert (nanmax (ones (2), [], 3), ones (2, 2))
%!assert (nanmax (ones (2, 2, 2), [], 99), ones (2, 2, 2))
%!assert (nanmax (magic (3), [], 3), magic (3))
%!assert (nanmax (magic (3), [], [1, 3]), [8, 9, 7])
%!assert (nanmax (magic (3), [], [1, 99]), [8, 9, 7])

## Test comparisons
%!assert (nanmax (ones (2), 3), 3 * ones (2,2))

## Test input validation
%!error <nanmax: VECDIM must contain non-repeating positive integers.> ...
%! nanmax (y, [], [1, 1, 2])
%!error <nanmax: a second output is not supported with this syntax.> ...
%! [v, idx] = nanmax(x, y, [1 2])
