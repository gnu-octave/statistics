## Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
## Copyright (C) 2003 Alois Schloegl <alois.schloegl@ist.ac.at>
## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{v} =} nanmin (@var{x})
## @deftypefnx {statistics} {@var{v} =} nanmin (@var{x}, [], @var{dim})
## @deftypefnx {statistics} {[@var{v}, @var{idx}] =} nanmin (@dots{})
## @deftypefnx {statistics} {@var{v} =} nanmin (@var{x}, [], @qcode{'all'})
## @deftypefnx {statistics} {@var{v} =} nanmin (@var{x}, [], @var{vecdim})
## @deftypefnx {statistics} {@var{v} =} nanmin (@var{x}, @var{y})
##
## Find the minimum while ignoring NaN values.
##
## @code{@var{v} = nanmin (@var{x})} returns the minimum of @var{x}, after
## removing @qcode{NaN} values.  If @var{x} is a vector, a scalar minimum value
## is returned.  If @var{x} is a matrix, a row vector of column maxima is
## returned.  If @var{x} is a multidimensional array, the @code{nanmin} operates
## along the first nonsigleton dimension.  If all values in a column are
## @qcode{NaN}, the minimum is returned as @qcode{NaN} rather than @qcode{[]}.
##
## @code{@var{v} = nanmin (@var{x}, [], @var{dim})} operates along the dimension
## @var{dim} of @var{x}.
##
## @code{[@var{v}, @var{idx}] = nanmin (@dots{})} also returns the row indices
## of the minimum values for each column in the vector @var{idx}.  When @var{x}
## is a vector, then @var{idx} is a scalar value as @var{v}.
##
## @code{@var{v} = nanmin (@var{x}, [], @qcode{'all'})} returns the minimum of
## all elements of @var{x}, after removing @qcode{NaN} values.  It is the
## equivalent of @code{nanmin (@var{x}(:))}.   The optional flag @qcode{'all'}
## cannot be used together with @var{dim} or @var{vecdim} input arguments.
##
## @code{@var{v} = nanmin (@var{x}, [], @var{vecdim})} returns the minimum over
## the dimensions specified in the vector @var{vecdim}.  Each element of
## @var{vecdim} represents a dimension of the input array @var{x} and the output
## @var{v} has length 1 in the specified operating dimensions.  The lengths of
## the other dimensions are the same for @var{x} and @var{y}.  For example, if
## @var{x} is a 2-by-3-by-4 array, then @code{nanmin (@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the minimum of the
## elements on the corresponding page of @var{x}.  If @var{vecdim} indexes all
## dimensions of @var{x}, then it is equivalent to
## @code{nanmin (@var{x}, @qcode{'all'})}.  Any dimension in @var{vecdim}
## greater than @code{ndims (@var{x})} is ignored.
##
## @seealso{min, nanmax, nansum}
## @end deftypefn

function [v, idx] = nanmin (x, y, dim)
  if (nargin < 1 || nargin > 3)
    print_usage;
  elseif (nargin == 1 || (nargin == 2 && isempty (y)))
    nanvals = isnan (x);
    x(nanvals) = Inf;
    [v, idx] = min (x);
    v(all (nanvals)) = NaN;
  elseif (nargin == 3 && strcmpi (dim, "all") && isempty (y))
    x = x(:);
    nanvals = isnan (x);
    x(nanvals) = Inf;
    [v, idx] = min (x);
    v(all (nanvals)) = NaN;
  elseif (nargin == 3 && isempty (y))
    if (isscalar (dim))
      nanvals = isnan (x);
      x(nanvals) = Inf;
      [v, idx] = min (x, [], dim);
      v(all (nanvals, dim)) = NaN;
    else
      vecdim = sort (dim);
      if (! all (diff (vecdim)))
         error ("nanmin: VECDIM must contain non-repeating positive integers.");
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
          x(nanvals) = Inf;
          [v, idx] = min (x);
          v(all (nanvals)) = NaN;

        else
          ## Permute to push vecdims to back
          perm = [remdims, vecdim];
          x = permute (x, perm);

          ## Reshape to squash all vecdims in final dimension
          sznew = [szx(remdims), prod(szx(vecdim))];
          x = reshape (x, sznew);

          ## Calculate nanmin on final dimension
          dim = nremd + 1;
          nanvals = isnan (x);
          x(nanvals) = Inf;
          [v, idx] = min (x, [], dim);
          v(all (nanvals, dim)) = NaN;

          ## Inverse permute back to correct dimensions
          v = ipermute (v, perm);
          idx = ipermute (idx, perm);
        endif
      endif
    endif
  else
    if (nargout > 1)
      error ("nanmin: a second output is not supported with this syntax.");
    endif
    Xnan = isnan (x);
    Ynan = isnan (y);
    x(Xnan) = Inf;
    y(Ynan) = Inf;
    v = min (x, y);
    v(Xnan & Ynan) = NaN;
  endif
endfunction

%!demo
%! ## Find the column minimum values and their indices
%! ## for matrix data with missing values.
%!
%! x = magic (3);
%! x([1, 6:9]) = NaN
%! [y, ind] = nanmin (x)

%!demo
%! ## Find the minimum of all the values in an array, ignoring missing values.
%! ## Create a 2-by-5-by-3 array x with some missing values.
%!
%! x = reshape (1:30, [2, 5, 3]);
%! x([10:12, 25]) = NaN
%!
%! ## Find the minimum of the elements of x.
%!
%! y = nanmin (x, [], 'all')

## Test output
%!assert (nanmin ([2, 4, NaN, 7]), 2)
%!assert (nanmin ([2, 4, NaN, -Inf]), -Inf)
%!assert (nanmin ([1, NaN, 3; NaN, 5, 6; 7, 8, NaN]), [1, 5, 3])
%!assert (nanmin ([1, NaN, 3; NaN, 5, 6; 7, 8, NaN]'), [1, 5, 7])
%!assert (nanmin (single ([1, NaN, 3; NaN, 5, 6; 7, 8, NaN])), single ([1, 5, 3]))
%!shared x, y
%! x(:,:,1) = [1.77, -0.005, NaN, -2.95; NaN, 0.34, NaN, 0.19];
%! x(:,:,2) = [1.77, -0.005, NaN, -2.95; NaN, 0.34, NaN, 0.19] + 5;
%! y = x;
%! y(2,3,1) = 0.51;
%!assert (nanmin (x, [], [1, 2])(:), [-2.95; 2.05])
%!assert (nanmin (x, [], [1, 3])(:), [1.77; -0.005; NaN; -2.95])
%!assert (nanmin (x, [], [2, 3])(:), [-2.95; 0.19])
%!assert (nanmin (x, [], [1, 2, 3]), -2.95)
%!assert (nanmin (x, [], 'all'), -2.95)
%!assert (nanmin (y, [], [1, 3])(:), [1.77; -0.005; 0.51; -2.95])
%!assert (nanmin (x(1,:,1), x(2,:,1)), [1.77, -0.005, NaN, -2.95])
%!assert (nanmin (x(1,:,2), x(2,:,2)), [6.77, 4.995, NaN, 2.05])
%!assert (nanmin (y(1,:,1), y(2,:,1)), [1.77, -0.005, 0.51, -2.95])
%!assert (nanmin (y(1,:,2), y(2,:,2)), [6.77, 4.995, NaN, 2.05])

## Test dimension indexing with vecdim in N-dimensional arrays
%!test
%! xx = repmat ([1:20;6:25], [5 2 6 3]);
%! assert (size (nanmin (xx, [], [3, 2])), [10, 1, 1, 3]);
%! assert (size (nanmin (xx, [], [1, 2])), [1, 1, 6, 3]);
%! assert (size (nanmin (xx, [], [1, 2, 4])), [1, 1, 6]);
%! assert (size (nanmin (xx, [], [1, 4, 3])), [1, 40]);
%! assert (size (nanmin (xx, [], [1, 2, 3, 4])), [1, 1]);

## Test exceeding dimensions
%!assert (nanmin (ones (2), [], 3), ones (2, 2))
%!assert (nanmin (ones (2, 2, 2), [], 99), ones (2, 2, 2))
%!assert (nanmin (magic (3), [], 3), magic (3))
%!assert (nanmin (magic (3), [], [1, 3]), [3, 1, 2])
%!assert (nanmin (magic (3), [], [1, 99]), [3, 1, 2])

## Test comparisons
%!assert (nanmin (ones (2), 3), ones (2,2))

## Test input validation
%!error <nanmin: VECDIM must contain non-repeating positive integers.> ...
%! nanmin (y, [], [1, 1, 2])
%!error <nanmin: a second output is not supported with this syntax.> ...
%! [v, idx] = nanmin(x, y, [1 2])
