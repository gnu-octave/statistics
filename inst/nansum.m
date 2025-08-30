## Copyright (C) 2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{s} =} nansum (@var{x})
## @deftypefnx {statistics} {@var{s} =} nanmax (@var{x}, @qcode{'all'})
## @deftypefnx {statistics} {@var{s} =} nanmax (@var{x}, @var{dim})
## @deftypefnx {statistics} {@var{s} =} nanmax (@var{x}, @var{vecdim})
##
## Compute the sum while ignoring NaN values.
##
## @code{@var{s} = nansum (@var{x})} returns the sum of @var{x}, after removing
## @qcode{NaN} values.  If @var{x} is a vector, a scalar value is returned.  If
## @var{x} is a matrix, a row vector of column sums is returned.  If @var{x} is
## a multidimensional array, the @code{nansum} operates along the first
## nonsingleton dimension.  If all values along a dimesion are @qcode{NaN}, the
## sum is returned returned as 0.
##
## @code{@var{s} = nansum (@var{x}, @qcode{'all'})} returns the sum of all
## elements of @var{x}, after removing @qcode{NaN} values.  It is the equivalent
## of @code{nansum (@var{x}(:))}.
##
## @code{@var{s} = nansum (@var{x}, @var{dim})} operates along the dimension
## @var{dim} of @var{x}.
##
## @code{@var{s} = nansum (@var{x}, @var{vecdim})} returns the sum over the
## dimensions specified in the vector @var{vecdim}.  Each element of
## @var{vecdim} represents a dimension of the input array @var{x} and the output
## @var{s} has length 1 in the specified operating dimensions.  The lengths of
## the other dimensions are the same for @var{x} and @var{y}.  For example, if
## @var{x} is a 2-by-3-by-4 array, then @code{nanmax (@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the maximum of the
## elements on the corresponding page of @var{x}.  If @var{vecdim} indexes all
## dimensions of @var{x}, then it is equivalent to
## @code{nanmax (@var{x}, @qcode{'all'})}.  Any dimension in @var{vecdim}
## greater than @code{ndims (@var{x})} is ignored.
##
## @seealso{sum, nanmin, nanmax}
## @end deftypefn

function s = nansum (x, dim)
  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (! isnumeric (x))
    error ("nansum: X must be numeric.");
  elseif (isempty (x))
    s = 0;
  elseif (nargin == 1)
    nanvals = isnan (x);
    x(nanvals) = 0;
    s = sum (x);
    s(all (nanvals)) = 0;
  elseif (nargin == 2 && strcmpi (dim, "all"))
    x = x(:);
    nanvals = isnan (x);
    x(nanvals) = 0;
    s = sum (x);
    s(all (nanvals)) = 0;
  else  # DIM must be a numeric scalar or vector
    if (isscalar (dim))
      if (! isnumeric (dim) || fix (dim) != dim || dim <= 0)
        error ("nansum: DIM must be a positive integer.");
      endif
      nanvals = isnan (x);
      x(nanvals) = 0;
      s = sum (x, dim);
      s(all (nanvals, dim)) = 0;
    else
      if (! isvector (dim) || any (fix (dim) != dim) || any (dim <= 0))
        error ("nansum: VECDIM must be a vector of positive integer.");
      endif
      vecdim = sort (dim);
      if (! all (diff (vecdim)))
         error ("nansum: VECDIM must contain non-repeating positive integers.");
      endif
      ## Ignore dimensions in VECDIM larger than actual array
      vecdim(find (vecdim > ndims (x))) = [];

      if (isempty (vecdim))
        s = x;
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
          x(nanvals) = 0;
          s = sum (x);
          s(all (nanvals)) = 0;

        else
          ## Permute to push vecdims to back
          perm = [remdims, vecdim];
          x = permute (x, perm);

          ## Reshape to squash all vecdims in final dimension
          sznew = [szx(remdims), prod(szx(vecdim))];
          x = reshape (x, sznew);

          ## Calculate nansum on final dimension
          dim = nremd + 1;
          nanvals = isnan (x);
          x(nanvals) = 0;
          s = sum (x, dim);
          s(all (nanvals, dim)) = 0;

          ## Inverse permute back to correct dimensions
          s = ipermute (s, perm);
        endif
      endif
    endif
  endif
endfunction

%!demo
%! ## Find the column sums for a matrix with missing values.,
%!
%! x = magic (3);
%! x([1, 4, 7:9]) = NaN
%! s = nansum (x)

%!demo
%! ## Find the row sums for a matrix with missing values.,
%!
%! x = magic (3);
%! x([1, 4, 7:9]) = NaN
%! s = nansum (x, 2)

%!demo
%! ## Find the sum of all the values in a multidimensional array
%! ## with missing values.
%!
%! x = reshape (1:30, [2, 5, 3]);
%! x([10:12, 25]) = NaN
%! s = nansum (x, "all")

%!demo
%! ## Find the sum of a multidimensional array with missing values over
%! ## multiple dimensions.
%!
%! x = reshape (1:30, [2, 5, 3]);
%! x([10:12, 25]) = NaN
%! s = nansum (x, [2, 3])

## Test output
%!assert (nansum ([]), 0)
%!assert (nansum (NaN), 0)
%!assert (nansum (NaN(3)), [0, 0, 0])
%!assert (nansum ([2 4 NaN 7]), 13)
%!assert (nansum ([2 4 NaN Inf]), Inf)
%!assert (nansum ([1 NaN 3; NaN 5 6; 7 8 NaN]), [8 13 9])
%!assert (nansum ([1 NaN 3; NaN 5 6; 7 8 NaN], 2), [4; 11; 15])
%!assert (nansum (uint8 ([2 4 1 7])), 14)
%!test
%! x = magic(3);
%! x([1 6:9]) = NaN;
%! assert (nansum (x), [7, 6, 0])
%! assert (nansum (x, 2), [1; 8; 4])
%!test
%! x = reshape(1:24, [2, 4, 3]);
%! x([5:6, 20]) = NaN;
%! assert (nansum (x, "all"), 269)
%!test
%! x = reshape(1:24,[2, 4, 3]);
%! x([5:6, 20]) = NaN;
%! assert (squeeze (nansum (x, [1, 2])), [25; 100; 144])
%! assert (nansum (x, [2, 3]), [139; 130])

## Test input validation
%!error <nansum: X must be numeric.> nansum ({3})
%!error <nansum: DIM must be a positive integer.> nansum (ones (3), 0)
%!error <nansum: DIM must be a positive integer.> nansum (ones (3), 1.5)
%!error <nansum: DIM must be a positive integer.> nansum (ones (3), 1.5)
%!error <nansum: VECDIM must be a vector of positive integer.> ...
%! nansum (ones (3, 3, 3), [2, 2.5])
%!error <nansum: VECDIM must be a vector of positive integer.> ...
%! nansum (ones (3, 3, 3), [-1, 2])
%!error <nansum: VECDIM must contain non-repeating positive integers.> ...
%! nansum (ones (3, 3, 3), [2, 2, 3])
