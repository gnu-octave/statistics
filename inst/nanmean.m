## Copyright (C) 2025 Leonardo Araujos <leolca@gmail.com>
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
## @deftypefn  {statistics} {@var{s} =} nanmean (@var{x})
## @deftypefnx {statistics} {@var{s} =} nanmean (@var{x}, @qcode{'all'})
## @deftypefnx {statistics} {@var{s} =} nanmean (@var{x}, @var{dim})
## @deftypefnx {statistics} {@var{s} =} nanmean (@var{x}, @var{vecdim})
##
## Compute the mean while ignoring NaN values.
##
## @code{@var{s} = nanmean (@var{x})} returns the mean of @var{x}, after removing
## @qcode{NaN} values.  If @var{x} is a vector, a scalar value is returned.  If
## @var{x} is a matrix, a row vector of column means is returned.  If @var{x} is
## a multidimensional array, the @code{nanmean} operates along the first
## nonsingleton dimension.  If all values along a dimesion are @qcode{NaN}, the
## mean is returned returned as NaN.
##
## @code{@var{s} = nanmean (@var{x}, @qcode{'all'})} returns the mean of all
## elements of @var{x}, after removing @qcode{NaN} values.  It is the equivalent
## of @code{nanmean (@var{x}(:))}.
##
## @code{@var{s} = nanmean (@var{x}, @var{dim})} operates along the dimension
## @var{dim} of @var{x}.
##
## @code{@var{s} = nanmean (@var{x}, @var{vecdim})} returns the mean over the
## dimensions specified in the vector @var{vecdim}.  Each element of
## @var{vecdim} represents a dimension of the input array @var{x} and the output
## @var{s} has length 1 in the specified operating dimensions.  The lengths of
## the other dimensions are the same for @var{x} and @var{y}.  For example, if
## @var{x} is a 2-by-3-by-4 array, then @code{nanmean (@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the mean of the
## elements on the corresponding page of @var{x}.  If @var{vecdim} indexes all
## dimensions of @var{x}, then it is equivalent to
## @code{nanmean (@var{x}, @qcode{'all'})}.  Any dimension in @var{vecdim}
## greater than @code{ndims (@var{x})} is ignored.
##
## @seealso{mean, nanmean, NaN}
## @end deftypefn

function y = nanmean (x, dim)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (! isnumeric (x))
    error ("nanmean: X must be numeric.");
  elseif (isempty (x))
    y = NaN;
  else
    % size of the input
    sx = size(x);
    % Determine the first nonsingleton dimension to operate on
    if nargin < 2
        dim = find(sx != 1, 1);
        if isempty(dim) % scalar
            dim = 1;
        end
    else
      if (isscalar (dim))
        if (! isnumeric (dim) || fix (dim) != dim || dim <= 0)
          error ("nanmean: DIM must be a positive integer.");
        endif
      elseif (isvector (dim)) 
        if (ischar (dim))
          if strcmpi (dim, 'all')
             x = x(:);
             dim = []; % type argument 'all' not yet implemented in Octave
          else
             error ("nanmean: Invalid option.");
          endif
        else
          if (! isnumeric (dim) || any (fix (dim) != dim) || any (dim <= 0))
            error ("nanmean: VECDIM must be a vector of positive integer.");
          endif
        endif
      endif
    endif

    na = isnan (x); 
    x(na) = 0;
    if isempty (dim)
      S = sum(x);
      C = sum(!na);
    else
      if numel (dim) > 1, % correction since octave's sum gives wrong result for dim vector
        S = vsum (x, dim);
        C = vsum (!na, dim);
      else
        S = sum (x, dim);
        C = sum (!na, dim);
      endif
    endif
    y = S./C;
  endif
endfunction

function S = vsum(x, dim)
  x(isnan(x)) = 0;
  for i=numel (dim):-1:1,
    x = sum (x, dim(i));
  endfor
  S = x;
endfunction

%!demo
%! ## Find the column means for a matrix with missing values.,
%! 
%! x = magic (3);
%! x([1, 4, 7:9]) = NaN
%! y = nanmean (x)

%!demo
%! ## Find the row means for a matrix with missing values.,
%!
%! x = magic (3);
%! x([1, 4, 7:9]) = NaN
%! y = nanmean (x, 2)

%!demo
%! ## Find the mean of all the values in a multidimensional array
%! ## with missing values.
%!
%! x = reshape (1:30, [2, 5, 3]);
%! x([10:12, 25]) = NaN
%! y = nanmean (x, "all")

%!demo
%! ## Find the mean of a multidimensional array with missing values over
%! ## multiple dimensions.
%!
%! x = reshape (1:30, [2, 5, 3]);
%! x([10:12, 25]) = NaN
%! y = nanmean (x, [2, 3])

## Test output
%!assert (nanmean ([]), NaN)
%!assert (nanmean (NaN), NaN)
%!assert (nanmean (NaN(3)), [NaN, NaN, NaN])
%!assert (nanmean ([3 2 NaN 7]), 4)
%!assert (nanmean ([2 4 NaN Inf]), Inf)
%!assert (nanmean ([1 NaN 3; NaN 4 6; 7 8 NaN]), [4 6 4.5])
%!assert (nanmean ([1 NaN 3; NaN 5 6; 7 8 NaN], 2), [2; 5.5; 7.5])
%!assert (nanmean (uint8 ([2 4 1 7])), 3.5)
%!test
%! x = magic(3);
%! x([1 6:9]) = NaN;
%! assert (nanmean (x), [3.5, 3, NaN])
%! assert (nanmean (x, 2), [1; 4; 4])
%!test
%! x = reshape(1:24, [2, 4, 3]);
%! x([5:6, 20]) = NaN;
%! assert (nanmean (x, "all"), 269/21)
%!test
%! x = reshape(1:24,[2, 4, 3]);
%! x([5:6, 20]) = NaN;
%! assert (squeeze (nanmean (x, [1, 2])), [25/6; 100/8; 144/7])
%! assert (nanmean (x, [2, 3]), [139/11; 13])

