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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
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
## @qcode{NaN}, the median is returned as @qcode{NaN}.  As in @code{median},
## @var{m} takes the class of @var{x}, so an integer median is rounded.
##
## @code{@var{m} = nanmedian (@var{x}, @qcode{'all'})} returns the median of all
## elements of @var{x}, after removing @qcode{NaN} values.  It is the equivalent
## of @code{nanmedian (@var{x}(:))}.
##
## @code{@var{m} = nanmedian (@var{x}, @var{dim})} operates along the dimension
## @var{dim} of @var{x}.  A logical @var{dim} or @var{vecdim} is an error, as in
## MATLAB.
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
    error ("nanmedian: X must be numeric or logical.");
  endif
  ## Operating dimension(s)
  if (nargin < 2)
    m = median (x, 'omitnan');
    return;
  elseif (ischar (dim) && strcmpi (dim, 'all'))
    dim = 'all';
  elseif (isscalar (dim))
    if (! isnumeric (dim) || ! isfinite (dim) || fix (dim) != dim || dim <= 0)
      error ("nanmedian: DIM must be a positive integer.");
    endif
  elseif (isnumeric (dim) && isvector (dim))
    if (any (! isfinite (dim)) || any (fix (dim) != dim) || any (dim <= 0))
      error ("nanmedian: VECDIM must be a vector of positive integers.");
    endif
    if (! all (diff (sort (dim))))
      error ("nanmedian: VECDIM must contain non-repeating positive integers.");
    endif
  else
    error ("nanmedian: invalid DIM argument.");
  endif

  m = median (x, dim, 'omitnan');

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
%!assert_equal (nanmedian (zeros (0, 3)), [NaN, NaN, NaN])
%!assert_equal (nanmedian (zeros (3, 0)), zeros (1, 0))
%!assert_equal (nanmedian ([], 1), zeros (1, 0))
%!assert_equal (nanmedian ([], 2), zeros (0, 1))
%!assert_equal (size (nanmedian (ones (2, 0, 3, 2), 2)), [2, 1, 3, 2])
%!assert_equal (nanmedian (ones (2, 0, 3, 2), 'all'), NaN)
%!assert_equal (size (nanmedian (ones (2, 0, 3, 2), [1, 2])), [1, 1, 3, 2])
%!assert_equal (nanmedian (NaN), NaN)
%!assert_equal (nanmedian (5), 5)
%!assert_equal (nanmedian ([2, 4, NaN, 8]), 4)
%!assert_equal (nanmedian ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12]), [5.5, 8, 9])
%!assert_equal (nanmedian ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12], 2), ...
%!              [1.5; 5; 8; 11])
%!assert_equal (nanmedian ([1 NaN; NaN NaN; 3 NaN]), [2, NaN])
%!assert_equal (nanmedian (reshape (1:12, [2, 3, 2])(:)), 6.5)
%!assert_equal (nanmedian (single ([1 2 NaN 4])), single (2))
%!assert_equal (nanmedian (int8 ([2 5; 3 9])), int8 ([3, 7]))
%!assert_equal (nanmedian (uint8 ([200 250])), uint8 (225))
%!assert_equal (nanmedian (intmax ('int64') - int64 ([1 3])), ...
%!              intmax ('int64') - int64 (2))
%!assert_equal (nanmedian ([true false]), true)
%!assert_equal (nanmedian (int8 ([])), int8 (0))
%!assert_equal (nanmedian (sparse ([1 NaN 3])), sparse (2))
%!assert_equal (nanmedian (int8 ([1 2; 3 4]), 'all'), int8 (3))
%!assert_equal (nanmedian (single ([1 NaN; 3 4]), [1 2]), single (3))
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
%!error <nanmedian: X must be numeric or logical.> nanmedian ({3})
%!error <nanmedian: DIM must be a positive integer.> nanmedian (ones (3), 0)
%!error <nanmedian: DIM must be a positive integer.> nanmedian (ones (3), 1.5)
%!error <nanmedian: DIM must be a positive integer.> nanmedian (ones (3), true)
%!error <nanmedian: DIM must be a positive integer.> nanmedian (ones (3), Inf)
%!error <nanmedian: VECDIM must be a vector of positive integers.> ...
%! nanmedian (ones (3), [1, Inf])
%!error <nanmedian: VECDIM must contain non-repeating positive integers.> ...
%! nanmedian (ones (3, 3, 3), [2, 2, 3])
