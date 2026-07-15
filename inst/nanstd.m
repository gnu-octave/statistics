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
## @deftypefn  {statistics} {@var{s} =} nanstd (@var{x})
## @deftypefnx {statistics} {@var{s} =} nanstd (@var{x}, @var{w})
## @deftypefnx {statistics} {@var{s} =} nanstd (@var{x}, @var{w}, @qcode{'all'})
## @deftypefnx {statistics} {@var{s} =} nanstd (@var{x}, @var{w}, @var{dim})
## @deftypefnx {statistics} {@var{s} =} nanstd (@var{x}, @var{w}, @var{vecdim})
##
## Compute the standard deviation while ignoring NaN values.
##
## @code{@var{s} = nanstd (@var{x})} returns the standard deviation of @var{x},
## after removing @qcode{NaN} values.  If @var{x} is a vector, a scalar value is
## returned.  If @var{x} is a matrix, a row vector of column standard deviations
## is returned.  If @var{x} is a multidimensional array, @code{nanstd} operates
## along the first nonsingleton dimension.  If a dimension contains fewer than
## two non-@qcode{NaN} values, the standard deviation is returned as 0 for a
## single value and as @qcode{NaN} when all values are @qcode{NaN}.
##
## @code{@var{s} = nanstd (@var{x}, @var{w})} specifies the normalization.  When
## @var{w} is 0 (default), the standard deviation is normalized by @math{N-1},
## where @math{N} is the number of non-@qcode{NaN} observations.  When @var{w}
## is 1, it is normalized by @math{N}.  @var{w} may also be a vector of
## nonnegative weights whose length matches the operating dimension, in which
## case the weighted standard deviation normalized by the sum of the weights is
## returned.
##
## @code{@var{s} = nanstd (@var{x}, @var{w}, @qcode{'all'})} returns the
## standard deviation of all elements of @var{x}, after removing @qcode{NaN}
## values.  Use an empty value, @code{@var{w} = []}, to pass the default
## normalization.
##
## @code{@var{s} = nanstd (@var{x}, @var{w}, @var{dim})} operates along the
## dimension @var{dim} of @var{x}.
##
## @code{@var{s} = nanstd (@var{x}, @var{w}, @var{vecdim})} returns the standard
## deviation over the dimensions specified in the vector @var{vecdim}.  A weight
## vector is not supported together with @qcode{'all'} or @var{vecdim}.  Any
## dimension in @var{vecdim} greater than @code{ndims (@var{x})} is ignored.
##
## @seealso{std, nanvar, nanmean, nansum}
## @end deftypefn

function s = nanstd (x, varargin)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  elseif (! isnumeric (x) && ! islogical (x))
    error ("nanstd: X must be numeric.");
  endif
  s = sqrt (nanvar (x, varargin{:}));

endfunction

%!demo
%! ## Find the column standard deviations for a matrix with missing values.
%!
%! x = magic (3);
%! x([1, 6:9]) = NaN
%! s = nanstd (x)

%!demo
%! ## Find the row standard deviations, normalized by N instead of N-1.
%!
%! x = magic (3);
%! x([1, 6:9]) = NaN
%! s = nanstd (x, 1, 2)

## Test output
%!assert_equal (nanstd ([]), NaN)
%!assert_equal (nanstd (NaN), NaN)
%!assert_equal (nanstd (5), 0)
%!assert (nanstd ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12]), ...
%!        sqrt ([15, 21, 9]), 1e-14)
%!assert (nanstd ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12], 1), ...
%!        sqrt ([11.25, 14, 6]), 1e-14)
%!assert (nanstd ([1 2 NaN; 4 NaN 6; 7 8 9; 10 11 12], 0, 2), ...
%!        sqrt ([0.5; 2; 1; 1]), 1e-14)
%!assert_equal (nanstd (NaN (2, 3)), [NaN, NaN, NaN])

## Test input validation
%!error <Invalid call to nanstd> nanstd ()
%!error <nanstd: X must be numeric.> nanstd ({3})
