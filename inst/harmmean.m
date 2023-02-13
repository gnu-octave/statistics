## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{m} =} harmmean (@var{x})
## @deftypefnx {statistics} {@var{m} =} harmmean (@var{x}, "all")
## @deftypefnx {statistics} {@var{m} =} harmmean (@var{x}, @var{dim})
## @deftypefnx {statistics} {@var{m} =} harmmean (@var{x}, @var{vecdim})
## @deftypefnx {statistics} {@var{m} =} harmmean (@dots{}, @var{nanflag})
##
## Compute the harmonic mean of @var{x}.
##
## @itemize
## @item If @var{x} is a vector, then @code{harmmean(@var{x})} returns the
## harmonic mean of the elements in @var{x} defined as
## @tex
## $$ {\rm harmmean}(x) = \frac{N}{\sum_{i=1}^N \frac{1}{x_i}} $$
## where $N$ is the number of elements of @var{x}.
## @end tex
##
## @ifnottex
## @example
## harmmean (@var{x}) = N / SUM_i @var{x}(i)^-1
## @end example
##
## @noindent
## where @math{N} is the length of the @var{x} vector.
## @end ifnottex
##
## @item If @var{x} is a matrix, then @code{harmmean(@var{x})} returns a row
## vector with the harmonic mean of each columns in @var{x}.
##
## @item If @var{x} is a multidimensional array, then @code{harmmean(@var{x})}
## operates along the first nonsingleton dimension of @var{x}.
##
## @item If @var{x} contains any negative values, then @code{harmmean(@var{x})}
## returns an error.
## @end itemize
##
## @code{harmmean(@var{x}, "all")} returns the harmonic mean of all the elements
## in @var{x}.  If @var{x} contains any 0, then the returned value is 0.
##
## @code{harmmean(@var{x}, @var{dim})} returns the harmonic mean along the
## operating dimension @var{dim} of @var{x}.  Calculating the harmonic mean of
## any subarray containing any 0 will return 0.
##
## @code{harmmean(@var{x}, @var{vecdim})} returns the harmonic mean over the
## dimensions specified in the vector @var{vecdim}.  For example, if @var{x} is
## a 2-by-3-by-4 array, then @code{harmmean(@var{x}, [1 2])} returns a 1-by-4
## array. Each element of the output array is the harmonic mean of the elements
## on the corresponding page of @var{x}.  NOTE! @var{vecdim} MUST index at least
## N-2 dimensions of @var{x}, where @code{N = length (size (@var{x}))} and N < 8.
## If @var{vecdim} indexes all dimensions of @var{x}, then it is equivalent to
## @code{geomean(@var{x}, "all")}.
##
## @code{harmmean(@dots{}, @var{nanflag})} specifies whether to exclude NaN
## values from the calculation, using any of the input argument combinations in
## previous syntaxes. By default, harmmean includes NaN values in the calculation
## (@var{nanflag} has the value "includenan").  To exclude NaN values, set the
## value of @var{nanflag} to "omitnan".
##
## @seealso{geomean, mean}
## @end deftypefn

function m = harmmean (x, varargin)
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif
  if (! isnumeric (x) && ! isbool (x))
    error ("X must be either numeric or boolean vector or matrix");
  endif
  if (any (x(:) < 0))
    error ("X must not contain any negative values");
  endif
  ## check for omitnan option
  omitnan = false;
  nanflag_option = false;
  if nargin > 1
    for i = 1:nargin - 1
      if (ischar (varargin{i}) && strcmpi (varargin{i}, "omitnan"))
        omitnan = true;
        nanflag_option = true;
      elseif (ischar (varargin{i}) && strcmpi (varargin{i}, "includenan"))
        nanflag_option = true;
      endif
    endfor
  endif
  ## for single input argument or with option omitnan
  if (nargin == 1 || (nargin == 2 && nanflag_option))
    sz = size (x);
    dim = find (sz > 1, 1);
    if length (dim) == 0
      dim = 1;
    endif
    if omitnan
      m = sum (! isnan (x), dim) ./ nansum (1 ./ x, dim);
      m(m == Inf) = 0;
    else
      m = size (x, dim) ./ sum (1 ./ x, dim);
      m(m == Inf) = 0;
    endif
  endif
  ## for option "all"
  if ((nargin == 2 || nargin == 3) && ischar (varargin{1}) && ...
       strcmpi (varargin{1}, "all"))
    if omitnan
      m = length (x(! isnan (x))) ./ nansum (1 ./ x(:), 1);
      m(m == Inf) = 0;
    else
      m = length (x(:)) ./ sum (1 ./ x(:), 1);
      m(m == Inf) = 0;
    endif
  endif
  ## for option DIM
  if ((nargin == 2 || nargin == 3) && isnumeric (varargin{1}) && ...
       isscalar (varargin{1}))
    dim = varargin{1};
    if omitnan
      m = sum (! isnan (x), dim) ./ nansum (1 ./ x, dim);
      m(m == Inf) = 0;
    else
      m = size (x, dim) ./ sum (1 ./ x, dim);
      m(m == Inf) = 0;
    endif
  endif
  ## resolve the pages of X when vecdim argument is provided
  if (nargin == 2 || nargin == 3 ) && isnumeric (varargin{1}) && ...
      ! isscalar (varargin{1}) && isvector (varargin{1})
    vecdim = varargin{1};
    sz = size (x);
    ndims = length (sz);
    misdim = [1:ndims];
    ## keep remaining dimensions
    for i = 1:length (vecdim)
      misdim(misdim == vecdim(i)) = [];
    endfor
    ## if all dimensions are given, compute x(:)
    if length (misdim) == 0
      if omitnan
        m = length (x(! isnan (x))) ./ nansum (1 ./ x(:), 1);
        m(m == Inf) = 0;
      else
        m = length (x(:)) ./ sum (1 ./ x(:), 1);
        m(m == Inf) = 0;
      endif
    ## for 1 dimension left, return column vector
    elseif length (misdim) == 1
      x = permute (x, [misdim, vecdim]);
      for i = 1:size (x, 1)
        x_vec = x(i,:,:,:,:,:,:)(:);
        if omitnan
          x_vec = x_vec(! isnan (x_vec));
        endif
        m(i) = length (x_vec) ./ sum (1 ./ x_vec, 1);
      endfor
      m(m == Inf) = 0;
    ## for 2 dimensions left, return matrix
    elseif length (misdim) == 2
      x = permute (x, [misdim, vecdim]);
      for i = 1:size (x, 1)
        for j = 1:size (x, 2)
          x_vec = x(i,j,:,:,:,:,:)(:);
          if omitnan
            x_vec = x_vec(! isnan (x_vec));
          endif
          m(i,j) = length (x_vec) ./ sum (1 ./ x_vec, 1);
        endfor
      endfor
      m(m == Inf) = 0;
    ## for more that 2 dimensions left, print usage
    else
      error ("vecdim must index at least N-2 dimensions of X");
    endif
  endif
endfunction


## Test single input and optional arguments "all", DIM, "omitnan")
%!test
%! x = [0:10];
%! y = [x;x+5;x+10];
%! assert (harmmean (x), 0);
%! m = [0 8.907635160795225 14.30854471766802];
%! assert (harmmean (y, 2), m', 4e-14);
%! assert (harmmean (y, "all"), 0);
%! y(2,4) = NaN;
%! m(2) = 9.009855936313949;
%! assert (harmmean (y', "omitnan"), m, 4e-14);
%! z = y + 20;
%! assert (harmmean (z, "all"), NaN);
%! m = [24.59488458841874 NaN 34.71244385944397];
%! assert (harmmean (z'), m, 4e-14);
%! assert (harmmean (z', "includenan"), m, 4e-14);
%! m = [24.59488458841874 29.84104075528276 34.71244385944397];
%! assert (harmmean (z', "omitnan"), m, 4e-14);
%! assert (harmmean (z, 2, "omitnan"), m', 4e-14);

## Test boolean input
%!test
%! assert (harmmean (true, "all"), 1);
%! assert (harmmean (false), 0);
%! assert (harmmean ([true false true]), 0);
%! assert (harmmean ([true false true], 1), [1 0 1]);
%! assert (harmmean ([true false NaN], 1), [1 0 NaN]);
%! assert (harmmean ([true false NaN], 2), NaN);
%! assert (harmmean ([true false NaN], 2, "omitnan"), 0);

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! assert (size (harmmean (x, [3 2])), [10 3]);
%! assert (size (harmmean (x, [1 2])), [6 3]);
%! assert (size (harmmean (x, [1 2 4])), [1 6]);
%! assert (size (harmmean (x, [1 4 3])), [1 40]);
%! assert (size (harmmean (x, [1 2 3 4])), [1 1]);

## Test results with vecdim in n-dimensional arrays and "omitnan"
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! m = repmat ([5.559045930488016;13.04950789021461], [5,3]);
%! assert (harmmean (x, [3 2]), m, 4e-14);
%! x(2,5,6,3) = NaN;
%! m(2,3) = NaN;
%! assert (harmmean (x, [3 2]), m, 4e-14);
%! m(2,3) = 13.06617961315406;
%! assert (harmmean (x, [3 2], "omitnan"), m, 4e-14);

## Test errors
%!error <X must be either numeric or boolean vector or matrix> harmmean ("char")
%!error <vecdim must index at least N-2 dimensions of X> harmmean ...
%!       (repmat ([1:20;6:25], [5 2 6 3 5]), [1 2])
