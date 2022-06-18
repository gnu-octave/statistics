## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{m} = geomean (@var{x})
## @deftypefnx{Function File} @var{m} = geomean (@var{x}, "all")
## @deftypefnx{Function File} @var{m} = geomean (@var{x}, @var{dim})
## @deftypefnx{Function File} @var{m} = geomean (@var{x}, @var{vecdim})
## @deftypefnx{Function File} @var{m} = geomean (@dots{}, @var{nanflag})
##
## Compute the geometric mean of @var{x}.
##
## @itemize
## @item If @var{x} is a vector, then @code{geomean(@var{x})} returns the
## geometric mean of the elements in @var{x} defined as
## @tex
## $$ {\rm geomean}(x) = \left( \prod_{i=1}^N x_i \right)^\frac{1}{N} 
## = exp \left({1\over N} \sum_{i=1}^N log x_i \right) $$
## where $N$ is the number of elements of @var{x}.
## @end tex
##
## @ifnottex
## @example
## geomean (@var{x}) = PROD_i @var{x}(i) ^ (1/N)
## @end example
##
## @noindent
## where @math{N} is the length of the @var{x} vector.
## @end ifnottex
##
## @item If @var{x} is a matrix, then @code{geomean(@var{x})} returns a row
## vector with the geometric mean of each columns in @var{x}.
##
## @item If @var{x} is a multidimensional array, then @code{geomean(@var{x})}
## operates along the first nonsingleton dimension of @var{x}.
##
## @item If @var{x} contains any negative values, then @code{geomean(@var{x})}
## returns complex values.
## @end itemize
##
## @code{geomean(@var{x}, "all")} returns the geometric mean of all the elements
## in @var{x}.
##
## @code{geomean(@var{x}, @var{dim})} returns the geometric mean along the
## operating dimension @var{dim} of @var{x}.
##
## @code{geomean(@var{x}, @var{vecdim})} returns the geometric mean over the
## dimensions specified in the vector @var{vecdim}.  For example, if @var{x} is
## a 2-by-3-by-4 array, then @code{geomean(@var{x}, [1 2])} returns a 1-by-4
## array. Each element of the output array is the geometric mean of the elements
## on the corresponding page of @var{x}.  NOTE! @var{vecdim} MUST index at least
## N-2 dimensions of @var{x}, where @code{N = length (size (@var{x}))} and N < 8.
## If @var{vecdim} indexes all dimensions of @var{x}, then it is equivalent to
## @code{geomean(@var{x}, "all")}.
##
## @code{geomean(@dots{}, @var{nanflag})} specifies whether to exclude NaN
## values from the calculation, using any of the input argument combinations in
## previous syntaxes. By default, geomean includes NaN values in the calculation
## (@var{nanflag} has the value "includenan").  To exclude NaN values, set the
## value of @var{nanflag} to "omitnan".
##
## @seealso{harmmean, mean}
## @end deftypefn

function m = geomean (x, varargin)
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif
  if (! isnumeric (x) && ! isbool (x))
    error ("X must be either numeric or boolean vector or matrix");
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
      m = exp (nansum (log (x), dim) ./ sum (! isnan (x), dim));
    else
      m = exp (sum (log (x), dim) ./ size (x, dim));
    endif
  endif
  ## for option "all"
  if ((nargin == 2 || nargin == 3) && ischar (varargin{1}) && ...
       strcmpi (varargin{1}, "all"))
    if omitnan
      m = exp (nansum (log (x(:)), 1) ./ length (x(! isnan (x))));
    else
      m = exp (sum (log (x(:)), 1) ./ length (x(:)));
    endif
  endif
  ## for option DIM
  if ((nargin == 2 || nargin == 3) && isnumeric (varargin{1}) && ...
       isscalar (varargin{1}))
    dim = varargin{1};
    if omitnan
      m = exp (nansum (log (x), dim) ./ sum (! isnan (x), dim));
    else
      m = exp (sum (log (x), dim) ./ size (x, dim));
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
        m = exp (nansum (log (x(:)), 1) ./ length (x(! isnan (x))));
      else
        m = exp (sum (log (x(:)), 1) ./ length (x(:)));
      endif
    ## for 1 dimension left, return column vector
    elseif length (misdim) == 1
      x = permute (x, [misdim, vecdim]);
      for i = 1:size (x, 1)
        x_vec = x(i,:,:,:,:,:,:)(:);
        if omitnan
          x_vec = x_vec(! isnan (x_vec));
        endif
        m(i) = exp (sum (log (x_vec), 1) ./ length (x_vec));
      endfor
    ## for 2 dimensions left, return matrix
    elseif length (misdim) == 2
      x = permute (x, [misdim, vecdim]);
      for i = 1:size (x, 1)
        for j = 1:size (x, 2)
          x_vec = x(i,j,:,:,:,:,:)(:);
          if omitnan
            x_vec = x_vec(! isnan (x_vec));
          endif
          m(i,j) = exp (sum (log (x_vec), 1) ./ length (x_vec));
        endfor
      endfor
    ## for more that 2 dimensions left, print usage
    else
      error ("vecdim must index at least N-2 dimensions of X");
    endif
  endif
endfunction


## Test single input and optional arguments "all", DIM, "omitnan")
%!test
%! x = [-10:10];
%! y = [x;x+5;x-5];
%! assert (geomean (x), 0);
%! assert (geomean (y, 2), [0, 0, 0]');
%! assert (geomean (y, "all"), 0);
%! y(2,4) = NaN;
%! assert (geomean (y', "omitnan"), [0 0 0]);
%! z = y + 20;
%! assert (geomean (z, "all"), NaN);
%! m = [19.02099329497543 NaN 13.60912525683438];
%! assert (geomean (z'), m, 4e-14);
%! assert (geomean (z', "includenan"), m, 4e-14);
%! m = [19.02099329497543 24.59957418295707 13.60912525683438];
%! assert (geomean (z', "omitnan"), m, 4e-14);
%! assert (geomean (z, 2, "omitnan"), m', 4e-14);

## Test boolean input
%!test
%! assert (geomean (true, "all"), 1);
%! assert (geomean (false), 0);
%! assert (geomean ([true false true]), 0);
%! assert (geomean ([true false true], 1), [1 0 1]);
%! assert (geomean ([true false NaN], 1), [1 0 NaN]);
%! assert (geomean ([true false NaN], 2), NaN);
%! assert (geomean ([true false NaN], 2, "omitnan"), 0);

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! assert (size (geomean (x, [3 2])), [10 3]);
%! assert (size (geomean (x, [1 2])), [6 3]);
%! assert (size (geomean (x, [1 2 4])), [1 6]);
%! assert (size (geomean (x, [1 4 3])), [1 40]);
%! assert (size (geomean (x, [1 2 3 4])), [1 1]);

## Test results with vecdim in n-dimensional arrays and "omitnan"
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! m = repmat ([8.304361203739333;14.3078118884256], [5,3]);
%! assert (geomean (x, [3 2]), m, 4e-14);
%! x(2,5,6,3) = NaN;
%! m(2,3) = NaN;
%! assert (geomean (x, [3 2]), m, 4e-14);
%! m(2,3) = 14.3292729579901;
%! assert (geomean (x, [3 2], "omitnan"), m, 4e-14);

## Test errors
%!error <X must be either numeric or boolean vector or matrix> geomean ("char")
%!error <vecdim must index at least N-2 dimensions of X> geomean ...
%!       (repmat ([1:20;6:25], [5 2 6 3 5]), [1 2])