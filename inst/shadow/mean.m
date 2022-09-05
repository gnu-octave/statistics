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
## @deftypefn {Function File} @var{m} = mean (@var{x})
## @deftypefnx{Function File} @var{m} = mean (@var{x}, "all")
## @deftypefnx{Function File} @var{m} = mean (@var{x}, @var{dim})
## @deftypefnx{Function File} @var{m} = mean (@var{x}, @var{vecdim})
## @deftypefnx{Function File} @var{m} = mean (@dots{}, @var{outtype})
## @deftypefnx{Function File} @var{m} = mean (@dots{}, @var{nanflag})
##
## Compute the geometric mean of @var{x}.
##
## @itemize
## @item If @var{x} is a vector, then @code{geomean(@var{x})} returns the
## mean of the elements in @var{x} defined as
## @tex
## $$ {\rm mean}(x) = {1\over N} \sum_{i=1}^N x_i $$
## where $N$ is the number of elements of @var{x}.
## @end tex
##
## @ifnottex
## @example
## mean (@var{x}) = SUM_i @var{x}(i) * (1/N)
## @end example
##
## @noindent
## where @math{N} is the length of the @var{x} vector.
## @end ifnottex
##
## @item If @var{x} is a matrix, then @code{mean(@var{x})} returns a row vector
## with the mean of each columns in @var{x}.
##
## @item If @var{x} is a multidimensional array, then @code{mean(@var{x})}
## operates along the first nonsingleton dimension of @var{x}.
## @end itemize
##
## @code{mean(@var{x}, "all")} returns the geometric mean of all the elements
## in @var{x}.
##
## @code{mean(@var{x}, @var{dim})} returns the geometric mean along the
## operating dimension @var{dim} of @var{x}.
##
## @code{mean(@var{x}, @var{vecdim})} returns the geometric mean over the
## dimensions specified in the vector @var{vecdim}.  For example, if @var{x} is
## a 2-by-3-by-4 array, then @code{geomean(@var{x}, [1 2])} returns a 1-by-4
## array. Each element of the output array is the geometric mean of the elements
## on the corresponding page of @var{x}.  NOTE! @var{vecdim} MUST index at least
## N-2 dimensions of @var{x}, where @code{N = length (size (@var{x}))} and N < 8.
## If @var{vecdim} indexes all dimensions of @var{x}, then it is equivalent to
## @code{geomean(@var{x}, "all")}.
##
## @code{mean(@dots{}, @var{outtype})} returns the mean with a specified data
## type, using any of the input arguments in the previous syntaxes.
## @var{outtype} can be "default", "double", or "native".
##
## @code{mean(@dots{}, @var{nanflag})} specifies whether to exclude NaN values
## rom the calculation, using any of the input argument combinations in previous
## syntaxes. By default, geomean includes NaN values in the calculation
## (@var{nanflag} has the value "includenan").  To exclude NaN values, set the
## value of @var{nanflag} to "omitnan".
##
## @seealso{harmmean, mean}
## @end deftypefn

function m = mean (x, varargin)
  if (nargin < 1 || nargin > 4)
    print_usage ();
  endif
  if (! isnumeric (x) && ! isbool (x))
    error ("X must be either numeric or boolean vector or matrix");
  endif
  ## check for omitnan and outtype options
  omitnan = false;
  outtype = "default";
  nanflag_option = false;
  outtype_option = false;
  if nargin > 1
    for i = 1:nargin - 1
      if (ischar (varargin{i}) && strcmpi (varargin{i}, "omitnan"))
        omitnan = true;
        nanflag_option = true;
      elseif (ischar (varargin{i}) && strcmpi (varargin{i}, "includenan"))
        nanflag_option = true;
      endif
      if (ischar (varargin{i}) && ...
          any (strcmpi (varargin{i}, {"default", "double", "native"})))
        outtype = varargin{i};
        outtype_option = true;
      endif
    endfor
  endif
  ## for single input argument or with option omitnan
  if (nargin == 1 || (nargin == 2 && (nanflag_option || outtype_option)) || ...
     (nargin == 3 && nanflag_option && outtype_option))
    sz = size (x);
    dim = find (sz > 1, 1);
    if length (dim) == 0
      dim = 1;
    endif
    n = size (x, dim);
    if omitnan
      n = sum (! isnan (x), dim);
      x(isnan (x)) = 0;
    endif
    m = sum (x, dim) ./ n;
  endif
  ## for option "all"
  if ((nargin == 2 || nargin == 3 || nargin == 4) && ...
       ischar (varargin{1}) && strcmpi (varargin{1}, "all"))
    n = length (x(:));
    if omitnan
      n = length (x(! isnan (x)));
      x(isnan (x)) = 0;
    endif
    m = sum (x(:), 1) ./ n;
  endif
  ## for option DIM
  if ((nargin == 2 || nargin == 3 || nargin == 4) && ...
       isnumeric (varargin{1}) && isscalar (varargin{1}))
    dim = varargin{1};
    n = size (x, dim);
    if omitnan
      n = sum (! isnan (x), dim);
      x(isnan (x)) = 0;
    endif
    m = sum (x, dim) ./ n;
  endif
  ## resolve the pages of X when vecdim argument is provided
  if (nargin == 2 || nargin == 3 || nargin == 4) && isnumeric (varargin{1}) ...
      && ! isscalar (varargin{1}) && isvector (varargin{1})
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
      n = length (x(:));
      if omitnan
        n = length (x(! isnan (x)));
        x(isnan (x)) = 0;
      endif
      m = sum (x(:), 1) ./ n;
    ## for 1 dimension left, return column vector
    elseif length (misdim) == 1
      x = permute (x, [misdim, vecdim]);
      for i = 1:size (x, 1)
        x_vec = x(i,:,:,:,:,:,:)(:);
        if omitnan
          x_vec = x_vec(! isnan (x_vec));
        endif
        m(i) = sum (x_vec, 1) ./ length (x_vec);
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
          m(i,j) = sum (x_vec, 1) ./ length (x_vec);
        endfor
      endfor
    ## for more that 2 dimensions left, print usage
    else
      error ("vecdim must index at least N-2 dimensions of X");
    endif
  endif
  ## convert output as requested 
  switch (outtype)
    case "default"
      ## do nothing, the operators already do the right thing
    case "double"
      m = double (m);
    case "native"
      if (! islogical (x))
        m = cast (m, class (x));
      endif
  endswitch
endfunction


## Test single input and optional arguments "all", DIM, "omitnan")
%!test
%! x = [-10:10];
%! y = [x;x+5;x-5];
%! assert (mean (x), 0);
%! assert (mean (y, 2), [0, 5, -5]');
%! assert (mean (y, "all"), 0);
%! y(2,4) = NaN;
%! assert (mean (y', "omitnan"), [0 5.35 -5]);
%! z = y + 20;
%! assert (mean (z, "all"), NaN);
%! m = [20 NaN 15];
%! assert (mean (z'), m);
%! assert (mean (z', "includenan"), m);
%! m = [20 25.35 15];
%! assert (mean (z', "omitnan"), m);
%! assert (mean (z, 2, "omitnan"), m');
%! assert (mean (z, 2, "native", "omitnan"), m');
%! assert (mean (z, 2, "omitnan", "native"), m');

## Test boolean input
%!test
%! assert (mean (true, "all"), 1);
%! assert (mean (false), 0);
%! assert (mean ([true false true]), 0.6666666666666666, 4e-14);
%! assert (mean ([true false true], 1), [1 0 1]);
%! assert (mean ([true false NaN], 1), [1 0 NaN]);
%! assert (mean ([true false NaN], 2), NaN);
%! assert (mean ([true false NaN], 2, "omitnan"), 0.5);
%! assert (mean ([true false NaN], 2, "omitnan", "native"), 0.5);

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! assert (size (mean (x, [3 2])), [10 3]);
%! assert (size (mean (x, [1 2])), [6 3]);
%! assert (size (mean (x, [1 2 4])), [1 6]);
%! assert (size (mean (x, [1 4 3])), [1 40]);
%! assert (size (mean (x, [1 2 3 4])), [1 1]);

## Test results with vecdim in n-dimensional arrays and "omitnan"
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! m = repmat ([10.5;15.5], [5,3]);
%! assert (mean (x, [3 2]), m, 4e-14);
%! x(2,5,6,3) = NaN;
%! m(2,3) = NaN;
%! assert (mean (x, [3 2]), m, 4e-14);
%! m(2,3) = 15.52301255230125;
%! assert (mean (x, [3 2], "omitnan"), m, 4e-14);

## Test errors
%!error <X must be either numeric or boolean vector or matrix> mean ("char")
%!error <vecdim must index at least N-2 dimensions of X> mean ...
%!       (repmat ([1:20;6:25], [5 2 6 3 5]), [1 2])