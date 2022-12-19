## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{s} = std (@var{x})
## @deftypefnx {statistics} @var{s} = std (@var{x}, @var{w})
## @deftypefnx {statistics} @var{s} = std (@var{x}, @var{w}, "all")
## @deftypefnx {statistics} @var{s} = std (@var{x}, @var{w}, @var{dim})
## @deftypefnx {statistics} @var{s} = std (@var{x}, @var{w}, @var{vecdim})
## @deftypefnx {statistics} @var{s} = std (@dots{}, @var{nanflag})
## @deftypefnx {statistics} [@var{s}, @var{m}] = std (@dots{})
##
## Compute the standard deviation of the elements of @var{x}.
##
## @itemize
## @item
## If @var{x} is a vector, then @code{var(@var{x})} returns the standard
## deviation of the elements in @var{x} defined as
## @tex
## $$ {\rm std}(x) = \sqrt{{1\over N-1} \sum_{i=1}^N |x_i - \bar x |^2} $$
## where $N$ is the number of elements of @var{x}.
##
## @end tex
## @ifnottex
##
## @example
## std (@var{s}) = sqrt ((1 / (N-1)) * SUM_i (|@var{x}(i) - mean (@var{x})|^2))
## @end example
##
## @noindent
## where @math{N} is the length of the @var{x} vector.
##
## @end ifnottex
##
## @item
## If @var{x} is a matrix, then @code{std(@var{x})} returns a row vector with
## the standard deviation of each columns in @var{x}.
##
## @item
## If @var{x} is a multidimensional array, then @code{var(@var{x})} operates
## along the first nonsingleton dimension of @var{x}.
## @end itemize
##
## @code{var (@var{x}, @var{w})} specifies a weighting scheme.   When @var{w} =
## 0 (default), the standard deviation is normalized by N-1 (population standard
## deviation), where N is the number of observations.  When @var{w} = 1, the
## standard deviation is normalized by the number of observations (sample
## standard deviation).  To use the default value you may pass an empty input
## argument [] before entering other options.
##
## @var{w} can also be a weight vector, matrix or N-D array containing
## nonnegative elements.  When @var{w} is a vector, its length must equal the
## length of the dimension over which var is operating.  When "all" flag is
## used, the length of @var{w} must equal the elements in @var{x}.  When @var{w}
## is a matrix or N-D array, its size must equal the size of @var{x}.  NaN
## values in @var{w} are treated accordingly to those in @var{x}.
##
## @code{std (@var{x}, "all")} returns the standard deviation of all the
## elements in @var{x}.
##
## @code{std (@var{x}, [], @var{dim})} returns the standard deviation along the
## operating dimension @var{dim} of @var{x}.  For @var{dim} greater than
## @code{ndims (x)}, @var{v} is returned as zeros of the same size as @var{x}
## and @code{@var{m} = @var{x}}.
##
## @code{std (@var{x}, [], @var{vecdim})} returns the standard deviation over
## the dimensions specified in the vector @var{vecdim}.  For example, if @var{x}
## is a 2-by-3-by-4 array, then @code{var (@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the standard
## deviation of the elements on the corresponding page of @var{x}.
## If @var{vecdim} indexes all dimensions of @var{x}, then it is equivalent to
## @code{std (@var{x}, "all")}.  Any dimension in @var{vecdim} greater than
## @code{ndims (@var{x})} is ignored.
##
## @code{std (@var{x}, "all")} returns the standard deviation of all the
## elements in @var{x}.  The optional flag "all" cannot be used together with
## @var{dim} or @var{vecdim} input arguments.
##
## @code{std (@dots{}, @var{nanflag})} specifies whether to exclude NaN values
## from the calculation, using any of the input argument combinations in
## previous syntaxes.  By default, NaN values are included in the calculation
## (@var{nanflag} has the value "includenan").  To exclude NaN values, set the
## value of @var{nanflag} to "omitnan".
##
## @code{[@var{s}, @var{m}] = std (@dots{})} also returns the mean of the
## elements of @var{x} used to calculate the standard deviation.  If @var{s} is
## the weighted standard deviation, then @var{m} is the weighted mean.
##
## @seealso{var, mean}
## @end deftypefn

function [y, m] = std (x, varargin)

  if (nargin < 1 || nargin > 4 || any (cellfun (@isnumeric, varargin(3:end))))
    print_usage ();
  endif

  ## Check all char arguments.
  all_flag = false;
  omitnan = false;

  for i = 1:length (varargin)
    if (ischar (varargin{i}))
      switch (varargin{i})
        case "all"
          all_flag = true;
        case "omitnan"
          omitnan = true;
        case "includenan"
          omitnan = false;
        otherwise
          print_usage ();
      endswitch
    endif
  endfor
  varargin(cellfun (@ischar, varargin)) = [];

  ## Check all numeric arguments
  if (((length (varargin) == 1) && ! (isnumeric (varargin{1}))) ...
      || ((length (varargin) == 2) && (! (isnumeric (varargin{1})) ...
          || ! (isnumeric (varargin{2})))) || (length (varargin) > 2))
    print_usage ();
  endif

  w = 0;
  weighted = false;
  if (length (varargin) > 0 && isscalar (varargin{1}))
    w = varargin{1};
    if (! (w == 0 || w == 1))
      error ("std: normalization scalar must be either 0 or 1");
    endif
  elseif (length (varargin) > 0 && numel (varargin{1}) > 1)
    weights = varargin{1};
    if (any (weights(:) < 0))
      error ("std: weights must not contain any negative values");
    endif
    weighted = true;
  endif

  vecdim = [];
  if (length (varargin) > 1)
    vecdim = varargin{2};
    if (! (isvector (vecdim) && all (vecdim)) || any (rem (vecdim, 1)))
      error ("std: DIM must be a positive integer scalar or vector");
    endif
    if (! isequal (vecdim, unique (vecdim, "stable")))
      error ("std: VECDIM must contain non-repeating positive integers");
    endif
    if (isscalar (vecdim) && vecdim > ndims (x))
      y = zeros (size (x));
      m = x;
      return;
    endif
  endif

  if (! (isnumeric (x)))
    error ("std: X must be a numeric vector or matrix");
  endif

  ## Check for conflicting input arguments
  if (all_flag && ! isempty (vecdim))
    error ("std: 'all' flag cannot be used with DIM or VECDIM options");
  endif
  if (weighted && isempty (vecdim) && ! all_flag)
    sz = size (x);
    dim = find (sz > 1, 1);
    if length (dim) == 0
      dim = 1;
    endif
    if (isvector (weights) && numel (weights) != size (x, dim))
      error ("std: weight vector does not match first operating dimension");
    endif
  elseif (weighted && isscalar (vecdim))
    if (isvector (weights) && numel (weights) != size (x, vecdim))
      error ("std: weight vector does not match given operating dimension");
    endif
  elseif (weighted && isvector (vecdim))
    if (! (isequal (size (weights), size (x))))
      error ("std: weight matrix or array does not match X in size");
    endif
  endif
  if (all_flag && weighted)
    if (isvector (weights) && numel (weights) != numel (x))
      error ("std: elements in weight vector do not match elements in X");
    endif
    if (! isvector (weights) && ! (isequal (size (weights), size (x))))
      error ("std: weight matrix or array does not match X in size");
    endif
  endif

  ## Force output for X being empty or scalar
  if (isempty (x))
    y = NaN;
    m = NaN;
    return;
  endif
  if (isnumeric (x) && isscalar (x))
    y = 0;
    m = x;
    return;
  endif

  if (length (varargin) == 0)

    ## Single numeric input argument, no dimensions or weights given.
    if (all_flag)
      x = x(:);
      if (omitnan)
        x = x(! isnan (x));
      endif
      n = length (x);
      m = sum (x) ./ n;
      y = sqrt (sum (abs (x - m) .^ 2) ./ (n - 1 + w));
    else
      sz = size (x);
      dim = find (sz > 1, 1);
      if length (dim) == 0
        dim = 1;
      endif
      n = size (x, dim);
      if (omitnan)
        n = sum (! isnan (x), dim);
        xn = isnan (x);
        x(xn) = 0;
      endif
      m = sum (x, dim) ./ n;
      dims = ones (1, ndims (x));
      dims(dim) = size (x, dim);
      m_exp = repmat (m, dims);
      if (omitnan)
        x(xn) = m_exp(xn);
      endif
      y = sqrt (sumsq (x - m_exp, dim) ./ (n - 1 + w));
    endif

  elseif (length (varargin) == 1)

    ## Two numeric input arguments, w or weights given.
    if (all_flag)
      xv = x(:);
      if (weighted)
        wv = weights(:);
      else
        wv = ones (prod (size (xv)), 1);
      endif
      wx = wv .* xv;
      if (omitnan)
        xn = wx;
        wx = wx(! isnan (xn));
        wv = wv(! isnan (xn));
        xv = xv(! isnan (xn));
      endif
      n = length (wx);
      m = sum (wx) ./ sum (wv);
      if (weighted)
        y = sqrt (sum (wv .* (abs (xv - m) .^ 2)) ./ sum (weights(:)));
      else
        y = sqrt (sum (wv .* (abs (xv - m) .^ 2)) ./ (n - 1 + w));
      endif
    else
      sz = size (x);
      dim = find (sz > 1, 1);
      if length (dim) == 0
        dim = 1;
      endif
      if (weighted)
        wv = weights(:);
      else
        wv = ones (size (x, dim), 1);
      endif
      wv = zeros (size (x)) + shiftdim (wv, 1 - dim);
      wx = wv .* x;
      n = size (wx, dim);
      if (omitnan)
        n = sum (! isnan (wx), dim);
        xn = isnan (wx);
        wx(xn) = 0;
        wv(xn) = 0;
      endif
      m = sum (wx, dim) ./ sum (wv, dim);
      dims = ones (1, ndims (wx));
      dims(dim) = size (wx, dim);
      m_exp = repmat (m, dims);
      if (omitnan)
        x(xn) = m_exp(xn);
      endif
      if (weighted)
        y = sqrt (sum (wv .* ((x - m_exp) .* (x - m_exp)), dim) ./ ...
                  sum (weights(:)));
      else
        y = sqrt (sumsq (x - m_exp, dim) ./ (n - 1 + w));
      endif
    endif

  elseif (length (varargin) == 2)

    ## Three numeric input arguments, both w or weights and dim or vecdim given.
    if (isscalar (vecdim))
      if (weighted)
        wv = weights(:);
      else
        wv = ones (size (x, vecdim), 1);
      endif
      wv = zeros (size (x)) + shiftdim (wv, 1 - vecdim);
      wx = wv .* x;
      n = size (wx, vecdim);
      if (omitnan)
        n = sum (! isnan (wx), vecdim);
        xn = isnan (wx);
        wx(xn) = 0;
        wv(xn) = 0;
      endif
      m = sum (wx, vecdim) ./ sum (wv, vecdim);
      dims = ones (1, ndims (wx));
      dims(vecdim) = size (wx, vecdim);
      m_exp = repmat (m, dims);
      if (omitnan)
        x(xn) = m_exp(xn);
      endif
      if (weighted)
        y = sqrt (sum (wv .* ((x - m_exp) .* (x - m_exp)), vecdim) ./ ...
                  sum (weights(:)));
      else
        y = sqrt (sumsq (x - m_exp, vecdim) ./ (n - 1 + w));
      endif

    else

      ## Ignore exceeding dimensions in VECDIM
      vecdim(find (vecdim > ndims (x))) = [];
      # Calculate permutation vector
      remdims = 1:ndims (x);    # all dimensions
      remdims(vecdim) = [];     # delete dimensions specified by vecdim
      nremd = numel (remdims);

      ## If all dimensions are given, it is similar to all flag
      if (nremd == 0)
        xv = x(:);
        if (weighted)
          wv = weights(:);
        else
          wv = ones (prod (size (xv)), 1);
        endif
        wx = wv .* xv;
        if (omitnan)
          xn = wx;
          wx = wx(! isnan (xn));
          wv = wv(! isnan (xn));
          xv = xv(! isnan (xn));
        endif
        n = length (wx);
        m = sum (wx) ./ sum (wv);
        if (weighted)
          y = sqrt (sum (wv .* (abs (xv - m) .^ 2)) ./ sum (weights(:)));
        else
          y = sqrt (sum (wv .* (abs (xv - m) .^ 2)) ./ (n - 1 + w));
        endif

      else

        ## Apply weights
        if (weighted)
          wv = weights;
        else
          wv = ones (size (x));
        endif
        wx = wv .* x;

        ## Permute to bring remaining dims forward
        perm = [remdims, vecdim];
        wx = permute (wx, perm);
        wv = permute (wv, perm);
        x = permute (x, perm);

        ## Reshape to put all vecdims in final dimension
        szwx = size (wx);
        sznew = [szwx(1:nremd), prod(szwx(nremd+1:end))];
        wx = reshape (wx, sznew);
        wv = reshape (wv, sznew);
        x = reshape (x, sznew);

        ## Calculate var on single, squashed dimension
        dim = nremd + 1;
        n = size (wx, dim);
        if (omitnan)
          n = sum (! isnan (wx), dim);
          xn = isnan (wx);
          wx(xn) = 0;
          wv(xn) = 0;
        endif
        m = sum (wx, dim) ./ sum (wv, dim);
        m_exp = zeros (size (wx)) + shiftdim (m, 0);
        if (omitnan)
          x(xn) = m_exp(xn);
        endif
        if (weighted)
          y = sqrt (sum (wv .* ((x - m_exp) .* (x - m_exp)), dim) ./ ...
                    sum (weights(:)));
        else
          y = sqrt (sumsq (x - m_exp, dim) ./ (n - 1 + w));
        endif

        ## Inverse permute back to correct dimensions
        y = ipermute (y, perm);
        m = ipermute (m, perm);
      endif
    endif
  endif

endfunction


## Test input validation
%!error <Invalid call to std.  Correct usage is> std ()
%!error <Invalid call to std.  Correct usage is> std (1, 2, "omitnan", 3)
%!error <Invalid call to std.  Correct usage is> std (1, 2, 3, 4)
%!error <Invalid call to std.  Correct usage is> std (1, 2, 3, 4, 5)
%!error <Invalid call to std.  Correct usage is> std (1, "foo")
%!error <Invalid call to std.  Correct usage is> std (1, [], "foo")
%!error <std: normalization scalar must be either 0 or 1> std (1, 2, "all")
%!error <std: normalization scalar must be either 0 or 1> std (1, 0.5, "all")
%!error <std: weights must not contain any negative values> ...
%! std ([1 2 3], [1 -1 0])
%!error <std: X must be a numeric vector or matrix> std ({1:5})
%!error <std: X must be a numeric vector or matrix> std ("char")
%!error <std: DIM must be a positive integer> std (1, [], ones (2,2))
%!error <std: DIM must be a positive integer> std (1, 0, 1.5)
%!error <std: DIM must be a positive integer> std (1, [], 0)
%!error <std: VECDIM must contain non-repeating positive integers> ...
%! std (repmat ([1:20;6:25], [5 2 6 3]), 0, [1 2 2 2])
%!error <std: weight vector does not match first operating dimension> ...
%! std ([1 2 3; 2 3 4], [1 3 4])
%!error <std: weight vector does not match given operating dimension> ...
%! std ([1 2 3; 2 3 4], [1 3 4], 1)
%!error <std: weight vector does not match given operating dimension> ...
%! std ([1 2 3; 2 3 4], [1 3], 2)
%!error <std: weight matrix or array does not match X in size> ...
%! std (repmat ([1:20;6:25], [5 2 6 3]), repmat ([1:20;6:25], [5 2 3]), [2 3])
%!error <std: 'all' flag cannot be used with DIM or VECDIM options> ...
%! std (1, [], 1, "all")
%!error <std: elements in weight vector do not match elements in X> ...
%! std ([1 2 3; 2 3 4], [1 3], "all")
%!error <std: weight matrix or array does not match X in size> ...
%! std (repmat ([1:20;6:25], [5 2 6 3]), repmat ([1:20;6:25], [5 2 3]), "all")

## Test single input and optional arguments "all", DIM, "omitnan")
%!test
%! x = [-10:10];
%! y = [x;x+5;x-5];
%! assert (std (x), sqrt (38.5), 1e-14);
%! assert (std (y, [], 2), sqrt ([38.5; 38.5; 38.5]), 1e-14);
%! assert (std (y, 0, 2), sqrt ([38.5; 38.5; 38.5]), 1e-14);
%! assert (std (y, 1, 2), ones (3,1) * sqrt (36.66666666666666), 1e-14);
%! assert (std (y, "all"), sqrt (54.19354838709678), 1e-14);
%! y(2,4) = NaN;
%! assert (std (y, "all"), NaN);
%! assert (std (y, "all", "includenan"), NaN);
%! assert (std (y, "all", "omitnan"), sqrt (55.01533580116342), 1e-14);
%! assert (std (y, 0, 2, "includenan"), sqrt ([38.5; NaN; 38.5]), 1e-14);
%! assert (std (y, [], 2), sqrt ([38.5; NaN; 38.5]), 1e-14);
%! assert (std (y, [], 2, "omitnan"), ...
%!         sqrt ([38.5; 37.81842105263158; 38.5]), 1e-14);

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5, 2, 6, 3]);
%! assert (size (std (x, 0, [3 2])), [10, 1, 1, 3]);
%! assert (size (std (x, 1, [1 2])), [1, 1, 6, 3]);
%! assert (size (std (x, [], [1 2 4])), [1, 1, 6]);
%! assert (size (std (x, 0, [1 4 3])), [1, 40]);
%! assert (size (std (x, [], [1 2 3 4])), [1, 1]);

## Test results with vecdim in n-dimensional arrays and "omitnan"
%!test
%! x = repmat ([1:20;6:25], [5, 2, 6, 3]);
%! v = repmat (sqrt (33.38912133891213), [10, 1, 1, 3]);
%! assert (std (x, 0, [3, 2]), v, 1e-14);
%! v = repmat (sqrt (33.250), [10, 1, 1, 3]);
%! assert (std (x, 1, [3, 2]), v, 1e-14);
%! x(2,5,6,3) = NaN;
%! v(2,1,1,3) = NaN;
%! assert (std (x, 1, [3, 2]), v, 1e-14);
%! v = repmat (sqrt (33.38912133891213), [10 1 1 3]);
%! v(2,1,1,3) = NaN;
%! assert (std (x, [], [3, 2]), v, 1e-14);
%! v(2,1,1,3) = sqrt (33.40177912169048);
%! assert (std (x, [], [3, 2], "omitnan"), v, 1e-14);

## Test mean output
%!test
%! x = repmat ([1:20;6:25], [5, 2, 6, 3]);
%! [v, m] = std (x, 0, [3 2]);
%! assert (m, mean (x, [3 2]));
%! [v, m] = std (x, 0, [1 2]);
%! assert (m, mean (x, [1 2]));
%! [v, m] = std (x, 0, [1 3 4]);
%! assert (m, mean (x, [1 3 4]));
%!test
%! x = repmat ([1:20;6:25], [5, 2, 6, 3]);
%! x(2,5,6,3) = NaN;
%! [v, m] = std (x, 0, [3 2], "omitnan");
%! assert (m, mean (x, [3 2], "omitnan"));

## Test weighted mean and variance output
%!test
%! [v, m] = std (4 * eye (2), [1, 3]);
%! assert (v, sqrt ([3, 3]), 1e-14);
%! assert (m, [1, 3]);

## Testing weights vector
%!assert (std (ones (2,2,2), [1:2], 3), [(zeros (2, 2))]);
%!assert (std (magic (3), [1:9], "all"), 2.581988897471611, 1e-14);

## Test exceeding dimensions
%!assert (std (ones (2,2), [], 3), zeros (2,2));
%!assert (std (ones (2,2,2), [], 99), zeros (2,2,2));
%!assert (std (magic (3), [], 3), zeros (3,3));
%!assert (std (magic (3), [], 1), sqrt ([7, 16, 7]));
%!assert (std (magic (3), [], [1 3]), sqrt ([7, 16, 7]));
%!assert (std (magic (3), [], [1 99]), sqrt ([7, 16, 7]));

## Test empty and scalar X
%!test
%! [v, m] = std ([]);
%! assert (v, NaN);
%! assert (m, NaN);
%! [v, m] = std (3);
%! assert (v, 0);
%! assert (m, 3);
