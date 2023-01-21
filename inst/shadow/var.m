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
## @deftypefn  {statistics} @var{v} = var (@var{x})
## @deftypefnx {statistics} @var{v} = var (@var{x}, @var{w})
## @deftypefnx {statistics} @var{v} = var (@var{x}, @var{w}, "all")
## @deftypefnx {statistics} @var{v} = var (@var{x}, @var{w}, @var{dim})
## @deftypefnx {statistics} @var{v} = var (@var{x}, @var{w}, @var{vecdim})
## @deftypefnx {statistics} @var{v} = var (@dots{}, @var{nanflag})
## @deftypefnx {statistics} [@var{v}, @var{m}] = var (@dots{})
##
## Compute the variance of the elements of @var{x}.
##
## @itemize
## @item
## If @var{x} is a vector, then @code{var(@var{x})} returns the variance of the
## elements in @var{x} defined as
## @tex
## $$ {\rm var}(x) = {1\over N-1} \sum_{i=1}^N |x_i - \bar x |^2 $$
## where $N$ is the number of elements of @var{x}.
##
## @end tex
## @ifnottex
##
## @example
## var (@var{v}) = (1 / (N-1)) * SUM_i (|@var{x}(i) - mean (@var{x})|^2)
## @end example
##
## @noindent
## where @math{N} is the length of the @var{x} vector.
##
## @end ifnottex
##
## @item
## If @var{x} is a matrix, then @code{var (@var{x})} returns a row vector with
## the variance of each column in @var{x}.
##
## @item
## If @var{x} is a multi-dimensional array, then @code{var (@var{x})} operates
## along the first non-singleton dimension of @var{x}.
## @end itemize
##
## @code{var (@var{x}, @var{w})} specifies a weighting scheme.  When @var{w} = 0
## (default), the variance is normalized by N-1 (population variance) where N is
## the number of observations.  When @var{w} = 1, the variance is normalized by
## the number of observations (sample variance).  To use the default value you
## may pass an empty input argument [] before entering other options.
##
## @var{w} can also be an array of non-negative numbers.  When @var{w} is a
## vector, it must have the same length as the number of elements in the
## operating dimension of @var{x}.  If @var{w} is a matrix or n-D array, or the
## operating dimension is supplied as a @var{vecdim} or "all", @var{w} must be
## the same size as @var{x}.  NaN values are permitted in @var{w}, will be
## multiplied with the associated values in @var{x}, and can be excluded by the
## @var{nanflag} option.
##
## @code{var (@var{x}, [], @var{dim})} returns the variance along the operating
## dimension @var{dim} of @var{x}.  For @var{dim} greater than
## @code{ndims (@var{x})} @var{v} is returned as zeros of the same size as
## @var{x} and @var{m} = @var{x}.
##
## @code{var (@var{x}, [], @var{vecdim})} returns the variance over the
## dimensions specified in the vector @var{vecdim}.  For example, if @var{x}
## is a 2-by-3-by-4 array, then @code{var (@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the variance of the
## elements on the corresponding page of @var{x}.  If @var{vecdim} indexes all
## dimensions of @var{x}, then it is equivalent to @code{var (@var{x}, "all")}.
## Any dimension in @var{vecdim} greater than @code{ndims (@var{x})} is ignored.
##
## @code{var (@var{x}, "all")} returns the variance of all the elements in
## @var{x}.  The optional flag "all" cannot be used together with @var{dim} or
## @var{vecdim} input arguments.
##
## @code{var (@dots{}, @var{nanflag})} specifies whether to exclude NaN values
## from the calculation using any of the input argument combinations in previous
## syntaxes.  The default value for @var{nanflag} is "includenan", and keeps NaN
## values in the calculation. To exclude NaN values, set the value of
## @var{nanflag} to "omitnan".
##
## @code{[@var{v}, @var{m}] = var (@dots{})} also returns the mean of the
## elements of @var{x} used to calculate the variance.  If @var{v} is the
## weighted variance, then @var{m} is the weighted mean.
##
## @seealso{std, mean}
## @end deftypefn

function [v, m] = var (x, varargin)

  if (nargin < 1 || nargin > 4)
    print_usage ();
  endif

  ## initialize variables
  all_flag = false;
  omitnan = false;
  nvarg = numel (varargin);
  varg_chars = cellfun ('ischar', varargin);

  ## Check all char arguments.
  if (nvarg == 3 && ! varg_chars(3))
    print_usage ();
  endif

  if (any (varg_chars))
    for i = varargin(varg_chars)
      switch (i{:})
        case "all"
          all_flag = true;
        case "omitnan"
          omitnan = true;
        case "includenan"
          omitnan = false;
        otherwise
          print_usage ();
      endswitch
    endfor
    varargin(varg_chars) = [];
    nvarg = numel (varargin);
  endif

  w = 0;
  weighted = false; # true if weight vector applied to individual values of x
  vecdim = [];
  vecdim_scalar_vector = [false, false]; # false false for empty vecdim
  szx = size (x);

  ## Check numeric arguments
  if (! (isnumeric (x)))
    error ("var: X must be a numeric vector or matrix");
  endif
  if (isa (x, "single"))
    outtype = "single";
  else
    outtype = "double";
  endif

  if (nvarg > 0)
    if (nvarg > 2 || any (! cellfun ('isnumeric', varargin)))
      print_usage ();
    endif
    ## Process weight input
    if (any (varargin{1} < 0))
      error ("var: weights must not contain any negative values");
    endif
    if (isscalar (varargin{1}))
      w = varargin{1};
      if (! (w == 0 || w == 1) && ! isscalar (x))
        error ("var: normalization scalar must be either 0 or 1");
      endif
    elseif (numel (varargin{1}) > 1)
      weights = varargin{1};
      weighted = true;
    endif
    if (nvarg > 1)
    ## Process dimension input
      vecdim = varargin{2};
      if (! (isvector (vecdim) && all (vecdim)) || any (rem (vecdim, 1)))
        error ("var: DIM must be a positive integer scalar or vector");
      endif
      if (! isequal (vecdim, unique (vecdim, "stable")))
        error ("var: VECDIM must contain non-repeating positive integers");
      endif
      if (! isempty (x) && isscalar (vecdim) && vecdim > ndims (x))
        v = zeros (szx, outtype);
        vn = ! isfinite (x);
        v(vn) = NaN;
        m = x;
        return;
      endif
    endif
  endif

  ## Check for conflicting input arguments
  if (all_flag && ! isempty (vecdim))
    error ("var: 'all' flag cannot be used with DIM or VECDIM options");
  endif
  if (weighted && isempty (vecdim) && ! all_flag)
    dim = find (szx > 1, 1);
    if length (dim) == 0
      dim = 1;
    endif
    if (isvector (weights) && numel (weights) != szx(dim))
      error ("var: weight vector does not match first operating dimension");
    endif
    if (! isvector (weights) && numel (weights) != szx(dim))
      error ("var: weight matrix or array does not match X in size");
    endif
  elseif (weighted && isscalar (vecdim))
    if (isvector (weights) && numel (weights) != szx(vecdim))
      error ("var: weight vector does not match given operating dimension");
    endif
  elseif (weighted && isvector (vecdim))
    if (! (isequal (size (weights), szx)))
      error ("var: weight matrix or array does not match X in size");
    endif
  endif
  if (all_flag && weighted)
    if (isvector (weights) && numel (weights) != numel (x))
      error ("var: elements in weight vector do not match elements in X");
    endif
    if (! isvector (weights) && ! (isequal (size (weights), szx)))
      error ("var: weight matrix or array does not match X in size");
    endif
  endif

  ## Force output for X being empty or scalar
  if (isempty (x))
    if (isempty (vecdim) && all ((szx) == 0))
      v = NaN;
      m = NaN;
      return;
    elseif (isempty (vecdim) && ndims (x) == 2)
      v = NaN;
      m = NaN;
      return;
    endif
    if (isscalar (vecdim))
      nanvec = szx;
      nanvec(vecdim) = 1;
      v = NaN(nanvec);
      m = NaN(nanvec);
      return;
    endif
  endif
  if (isscalar (x))
    if (isfinite (x))
      v = cast (0, outtype);
    else
      v = cast (NaN, outtype);
    endif
    m = x;
    return;
  endif

  if (nvarg == 0)

    ## Single numeric input argument, no dimensions or weights given.
    if (all_flag)
      x = x(:);
      if (omitnan)
        x = x(! isnan (x));
      endif
      n = length (x);
      m = sum (x) ./ n;
      v = sum (abs (x - m) .^ 2) ./ (n - 1 + w);
      if (n == 1)
        v = 0;
      endif
    else
      dim = find (szx > 1, 1);
      if length (dim) == 0
        dim = 1;
      endif
      n = szx(dim);
      if (omitnan)
        n = sum (! isnan (x), dim);
        xn = isnan (x);
        x(xn) = 0;
      endif
      m = sum (x, dim) ./ n;
      dims = ones (1, ndims (x));
      dims(dim) = szx(dim);
      m_exp = repmat (m, dims);
      if (omitnan)
        x(xn) = m_exp(xn);
      endif
      v = sumsq (x - m_exp, dim) ./ (n - 1 + w);
      if (numel (n) == 1)
        divby0 = repmat (n, size (v)) == 1;
      else
         divby0 = n == 1;
      endif
      v(divby0) = 0;
    endif

  elseif (nvarg == 1)

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
        v = sum (wv .* (abs (xv - m) .^ 2)) ./ sum (weights(:));
      else
        v = sum (wv .* (abs (xv - m) .^ 2)) ./ (n - 1 + w);
        if (n == 1)
          v = 0;
        endif
      endif
    else
      dim = find (szx > 1, 1);
      if length (dim) == 0
        dim = 1;
      endif
      if (weighted)
        wv = weights(:);
      else
        wv = ones (szx(dim), 1);
      endif
      wv = zeros (szx) + shiftdim (wv, 1 - dim);
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
        v = sum (wv .* ((x - m_exp) .* (x - m_exp)), dim) ./ sum (weights(:));
      else
        v = sumsq (x - m_exp, dim) ./ (n - 1 + w);
        if (numel (n) == 1)
          divby0 = repmat (n, size (v)) == 1;
        else
           divby0 = n == 1;
        endif
        v(divby0) = 0;
      endif
    endif

  elseif (nvarg == 2)

    ## Three numeric input arguments, both w or weights and dim or vecdim given.
    if (isscalar (vecdim))
      if (weighted)
        wv = weights(:);
      else
        wv = ones (szx(vecdim), 1);
      endif
      wv = zeros (szx) + shiftdim (wv, 1 - vecdim);
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
        v = sum (wv .* ((x - m_exp) .* (x - m_exp)), vecdim) ./ ...
            sum (weights(:));
      else
        v = sumsq (x - m_exp, vecdim);
        vn = isnan (v);
        v = v ./ (n - 1 + w);
        if (numel (n) == 1)
          divby0 = repmat (n, size (v)) == 1;
        else
          divby0 = n == 1;
        endif
        v(divby0) = 0;
        v(vn) = NaN;
      endif

    else

      ## Ignore exceeding dimensions in VECDIM
      vecdim(find (vecdim > ndims (x))) = [];
      ## Calculate permutation vector
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
          v = sum (wv .* (abs (xv - m) .^ 2)) ./ sum (weights(:));
        else
          v = sum (wv .* (abs (xv - m) .^ 2)) ./ (n - 1 + w);
          if (n == 1)
            v = 0;
          endif
        endif

      else

        ## Apply weights
        if (weighted)
          wv = weights;
        else
          wv = ones (szx);
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
          v = sum (wv .* ((x - m_exp) .* (x - m_exp)), dim) ./ sum (weights(:));
        else
          v = sumsq (x - m_exp, dim) ./ (n - 1 + w);
          if (numel (n) == 1)
            divby0 = repmat (n, size (v)) == 1;
          else
             divby0 = n == 1;
          endif
          v(divby0) = 0;
        endif

        ## Inverse permute back to correct dimensions
        v = ipermute (v, perm);
        m = ipermute (m, perm);
      endif
    endif
  endif

  ## Preserve class type
  v = cast (v, outtype);
  m = cast (m, outtype);

endfunction


## Test single input and optional arguments "all", DIM, "omitnan")
%!test
%! x = [-10:10];
%! y = [x;x+5;x-5];
%! assert (var (x), 38.5);
%! assert (var (y, [], 2), [38.5; 38.5; 38.5]);
%! assert (var (y, 0, 2), [38.5; 38.5; 38.5]);
%! assert (var (y, 1, 2), ones (3,1) * 36.66666666666666, 1e-14);
%! assert (var (y, "all"), 54.19354838709678, 1e-14);
%! y(2,4) = NaN;
%! assert (var (y, "all"), NaN);
%! assert (var (y, "all", "includenan"), NaN);
%! assert (var (y, "all", "omitnan"), 55.01533580116342, 1e-14);
%! assert (var (y, 0, 2, "includenan"), [38.5; NaN; 38.5]);
%! assert (var (y, [], 2), [38.5; NaN; 38.5]);
%! assert (var (y, [], 2, "omitnan"), [38.5; 37.81842105263158; 38.5], 1e-14);

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5, 2, 6, 3]);
%! assert (size (var (x, 0, [3 2])), [10, 1, 1, 3]);
%! assert (size (var (x, 1, [1 2])), [1, 1, 6, 3]);
%! assert (size (var (x, [], [1 2 4])), [1, 1, 6]);
%! assert (size (var (x, 0, [1 4 3])), [1, 40]);
%! assert (size (var (x, [], [1 2 3 4])), [1, 1]);

## Test results with vecdim in n-dimensional arrays and "omitnan"
%!test
%! x = repmat ([1:20;6:25], [5, 2, 6, 3]);
%! v = repmat (33.38912133891213, [10, 1, 1, 3]);
%! assert (var (x, 0, [3, 2]), v, 1e-14);
%! v = repmat (33.250, [10, 1, 1, 3]);
%! assert (var (x, 1, [3, 2]), v, 1e-14);
%! x(2,5,6,3) = NaN;
%! v(2,1,1,3) = NaN;
%! assert (var (x, 1, [3, 2]), v, 4e-14);
%! v = repmat (33.38912133891213, [10 1 1 3]);
%! v(2,1,1,3) = NaN;
%! assert (var (x, [], [3, 2]), v, 4e-14);
%! v(2,1,1,3) = 33.40177912169048;
%! assert (var (x, [], [3, 2], "omitnan"), v, 4e-14);

## Test mean output
%!test
%! x = repmat ([1:20;6:25], [5, 2, 6, 3]);
%! [v, m] = var (x, 0, [3 2]);
%! assert (m, mean (x, [3 2]));
%! [v, m] = var (x, 0, [1 2]);
%! assert (m, mean (x, [1 2]));
%! [v, m] = var (x, 0, [1 3 4]);
%! assert (m, mean (x, [1 3 4]));
%!test
%! x = repmat ([1:20;6:25], [5, 2, 6, 3]);
%! x(2,5,6,3) = NaN;
%! [v, m] = var (x, 0, [3 2], "omitnan");
%! assert (m, mean (x, [3 2], "omitnan"));

## Test weighted mean and variance output
%!test
%! [v, m] = var (4 * eye (2), [1, 3]);
%! assert (v, [3, 3]);
%! assert (m, [1, 3]);

## Testing weights vector
%!assert (var (ones (2,2,2), [1:2], 3), [(zeros (2, 2))]);
%!assert (var (magic (3), [1:9], "all"), 6.666666666666667, 1e-14);

## Test exceeding dimensions
%!assert (var (ones (2,2), [], 3), zeros (2,2));
%!assert (var (ones (2,2,2), [], 99), zeros (2,2,2));
%!assert (var (magic (3), [], 3), zeros (3,3));
%!assert (var (magic (3), [], 1), [7, 16, 7]);
%!assert (var (magic (3), [], [1 3]), [7, 16, 7]);
%!assert (var (magic (3), [], [1 99]), [7, 16, 7]);

## Test empty and scalar X
%!test
%! [v, m] = var ([]);
%! assert (v, NaN);
%! assert (m, NaN);
%! [v, m] = var (3);
%! assert (v, 0);
%! assert (m, 3);

####
#### BISTs from core Octave
####

%!assert (var (13), 0)
%!assert (var (single (13)), single (0))
%!assert (var ([1,2,3]), 1)
%!assert (var ([1,2,3], 1), 2/3, eps)
%!assert (var ([1,2,3], [], 1), [0,0,0])
%!assert (var ([1,2,3], [], 3), [0,0,0])
%!assert (var (5, 99), 0)
%!assert (var (5, 99, 1), 0)
%!assert (var (5, 99, 2), 0)
%!assert (var ([5 3], [99 99], 2), 1)
%!assert (var ([1:7], [1:7]), 3)
%!assert (var ([eye(3)], [1:3]), [5/36, 2/9, 1/4], eps)
%!assert (var (ones (2,2,2), [1:2], 3), [(zeros (2,2))])
%!assert (var ([1 2; 3 4], 0, 'all'), var ([1:4]))
%!assert (var (reshape ([1:8], 2, 2, 2), 0, [1 3]), [17/3 17/3], eps)
%!assert (var ([1 2 3;1 2 3], [], [1 2]), 0.8, eps)

## Test empty inputs
%!assert (var ([]), NaN)
%!assert (var ([],[],1), NaN(1,0))
%!assert (var ([],[],2), NaN(0,1))
%!assert (var ([],[],3), [])
%!assert (var (ones (1,0)), NaN)
%!assert (var (ones (1,0), [], 1), NaN(1,0))
%!assert (var (ones (1,0), [], 2), NaN)
%!assert (var (ones (1,0), [], 3), NaN(1,0))
%!assert (var (ones (0,1)), NaN)
%!assert (var (ones (0,1), [], 1), NaN)
%!assert (var (ones (0,1), [], 2), NaN(0,1))
%!assert (var (ones (0,1), [], 3), NaN(0,1))
%!assert (var (ones (1,3,0,2)), NaN(1,1,0,2))
%!assert (var (ones (1,3,0,2), [], 1), NaN(1,3,0,2))
%!assert (var (ones (1,3,0,2), [], 2), NaN(1,1,0,2))
%!assert (var (ones (1,3,0,2), [], 3), NaN(1,3,1,2))
%!assert (var (ones (1,3,0,2), [], 4), NaN(1,3,0))

## Test second output
%!test <*62395>
%! [~, m] = var (13);
%! assert (m, 13);
%! [~, m] = var (single(13));
%! assert (m, single(13));
%! [~, m] = var ([1, 2, 3; 3 2 1], []);
%! assert (m, [2 2 2]);
%! [~, m] = var ([1, 2, 3; 3 2 1], [], 1);
%! assert (m, [2 2 2]);
%! [~, m] = var ([1, 2, 3; 3 2 1], [], 2);
%! assert (m, [2 2]');
%! [~, m] = var ([1, 2, 3; 3 2 1], [], 3);
%! assert (m, [1 2 3; 3 2 1]);

## 2nd output, weighted inputs, vector dims
%!test <*62395>
%! [~, m] = var (5,99);
%! assert (m, 5);
%! [~, m] = var ([1:7], [1:7]);
%! assert (m, 5);
%! [~, m] = var ([eye(3)], [1:3]);
%! assert (m, [1/6, 1/3, 0.5], eps);
%! [~, m] = var (ones (2,2,2), [1:2], 3);
%! assert (m, ones (2,2));
%! [~, m] = var ([1 2; 3 4], 0, 'all');
%! assert (m, 2.5, eps);
%! [~, m] = var (reshape ([1:8], 2, 2, 2), 0, [1 3]);
%! assert (m, [3.5, 5.5], eps);

## 2nd output, empty inputs
%!test <*62395>
%! [~, m] = var ([]);
%! assert (m, NaN);
#%! [~, m] = var ([],[],1);
#%! assert (m, NaN(1,0));
#%! [~, m] = var ([],[],2);
#%! assert (m, NaN(0,1));
#%! [~, m] = var ([],[],3);
#%! assert (m, []);
#%! [~, m] = var (ones (1,3,0,2));
#%! assert (m, NaN(1,1,0,2));

## Test Inf and NaN inputs
%!test <*63203>
%! [v, m] = var (Inf);
%! assert (v, NaN);
%! assert (m, Inf);
%!test <*63203>
%! [v, m] = var (NaN);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([1, Inf, 3]);
%! assert (v, NaN);
%! assert (m, Inf);
%!test <*63203>
%! [v, m] = var ([1, Inf, 3]');
%! assert (v, NaN);
%! assert (m, Inf);
%!test <*63203>
%! [v, m] = var ([1, NaN, 3]);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([1, NaN, 3]');
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([1, Inf, 3], [], 1);
%! assert (v, [0, NaN, 0]);
%! assert (m, [1, Inf, 3]);
%!test <*63203>
%! [v, m] = var ([1, Inf, 3], [], 2);
%! assert (v, NaN);
%! assert (m, Inf);
%!test <*63203>
%! [v, m] = var ([1, Inf, 3], [], 3);
%! assert (v, [0, NaN, 0]);
%! assert (m, [1, Inf, 3]);
%!test <*63203>
%! [v, m] = var ([1, NaN, 3], [], 1);
%! assert (v, [0, NaN, 0]);
%! assert (m, [1, NaN, 3]);
%!test <*63203>
%! [v, m] = var ([1, NaN, 3], [], 2);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([1, NaN, 3], [], 3);
%! assert (v, [0, NaN, 0]);
%! assert (m, [1, NaN, 3]);
%!test <*63203>
%! [v, m] = var ([1, 2, 3; 3, Inf, 5]);
%! assert (v, [2, NaN, 2]);
%! assert (m, [2, Inf, 4]);
%!test <*63203>
%! [v, m] = var ([1, Inf, 3; 3, Inf, 5]);
%! assert (v, [2, NaN, 2]);
%! assert (m, [2, Inf, 4]);
%!test <*63203>
%! [v, m] = var ([1, 2, 3; 3, NaN, 5]);
%! assert (v, [2, NaN, 2]);
%! assert (m, [2, NaN, 4]);
%!test <*63203>
%! [v, m] = var ([1, NaN, 3; 3, NaN, 5]);
%! assert (v, [2, NaN, 2]);
%! assert (m, [2, NaN, 4]);
%!test <*63203>
%! [v, m] = var ([Inf, 2, NaN]);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([Inf, 2, NaN]');
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([NaN, 2, Inf]);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([NaN, 2, Inf]');
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([Inf, 2, NaN], [], 1);
%! assert (v, [NaN, 0, NaN]);
%! assert (m, [Inf, 2, NaN]);
%!test <*63203>
%! [v, m] = var ([Inf, 2, NaN], [], 2);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([NaN, 2, Inf], [], 1);
%! assert (v, [NaN, 0, NaN]);
%! assert (m, [NaN, 2, Inf]);
%!test <*63203>
%! [v, m] = var ([NaN, 2, Inf], [], 2);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = var ([1, 3, NaN; 3, 5, Inf]);
%! assert (v, [2, 2, NaN]);
%! assert (m, [2, 4, NaN]);
%!test <*63203>
%! [v, m] = var ([1, 3, Inf; 3, 5, NaN]);
%! assert (v, [2, 2, NaN]);
%! assert (m, [2, 4, NaN]);

## Test sparse/diagonal inputs
%!test <*63291>
%! [v, m] = var (2 * eye (2));
%! assert (v, [2, 2]);
%! assert (m, [1, 1]);
%!test <*63291>
%! [v, m] = var (4 * eye (2), [1, 3]);
%! assert (v, [3, 3]);
%! assert (m, [1, 3]);
%!test <*63291>
%! [v, m] = var (sparse (2 * eye (2)));
%! assert (full (v), [2, 2]);
%! assert (full (m), [1, 1]);
%!test <*63291>
%! [v, m] = var (sparse (4 * eye (2)), [1, 3]);
%! assert (full (v), [3, 3]);
%! assert (full (m), [1, 3]);

%!test <63291>
%! [v, m] = var (sparse (eye (2)));
%! assert (issparse (v));
%! assert (issparse (m));
%!test <63291>
%! [v, m] = var (sparse (eye (2)), [1, 3]);
%! assert (issparse (v));
%! assert (issparse (m));

## Test input validation
%!error <Invalid call> var ()
%!error <var: X must be a numeric vector or matrix> var (['A'; 'B'])
%!error <var: normalization scalar must be either 0 or 1> var ([1 2 3], 2)
%!error <var: weights must not contain any negative values> var ([1 2], [-1 0])
%!error <var: weight matrix or array does not match X in size> ...
%! var ([1 2], eye (2))
%!error <var: weight matrix or array does not match X in size> ...
%! var (ones (2, 2), [1 2], [1 2])
%!error <var: weight vector does not match first operating dimension> ...
%! var ([1 2], [1 2 3])
%!error <var: weight vector does not match first operating dimension> ...
%! var (1, [1 2])
%!error <var: weight vector does not match given operating dimension> ...
%! var ([1 2], [1 2], 1)
%!error <var: DIM must be a positive integer scalar or vector> ...
%! var (1, [], ones (2,2))
%!error <var: DIM must be a positive integer scalar or vector> var (1, [], 1.5)
%!error <var: DIM must be a positive integer scalar or vector> var (1, [], 0)

%!error <Invalid call to var.  Correct usage is> var ()
%!error <Invalid call to var.  Correct usage is> var (1, 2, "omitnan", 3)
%!error <Invalid call to var.  Correct usage is> var (1, 2, 3, 4)
%!error <Invalid call to var.  Correct usage is> var (1, 2, 3, 4, 5)
%!error <Invalid call to var.  Correct usage is> var (1, "foo")
%!error <Invalid call to var.  Correct usage is> var (1, [], "foo")
%!error <var: normalization scalar must be either 0 or 1> var ([1 2], 2, "all")
%!error <var: normalization scalar must be either 0 or 1> var ([1 2],0.5, "all")
%!error <var: weights must not contain any negative values> var (1, -1)
%!error <var: weights must not contain any negative values> var (1, [1 -1])
%!error <var: weights must not contain any negative values> ...
%! var ([1 2 3], [1 -1 0])
%!error <var: X must be a numeric vector or matrix> var ({1:5})
%!error <var: X must be a numeric vector or matrix> var ("char")
%!error <var: DIM must be a positive integer> var (1, [], ones (2,2))
%!error <var: DIM must be a positive integer> var (1, 0, 1.5)
%!error <var: DIM must be a positive integer> var (1, [], 0)
%!error <var: VECDIM must contain non-repeating positive integers> ...
%! var (repmat ([1:20;6:25], [5 2 6 3]), 0, [1 2 2 2])
%!error <var: weight vector does not match first operating dimension> ...
%! var ([1 2 3; 2 3 4], [1 3 4])
%!error <var: weight matrix or array does not match X in size> ...
%! var ([1 2], eye (2))
%!error <var: weight vector does not match given operating dimension> ...
%! var ([1 2 3; 2 3 4], [1 3 4], 1)
%!error <var: weight vector does not match given operating dimension> ...
%! var ([1 2 3; 2 3 4], [1 3], 2)
%!error <var: weight matrix or array does not match X in size> ...
%! var (repmat ([1:20;6:25], [5 2 6 3]), repmat ([1:20;6:25], [5 2 3]), [2 3])
%!error <var: 'all' flag cannot be used with DIM or VECDIM options> ...
%! var (1, [], 1, "all")
%!error <var: elements in weight vector do not match elements in X> ...
%! var ([1 2 3; 2 3 4], [1 3], "all")
%!error <var: weight matrix or array does not match X in size> ...
%! var (repmat ([1:20;6:25], [5 2 6 3]), repmat ([1:20;6:25], [5 2 3]), "all")

