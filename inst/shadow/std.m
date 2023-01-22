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
## If @var{x} is a vector, then @code{std (@var{x})} returns the standard
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
## If @var{x} is a matrix, then @code{std (@var{x})} returns a row vector with
## the standard deviation of each column in @var{x}.
##
## @item
## If @var{x} is a multi-dimensional array, then @code{std (@var{x})} operates
## along the first non-singleton dimension of @var{x}.
## @end itemize
##
## @code{std (@var{x}, @var{w})} specifies a weighting scheme.  When @var{w} = 0
## (default), the standard deviation is normalized by N-1 (population standard
## deviation), where N is the number of observations.  When @var{w} = 1, the
## standard deviation is normalized by the number of observations (sample
## standard deviation).  To use the default value you may pass an empty input
## argument [] before entering other options.
##
## @var{w} can also be an array of non-negative numbers.  When @var{w} is a
## vector, it must have the same length as the number of elements in the
## operating dimension of @var{x}.  If @var{w} is a matrix or n-D array, or the
## operating dimension is supplied as a @var{vecdim} or "all", @var{w} must be
## the same size as @var{x}.  NaN values are permitted in @var{w}, will be
## multiplied with the associated values in @var{x}, and can be excluded by the
## @var{nanflag} option.
##
## @code{std (@var{x}, [], @var{dim})} returns the standard deviation along the
## operating dimension @var{dim} of @var{x}.  For @var{dim} greater than
## @code{ndims (@var{x})}, then @var{s} is returned as zeros of the same size as
## @var{x} and @var{m} = @var{x}.
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
## from the calculation using any of the input argument combinations in previous
## syntaxes.  The default value for @var{nanflag} is "includenan", and keeps NaN
## values in the calculation. To exclude NaN values, set the value of
## @var{nanflag} to "omitnan".
##
## @code{[@var{s}, @var{m}] = std (@dots{})} also returns the mean of the
## elements of @var{x} used to calculate the standard deviation.  If @var{s} is
## the weighted standard deviation, then @var{m} is the weighted mean.
##
## @seealso{var, mean}
## @end deftypefn

function [s, m] = std (x, varargin)

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
      switch (lower (i{:}))
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
  weighted = false; # true if weight vector/array used
  vecdim = [];
  vecempty = true;
  vecdim_scalar_vector = [false, false]; # [false, false] for empty vecdim
  szx = size (x);
  ndx = ndims (x);

  ## Check numeric arguments
  if (! (isnumeric (x)))
    error ("std: X must be a numeric vector or matrix");
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
      error ("std: weights must not contain any negative values");
    endif

    if (isscalar (varargin{1}))
      w = varargin{1};
      if (! (w == 0 || w == 1) && ! isscalar (x))
        error ("std: normalization scalar must be either 0 or 1");
      endif
    elseif (numel (varargin{1}) > 1)
      weights = varargin{1};
      weighted = true;
    endif
    if (nvarg > 1)
      ## Process dimension input
      vecdim = varargin{2};
      if (! (vecempty = isempty (vecdim)))
        vecdim_scalar_vector = [isscalar(vecdim), isvector(vecdim)];
      endif
      if (! (vecdim_scalar_vector(2) && all (vecdim)) || any (rem (vecdim, 1)))
        error ("std: DIM must be a positive integer scalar or vector");
      endif
      if (numel (vecdim) != numel (unique (vecdim)))
        error ("std: VECDIM must contain non-repeating positive integers");
      endif
      if (! isempty (x) && vecdim_scalar_vector(1) && vecdim > ndx)
        s = zeros (szx, outtype);
        sn = ! isfinite (x);
        s(sn) = NaN;
        m = x;
        return;
      endif
    endif
  endif

  ## Check for conflicting input arguments
  if (all_flag && ! vecempty)
    error ("std: 'all' flag cannot be used with DIM or VECDIM options");
  endif
  if (weighted)
    if (all_flag)
      if (isvector (weights))
        if (numel (weights) != numel (x))
          error ("std: elements in weight vector do not match elements in X");
        endif
      elseif (! (isequal (size (weights), szx)))
        error ("std: weight matrix or array does not match X in size");
      endif

    elseif (vecempty)
      dim = find (szx > 1, 1);
      if length (dim) == 0
        dim = 1;
      endif
      if (numel (weights) != szx(dim))
        if (isvector (weights))
          error ("std: weight vector does not match first operating dimension");
        elseif (! isequal (size (weights), szx))
          error ("std: weight matrix or array does not match X in size");
        endif
      endif
    elseif (vecdim_scalar_vector(1))
      if (isvector (weights) && numel (weights) != szx(vecdim))
        error ("std: weight vector does not match given operating dimension");
      endif
    elseif (vecdim_scalar_vector(2) && ! (isequal (size (weights), szx)))
      error ("std: weight matrix or array does not match X in size");
    endif
  endif

  ## Force output for X being empty or scalar
  if (isempty (x))
    if (vecempty && (ndx == 2 || all ((szx) == 0)))
      s = NaN;
      m = NaN;
      return;
    endif
    if (vecdim_scalar_vector(1))
      szx(vecdim) = 1;
      s = NaN (szx);
      m = NaN (szx);
      return;
    endif
  endif
  if (isscalar (x))
    if (isfinite (x))
      s = zeros (outtype);
    else
      s = NaN (outtype);
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
      s = sqrt (sum (abs (x - m) .^ 2) ./ (n - 1 + w));
      if (n == 1)
        s = 0;
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
      dims = ones (1, ndx);
      dims(dim) = szx(dim);
      m_exp = repmat (m, dims);
      if (omitnan)
        x(xn) = m_exp(xn);
      endif
      s = sqrt (sumsq (x - m_exp, dim) ./ (n - 1 + w));
      if (numel (n) == 1)
        divby0 = repmat (n, size (s)) == 1;
      else
         divby0 = n == 1;
      endif
      s(divby0) = 0;
    endif

  elseif (nvarg == 1)
    ## Two numeric input arguments, w or weights given.
    if (all_flag)
      x = x(:);
      if (weighted)
        wv = weights(:);
        wx = wv .* x;
      else
        wv = ones (length (x), 1);
        wx = x;
      endif

      if (omitnan)
        xn = isnan (wx);
        wx = wx(! xn);
        wv = wv(! xn);
        x = x(! xn);
      endif
      n = length (wx);
      m = sum (wx) ./ sum (wv);
      if (weighted)
        s = sqrt (sum (wv .* (abs (x - m) .^ 2)) ./ sum (wv));
      else
        s = sqrt (sum (wv .* (abs (x - m) .^ 2)) ./ (n - 1 + w));
        if (n == 1)
          s = 0;
        endif
      endif
    else
      dim = find (szx > 1, 1);
      if length (dim) == 0
        dim = 1;
      endif
      if (! weighted)
        wv = ones (szx);
        wx = x;
      else
        if (isvector (weights))
          wv = zeros (szx) + shiftdim (weights(:), 1 - dim);
        else
          wv = weights;
        endif
        wx = wv .* x;
      endif
      n = size (wx, dim);
      if (omitnan)
        xn = isnan (wx);
        n = sum (! xn, dim);
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
        s = sqrt (sum (wv .* ((x - m_exp) .^ 2), dim) ./ sum (wv, dim));
      else
        s = sqrt (sumsq (x - m_exp, dim) ./ (n - 1 + w));
        if (numel (n) == 1)
          divby0 = repmat (n, size (s)) == 1;
        else
           divby0 = n == 1;
        endif
        s(divby0) = 0;
      endif
    endif

  elseif (nvarg == 2)
    ## Three numeric input arguments, both w or weights and dim or vecdim given.
    if (vecdim_scalar_vector(1))
      if (!weighted)
        wv = ones (szx);
        wx = x;
      else
        if (isvector (weights))
          wv = zeros (szx) + shiftdim (weights(:), 1 - vecdim);
        else
          wv = weights;
        endif
        wx = wv .* x;
      endif
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
        s = sqrt (sum (wv .* ((x - m_exp) .^ 2), vecdim) ./ sum (wv, vecdim));
      else
        s = sumsq (x - m_exp, vecdim);
        sn = isnan (s);
        s = sqrt (s ./ (n - 1 + w));
        if (numel (n) == 1)
          divby0 = repmat (n, size (s)) == 1;
        else
          divby0 = n == 1;
        endif
        s(divby0) = 0;
        s(sn) = NaN;
      endif

    else
    ## Weights and nonscalar vecdim specified

    ## Ignore exceeding dimensions in VECDIM
      remdims = 1:ndx;    # all dimensions
      vecdim(find (vecdim > ndx)) = [];
      ## Calculate permutation vector
      remdims(vecdim) = [];     # delete dimensions specified by vecdim
      nremd = numel (remdims);

      ## If all dimensions are given, it is similar to all flag
      if (nremd == 0)
        x = x(:);
        if (weighted)
          weights = weights(:);
          wx = weights .* x;
        else
          weights = ones (length (x), 1);
          wx = x;
        endif
        if (omitnan)
          xn = isnan (wx);
          wx = wx(! xn);
          weights = weights(! xn);
          x = x(! xn);
        endif
        n = length (wx);
        m = sum (wx) ./ sum (weights);
        if (weighted)
          s = sqrt (sum (weights .* (abs (x - m) .^ 2)) ./ sum (weights));
        else
          s = sqrt (sum (weights .* (abs (x - m) .^ 2)) ./ (n - 1 + w));
          if (n == 1)
            s = 0;
          endif
        endif

      else
        ## Apply weights
        if (weighted)
          wv = weights;
          wx = wv .* x;
        else
          wv = ones (szx);
          wx = x;
        endif

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
          xn = isnan (wx);
          n = sum (! xn, dim);
          wx(xn) = 0;
          wv(xn) = 0;
        endif
        m = sum (wx, dim) ./ sum (wv, dim);
        m_exp = zeros (size (wx)) + shiftdim (m, 0);
        if (omitnan)
          x(xn) = m_exp(xn);
        endif
        if (weighted)
          s = sqrt (sum (wv .* ((x - m_exp) .^ 2), dim) ./ sum (wv, dim));
        else
          s = sqrt (sumsq (x - m_exp, dim) ./ (n - 1 + w));
          if (numel (n) == 1)
            divby0 = repmat (n, size (s)) == 1;
          else
             divby0 = n == 1;
          endif
          s(divby0) = 0;
        endif

        ## Inverse permute back to correct dimensions
        s = ipermute (s, perm);
        m = ipermute (m, perm);
      endif
    endif
  endif

  ## Preserve class type
  switch outtype
    case 'single'
      s = single (s);
      m = single (m);
    case 'double'
      s = double (s);
      m = double (m);
  endswitch

endfunction


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

## Tests for different weight and omitnan code paths
%!assert (std ([4 NaN 6], [1 2 1], "omitnan"), 1, eps)
%!assert (std ([4 5 6], [1 NaN 1], "omitnan"), 1, eps)
%!assert (std (magic(3), [1 NaN 3], "omitnan"), sqrt(3)*[1 2 1], eps)
%!assert (std ([4 NaN 6], [1 2 1], "omitnan", "all"), 1, eps)
%!assert (std ([4 NaN 6], [1 2 1], "all", "omitnan"), 1, eps)
%!assert (std ([4 5 6], [1 NaN 1], "omitnan", "all"), 1, eps)
%!assert (std ([4 NaN 6], [1 2 1], 2, "omitnan"), 1, eps)
%!assert (std ([4 5 6], [1 NaN 1], 2, "omitnan"), 1, eps)
%!assert (std (magic(3), [1 NaN 3], 1, "omitnan"), sqrt(3)*[1 2 1], eps)
%!assert (std (magic(3), [1 NaN 3], 2, "omitnan"), sqrt(3)*[0.5;1;0.5], eps)
%!assert (std (4*[4 5; 6 7; 8 9], [1 3], 2, 'omitnan'), sqrt(3)*[1;1;1], eps)
%!assert (std ([4 NaN; 6 7; 8 9], [1 1 3], 1, 'omitnan'), [1.6 sqrt(3)/2], eps)
%!assert (std (4*[4 NaN; 6 7; 8 9], [1 3], 2, 'omitnan'), sqrt(3)*[0;1;1], eps)
%!assert (std (3*reshape(1:18, [3 3 2]), [1 2 3], 1, 'omitnan'), sqrt(5)*ones(1,3,2), eps)
%!assert (std (reshape(1:18, [3 3 2]), [1 2 3], 2, 'omitnan'), sqrt(5)*ones(3,1,2), eps)
%!assert (std (3*reshape(1:18, [3 3 2]), ones (3,3,2), [1 2], 'omitnan'), sqrt(60)*ones(1,1,2),eps)
%!assert (std (3*reshape(1:18, [3 3 2]), ones (3,3,2), [1 4], 'omitnan'), sqrt(6)*ones(1,3,2),eps)
%!assert (std (6*reshape(1:18, [3 3 2]), ones (3,3,2), [1:3], 'omitnan'), sqrt(969),eps)
%!test
%! x = reshape(1:18, [3 3 2]);
%! x([2, 14]) = NaN;
%! w = ones (3,3,2);
%! assert (std (16*x, w, [1:3], 'omitnan'), sqrt(6519), eps);
%!test
%! x = reshape(1:18, [3 3 2]);
%! w = ones (3,3,2);
%! w([2, 14]) = NaN;
%! assert (std (16*x, w, [1:3], 'omitnan'), sqrt(6519), eps);

## Test input case insensitivity
%!assert (std ([1 2 3], "aLl"), 1);
%!assert (std ([1 2 3], "OmitNan"), 1);
%!assert (std ([1 2 3], "IncludeNan"), 1);

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5, 2, 6, 3]);
%! assert (size (std (x, 0, [3 2])), [10, 1, 1, 3]);
%! assert (size (std (x, 1, [1 2])), [1, 1, 6, 3]);
%! assert (size (std (x, [], [1 2 4])), [1, 1, 6]);
%! assert (size (std (x, 0, [1 4 3])), [1, 40]);
%! assert (size (std (x, [], [1 2 3 4])), [1, 1]);

## Test matrix with vecdim, weighted, matrix weights, omitnan
%!assert (std (3*magic(3)), sqrt([63 144 63]), eps)
%!assert (std (3*magic(3), 'omitnan'), sqrt([63 144 63]), eps)
%!assert (std (3*magic(3), 1), sqrt([42 96 42]), eps)
%!assert (std (3*magic(3), 1, 'omitnan'), sqrt([42 96 42]), eps)
%!assert (std (3*magic(3), ones(1,3), 1), sqrt([42 96 42]), eps)
%!assert (std (3*magic(3), ones(1,3), 1, 'omitnan'), sqrt([42 96 42]), eps)
%!assert (std (2*magic(3), [1 1 NaN], 1, 'omitnan'), [5 4 1], eps)
%!assert (std (3*magic(3), ones(3,3)), sqrt([42 96 42]), eps)
%!assert (std (3*magic(3), ones(3,3), 'omitnan'), sqrt([42 96 42]), eps)
%!assert (std (3*magic(3), [1 1 1; 1 1 1; 1 NaN 1], 'omitnan'), sqrt([42 36 42]), eps)
%!assert (std (3*magic(3), ones(3,3), 1), sqrt([42 96 42]), eps)
%!assert (std (3*magic(3), ones(3,3), 1, 'omitnan'), sqrt([42 96 42]), eps)
%!assert (std (3*magic(3), [1 1 1; 1 1 1; 1 NaN 1], 1, 'omitnan'), sqrt([42 36 42]), eps)
%!assert (std (3*magic(3), ones(3,3), [1 4]), sqrt([42 96 42]), eps)
%!assert (std (3*magic(3), ones(3,3), [1 4], 'omitnan'), sqrt([42 96 42]), eps)
%!assert (std (3*magic(3), [1 1 1; 1 1 1; 1 NaN 1],[1 4],'omitnan'), sqrt([42 36 42]), eps)

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

####
#### BISTs from core Octave
####

%!assert (std (13), 0)
%!assert (std (single (13)), single (0))
%!assert (std ([1,2,3]), 1)
%!assert (std ([1,2,3], 1), sqrt (2/3), eps)
%!assert (std ([1,2,3], [], 1), [0,0,0])
%!assert (std ([1,2,3], [], 3), [0,0,0])
%!assert (std (5, 99), 0)
%!assert (std (5, 99, 1), 0)
%!assert (std (5, 99, 2), 0)
%!assert (std ([5 3], [99 99], 2), 1)
%!assert (std ([1:7], [1:7]), sqrt (3))
%!assert (std ([eye(3)], [1:3]), sqrt ([5/36, 2/9, 1/4]), eps)
%!assert (std (ones (2,2,2), [1:2], 3), [(zeros (2,2))])
%!assert (std ([1 2; 3 4], 0, 'all'), std ([1:4]))
%!assert (std (reshape ([1:8], 2, 2, 2), 0, [1 3]), sqrt ([17/3 17/3]), eps)
%!assert (std ([1 2 3;1 2 3], [], [1 2]), sqrt (0.8), eps)

## Test empty inputs
%!assert (std ([]), NaN)
%!assert (std ([],[],1), NaN(1,0))
%!assert (std ([],[],2), NaN(0,1))
%!assert (std ([],[],3), [])
%!assert (std (ones (1,0)), NaN)
%!assert (std (ones (1,0), [], 1), NaN(1,0))
%!assert (std (ones (1,0), [], 2), NaN)
%!assert (std (ones (1,0), [], 3), NaN(1,0))
%!assert (std (ones (0,1)), NaN)
%!assert (std (ones (0,1), [], 1), NaN)
%!assert (std (ones (0,1), [], 2), NaN(0,1))
%!assert (std (ones (0,1), [], 3), NaN(0,1))
%!assert (std (ones (1,3,0,2)), NaN(1,1,0,2))
%!assert (std (ones (1,3,0,2), [], 1), NaN(1,3,0,2))
%!assert (std (ones (1,3,0,2), [], 2), NaN(1,1,0,2))
%!assert (std (ones (1,3,0,2), [], 3), NaN(1,3,1,2))
%!assert (std (ones (1,3,0,2), [], 4), NaN(1,3,0))

## Test second output
%!test <*62395>
%! [~, m] = std (13);
%! assert (m, 13);
%! [~, m] = std (single(13));
%! assert (m, single(13));
%! [~, m] = std ([1, 2, 3; 3 2 1], []);
%! assert (m, [2 2 2]);
%! [~, m] = std ([1, 2, 3; 3 2 1], [], 1);
%! assert (m, [2 2 2]);
%! [~, m] = std ([1, 2, 3; 3 2 1], [], 2);
%! assert (m, [2 2]');
%! [~, m] = std ([1, 2, 3; 3 2 1], [], 3);
%! assert (m, [1 2 3; 3 2 1]);

## 2nd output, weighted inputs, vector dims
%!test <*62395>
%! [~, m] = std (5,99);
%! assert (m, 5);
%! [~, m] = std ([1:7], [1:7]);
%! assert (m, 5);
%! [~, m] = std ([eye(3)], [1:3]);
%! assert (m, [1/6, 1/3, 0.5], eps);
%! [~, m] = std (ones (2,2,2), [1:2], 3);
%! assert (m, ones (2,2));
%! [~, m] = std ([1 2; 3 4], 0, 'all');
%! assert (m, 2.5, eps);
%! [~, m] = std (reshape ([1:8], 2, 2, 2), 0, [1 3]);
%! assert (m, [3.5, 5.5], eps);

## 2nd output, empty inputs
%!test <*62395>
%! [~, m] = std ([]);
%! assert (m, NaN);
#%! [~, m] = std ([],[],1);
#%! assert (m, NaN(1,0));
#%! [~, m] = std ([],[],2);
#%! assert (m, NaN(0,1));
#%! [~, m] = std ([],[],3);
#%! assert (m, []);
#%! [~, m] = std (ones (1,3,0,2));
#%! assert (m, NaN(1,1,0,2));

## Test Inf and NaN inputs
%!test <*63203>
%! [v, m] = std (Inf);
%! assert (v, NaN);
%! assert (m, Inf);
%!test <*63203>
%! [v, m] = std (NaN);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([1, Inf, 3]);
%! assert (v, NaN);
%! assert (m, Inf);
%!test <*63203>
%! [v, m] = std ([1, Inf, 3]');
%! assert (v, NaN);
%! assert (m, Inf);
%!test <*63203>
%! [v, m] = std ([1, NaN, 3]);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([1, NaN, 3]');
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([1, Inf, 3], [], 1);
%! assert (v, [0, NaN, 0]);
%! assert (m, [1, Inf, 3]);
%!test <*63203>
%! [v, m] = std ([1, Inf, 3], [], 2);
%! assert (v, NaN);
%! assert (m, Inf);
%!test <*63203>
%! [v, m] = std ([1, Inf, 3], [], 3);
%! assert (v, [0, NaN, 0]);
%! assert (m, [1, Inf, 3]);
%!test <*63203>
%! [v, m] = std ([1, NaN, 3], [], 1);
%! assert (v, [0, NaN, 0]);
%! assert (m, [1, NaN, 3]);
%!test <*63203>
%! [v, m] = std ([1, NaN, 3], [], 2);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([1, NaN, 3], [], 3);
%! assert (v, [0, NaN, 0]);
%! assert (m, [1, NaN, 3]);
%!test <*63203>
%! [v, m] = std ([1, 2, 3; 3, Inf, 5]);
%! assert (v, sqrt ([2, NaN, 2]));
%! assert (m, [2, Inf, 4]);
%!test <*63203>
%! [v, m] = std ([1, Inf, 3; 3, Inf, 5]);
%! assert (v, sqrt ([2, NaN, 2]));
%! assert (m, [2, Inf, 4]);
%!test <*63203>
%! [v, m] = std ([1, 2, 3; 3, NaN, 5]);
%! assert (v, sqrt ([2, NaN, 2]));
%! assert (m, [2, NaN, 4]);
%!test <*63203>
%! [v, m] = std ([1, NaN, 3; 3, NaN, 5]);
%! assert (v, sqrt ([2, NaN, 2]));
%! assert (m, [2, NaN, 4]);
%!test <*63203>
%! [v, m] = std ([Inf, 2, NaN]);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([Inf, 2, NaN]');
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([NaN, 2, Inf]);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([NaN, 2, Inf]');
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([Inf, 2, NaN], [], 1);
%! assert (v, [NaN, 0, NaN]);
%! assert (m, [Inf, 2, NaN]);
%!test <*63203>
%! [v, m] = std ([Inf, 2, NaN], [], 2);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([NaN, 2, Inf], [], 1);
%! assert (v, [NaN, 0, NaN]);
%! assert (m, [NaN, 2, Inf]);
%!test <*63203>
%! [v, m] = std ([NaN, 2, Inf], [], 2);
%! assert (v, NaN);
%! assert (m, NaN);
%!test <*63203>
%! [v, m] = std ([1, 3, NaN; 3, 5, Inf]);
%! assert (v, sqrt ([2, 2, NaN]));
%! assert (m, [2, 4, NaN]);
%!test <*63203>
%! [v, m] = std ([1, 3, Inf; 3, 5, NaN]);
%! assert (v, sqrt ([2, 2, NaN]));
%! assert (m, [2, 4, NaN]);

## Test sparse/diagonal inputs
%!test <*63291>
%! [v, m] = std (2 * eye (2));
%! assert (v, sqrt ([2, 2]));
%! assert (m, [1, 1]);
%!test <*63291>
%! [v, m] = std (4 * eye (2), [1, 3]);
%! assert (v, sqrt ([3, 3]));
%! assert (m, [1, 3]);
%!test <*63291>
%! [v, m] = std (sparse (2 * eye (2)));
%! assert (full (v), sqrt ([2, 2]));
%! assert (full (m), [1, 1]);
%!test <*63291>
%! [v, m] = std (sparse (4 * eye (2)), [1, 3]);
%! assert (full (v), sqrt ([3, 3]));
%! assert (full (m), [1, 3]);

%!test <63291>
%! [v, m] = std (sparse (eye (2)));
%! assert (issparse (v));
%! assert (issparse (m));
%!test <63291>
%! [v, m] = std (sparse (eye (2)), [1, 3]);
%! assert (issparse (v));
%! assert (issparse (m));

## Test input validation
%!error <Invalid call> std ()
%!error <std: X must be a numeric vector or matrix> std (['A'; 'B'])
%!error <std: normalization scalar must be either 0 or 1> std ([1 2 3], 2)
%!error <std: weights must not contain any negative values> std ([1 2], [-1 0])
%!error <std: weight matrix or array does not match X in size> ...
%! std ([1 2], eye (2))
%!error <std: weight matrix or array does not match X in size> ...
%! std (ones (2, 2), [1 2], [1 2])
%!error <std: weight vector does not match first operating dimension> ...
%! std ([1 2], [1 2 3])
%!error <std: weight vector does not match first operating dimension> ...
%! std (1, [1 2])
%!error <std: weight vector does not match given operating dimension> ...
%! std ([1 2], [1 2], 1)
%!error <std: DIM must be a positive integer scalar or vector> ...
%! std (1, [], ones (2,2))
%!error <std: DIM must be a positive integer scalar or vector> std (1, [], 1.5)
%!error <std: DIM must be a positive integer scalar or vector> std (1, [], 0)

%!error <Invalid call to std.  Correct usage is> std ()
%!error <Invalid call to std.  Correct usage is> std (1, 2, "omitnan", 3)
%!error <Invalid call to std.  Correct usage is> std (1, 2, 3, 4)
%!error <Invalid call to std.  Correct usage is> std (1, 2, 3, 4, 5)
%!error <Invalid call to std.  Correct usage is> std (1, "foo")
%!error <Invalid call to std.  Correct usage is> std (1, [], "foo")
%!error <std: normalization scalar must be either 0 or 1> std ([1 2], 2, "all")
%!error <std: normalization scalar must be either 0 or 1> std ([1 2],0.5, "all")
%!error <std: weights must not contain any negative values> std (1, -1)
%!error <std: weights must not contain any negative values> std (1, [1 -1])
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
%!error <std: weight matrix or array does not match X in size> ...
%! std ([1 2], eye (2))
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
