## Copyright (C) 2001 Paul Kienzle <pkienzle@users.sf.net>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{m} =} trimmean (@var{x}, @var{p})
## @deftypefnx {statistics} {@var{m} =} trimmean (@var{x}, @var{p}, @var{flag})
## @deftypefnx {statistics} {@var{m} =} trimmean (@dots{}, @qcode{"all"})
## @deftypefnx {statistics} {@var{m} =} trimmean (@dots{}, @var{dim})
## @deftypefnx {statistics} {@var{m} =} trimmean (@dots{}, @var{vecdim})
##
## Compute the trimmed mean.
##
## The trimmed mean of @var{x} is defined as the mean of @var{x} excluding the
## highest and lowest @math{k} data values of @var{x}, calculated as
## @qcode{@var{k} = n * (@var{p} / 100) / 2)}, where @var{n} is the sample size.
##
## @code{@var{m} = trimmean (@var{x}, @var{p})} returns the mean of @var{x}
## after removing the outliers in @var{x} defined by @var{p} percent.
## @itemize
## @item If @var{x} is a vector, then @code{trimmean (@var{x}, @var{p})} is the
## mean of all the values of @var{x}, computed after removing the outliers.
## @item If @var{x} is a matrix, then @code{trimmean (@var{x}, @var{p})} is a
## row vector of column means, computed after removing the outliers.
## @item If @var{x} is a multidimensional array, then @code{trimmean} operates
## along the first nonsingleton dimension of @var{x}.
## @end itemize
##
## To specify the operating dimension(s) when @var{x} is a matrix or a
## multidimensional array, use the @var{dim} or @var{vecdim} input argument.
##
## @code{trimmean} treats @qcode{NaN} values in @var{x} as missing values and
## removes them.
##
## @code{@var{m} = trimmean (@var{x}, @var{p}, @var{flag})} specifies how to
## trim when @math{k}, i.e. half the number of outliers, is not an integer.
## @var{flag} can be specified as one of the following values:
## @multitable @columnfractions 0.2 0.05 0.75
## @headitem Value @tab @tab Description
## @item @qcode{"round"} @tab @tab Round @math{k} to the nearest integer.  This
## is the default.
## @item @qcode{"floor"} @tab @tab Round @math{k} down to the next smaller
## integer.
## @item @qcode{"weighted"} @tab @tab If @math{k = i + f}, where @math{i} is an
## integer and @math{f} is a fraction, compute a weighted mean with weight
## @math{(1 - f)} for the @math{(i + 1)}-th and @math{(n - i)}-th values, and
## full weight for the values between them.
## @end multitable
##
## @code{@var{m} = trimmean (@dots{}, @qcode{"all"})} returns the trimmed mean
## of all the values in @var{x} using any of the input argument combinations in
## the previous syntaxes.
##
## @code{@var{m} = trimmean (@dots{}, @var{dim})} returns the trimmed mean along
## the operating dimension @var{dim} specified as a positive integer scalar.  If
## not specified, then the default value is the first nonsingleton dimension of
## @var{x}, i.e. whose size does not equal 1.  If @var{dim} is greater than
## @qcode{ndims (@var{X})} or if @qcode{size (@var{x}, @var{dim})} is 1, then
## @code{trimmean} returns @var{x}.
##
## @code{@var{m} = trimmean (@dots{}, @var{vecdim})} returns the trimmed mean
## over the dimensions specified in the vector @var{vecdim}.  For example, if
## @var{x} is a 2-by-3-by-4 array, then @code{mean (@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the mean of the
## elements on the corresponding page of @var{x}.  If @var{vecdim} indexes all
## dimensions of @var{x}, then it is equivalent to @code{mean (@var{x}, "all")}.
## Any dimension in @var{vecdim} greater than @code{ndims (@var{x})} is ignored.
##
## @seealso{mean}
## @end deftypefn

function m = trimmean (x, p, varargin)

  if (nargin < 2 || nargin > 4)
    print_usage;
  endif

  if (p < 0 || p >= 100)
    error ("trimmean: invalid percent.");
  endif

  ## Parse extra arguments
  if (nargin < 3)
    flag = [];
    dim = [];
  elseif (nargin < 4)
    if (ischar (varargin{1}) && ! strcmpi (varargin{1}, "all"))
      flag = varargin{1};
      dim = [];
    elseif (isnumeric (varargin{1}) || strcmpi (varargin{1}, "all"))
      flag = [];
      dim = varargin{1};
    endif
  else
    flag = varargin{1};
    dim = varargin{2};
  endif

  ## Get size of X
  szx = size (x);
  ndx = numel (szx);

  ## Handle special case X = []
  if (isempty (x))
    m = NaN (class (x));
    return
  endif

  ## Check FLAG
  if (isempty (flag))
    flag = "round";
  endif
  if (! any (strcmpi (flag, {"round", "floor", "weighted"})))
    error ("trimmean: invalid FLAG argument.");
  endif

  ## Check DIM
  if (isempty (dim))
    (dim = find (szx != 1, 1)) || (dim = 1);
  endif
  if (strcmpi (dim, "all"))
    x = x(:);
    dim = 1;
  endif
  if (! (isvector (dim) && all (dim > 0) && all (rem (dim, 1) == 0)))
    error ("trimmean: DIM must be a positive integer scalar or vector.");
  endif
  vecdim_flag = false;
  if (numel (dim) > 1)
    dim = sort (dim);
    if (! all (diff (dim)))
      error ("trimmean: VECDIM must contain non-repeating positive integers.");
    endif
    vecdim_flag = true;
  endif

  ## If DIM is a scalar greater than ndims, return X
  if (isscalar (dim) && dim > ndx)
    m = x;
    return
  endif

  ## If DIM is a scalar and size (x, dim) == 1, return X
  if (isscalar (dim) && size (x, dim) == 1)
    m = x;
    return
  endif

  ## If DIM is a vector, ignore any value > ndims (x)
  if (numel (dim) > 1)
    dim(dim > ndx) = [];
  endif

  ## Permute dim to simplify all operations along dim1.  At func. end ipermute.
  if (numel (dim) > 1 || (dim != 1 && ! isvector (x)))
    perm = 1:ndx;
    if (! vecdim_flag)
      ## Move dim to dim 1
      perm([1, dim]) = [dim, 1];
      x = permute (x, perm);
      szx([1, dim]) = szx([dim, 1]);
      dim = 1;
    else
      ## Move vecdims to front
      perm(dim) = [];
      perm = [dim, perm];
      x = permute (x, perm);
      ## Reshape all vecdims into dim1
      num_dim = prod (szx(dim));
      szx(dim) = [];
      szx = [num_dim, ones(1, numel(dim)-1), szx];
      x = reshape (x, szx);
      dim = 1;
    endif
    perm_flag = true;
  else
    perm_flag = false;
  endif

  ## Create output matrix
  sizem = size (x);
  sizem(dim) = 1;

  ## Sort X along 1st dimensions
  x = sort (x);

  ## No missing data, all columns have the same length
  if (! any (isnan (x(:))))
    if (isempty (x))
        n = 0;
    else
        n = size (x, 1);
    end
    m = trim (x, n, p, flag, sizem);
    m = reshape (m, sizem);
  ## With missing data, each column is computed separately
  else
    m = NaN (sizem, class (x));
    for j = 1:prod (sizem(2:end))
      n = find (! isnan (x(:,j)), 1, "last");
      m(j) = trim (x(:,j), n, p, flag, [1, 1]);
    endfor
  endif

  ## Inverse permute back to correct dimensions (if necessary)
  if (perm_flag)
    m = ipermute (m, perm);
  endif

endfunction

## Help function for handling different flags
function m = trim (x, n, p, flag, sizem)
  switch (lower (flag))
    case "round"
      k = n * p / 200;
      k0 = round (k - eps (k));
      if (! isempty (n) && n > 0 && k0 < n / 2)
        m = mean (x((k0+1):(n-k0),:), 1);
      else
        m = NaN (sizem, class (x));
      endif
    case "floor"
      k0 = floor (n * p / 200);
      if (! isempty (n) && n > 0 && k0 < n / 2)
        m = mean (x((k0+1):(n-k0),:), 1);
      else
        m = NaN (sizem, class (x));
      endif
    case "weighted"
      k = n * p / 200;
      k0 = floor (k);
      fr = 1 + k0 - k;
      if (! isempty (n) && n > 0 && (k0 < n / 2 || fr > 0))
        m = (sum (x((k0+2):(n-k0-1),:),1) + fr * x(k0+1,:) + fr * x(n-k0,:)) ...
             / (max (0, n - 2 * k0 - 2) + 2 * fr);
      else
        m = NaN (sizem, class (x));
      endif
  endswitch
endfunction

## Test output
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! assert (trimmean (x, 10, "all"), 19.4722, 1e-4);
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! out = trimmean (x, 10, [1, 2]);
%! assert (out(1,1,1), 10.3889, 1e-4);
%! assert (out(1,1,2), 29.6111, 1e-4);
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! x([4, 38]) = NaN;
%! assert (trimmean (x, 10, "all"), 19.3824, 1e-4);
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! out = trimmean (x, 10, 1);
%! assert (out(:,:,1), [-17.6, 8, 13, 18]);
%! assert (out(:,:,2), [23, 28, 33, 10.6]);
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! x([4, 38]) = NaN;
%! out = trimmean (x, 10, 1);
%! assert (out(:,:,1), [-23, 8, 13, 18]);
%! assert (out(:,:,2), [23, 28, 33, 3.75]);
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! out = trimmean (x, 10, 2);
%! assert (out(:,:,1), [8.5; 9.5; -15.25; 11.5; 12.5]);
%! assert (out(:,:,2), [28.5; -4.75; 30.5; 31.5; 32.5]);
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! x([4, 38]) = NaN;
%! out = trimmean (x, 10, 2);
%! assert (out(:,:,1), [8.5; 9.5; -15.25; 14; 12.5]);
%! assert (out(:,:,2), [28.5; -4.75; 28; 31.5; 32.5]);
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! out = trimmean (x, 10, [1, 2, 3]);
%! assert (out, trimmean (x, 10, "all"));

## Test N-D array with NaNs
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! x([4, 38]) = NaN;
%! out = trimmean (x, 10, [1, 2]);
%! assert (out(1,1,1), 10.7647, 1e-4);
%! assert (out(1,1,2), 29.1176, 1e-4);
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! x([4, 38]) = NaN;
%! out = trimmean (x, 10, [1, 3]);
%! assert (out, [2.5556, 18, 23, 11.6667], 1e-4);
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! x([4, 38]) = NaN;
%! out = trimmean (x, 10, [2, 3]);
%! assert (out, [18.5; 2.3750; 3.2857; 24; 22.5], 1e-4);

%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! x([4, 38]) = NaN;
%! out = trimmean (x, 10, [1, 2, 3]);
%! assert (out, trimmean (x, 10, "all"));
%!test
%! x = reshape (1:40, [5, 4, 2]);
%! x([3, 37]) = -100;
%! x([4, 38]) = NaN;
%! out = trimmean (x, 10, [2, 3, 5]);
%! assert (out, [18.5; 2.3750; 3.2857; 24; 22.5], 1e-4);

## Test special cases
%!assert (trimmean (reshape (1:40, [5, 4, 2]), 10, 4), reshape(1:40, [5, 4, 2]))
%!assert (trimmean ([], 10), NaN)
%!assert (trimmean ([1;2;3;4;5], 10, 2), [1;2;3;4;5])

## Test input validation
%!error<Invalid call to trimmean.  Correct usage is:> trimmean (1)
%!error<Invalid call to trimmean.  Correct usage is:> trimmean (1,2,3,4,5)
%!error<trimmean: invalid percent.> trimmean ([1 2 3 4], -10)
%!error<trimmean: invalid percent.> trimmean ([1 2 3 4], 100)
%!error<trimmean: invalid FLAG argument.> trimmean ([1 2 3 4], 10, "flag")
%!error<trimmean: invalid FLAG argument.> trimmean ([1 2 3 4], 10, "flag", 1)
%!error<trimmean: DIM must be a positive integer scalar or vector.> ...
%! trimmean ([1 2 3 4], 10, -1)
%!error<trimmean: DIM must be a positive integer scalar or vector.> ...
%! trimmean ([1 2 3 4], 10, "floor", -1)
%!error<trimmean: DIM must be a positive integer scalar or vector.> ...
%! trimmean (reshape (1:40, [5, 4, 2]), 10, [-1, 2])
%!error<trimmean: VECDIM must contain non-repeating positive integers.> ...
%! trimmean (reshape (1:40, [5, 4, 2]), 10, [1, 2, 2])

