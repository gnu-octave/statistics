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
## @deftypefn  {statistics} {@var{m} =} mad (@var{x})
## @deftypefnx {statistics} {@var{m} =} mad (@var{x}, @var{flag})
## @deftypefnx {statistics} {@var{m} =} mad (@var{x}, @var{flag}, @qcode{"all"})
## @deftypefnx {statistics} {@var{m} =} mad (@var{x}, @var{flag}, @var{dim})
## @deftypefnx {statistics} {@var{m} =} mad (@var{x}, @var{flag}, @var{vecdim})
##
## Compute the mean or median absolute deviation (MAD).
##
## @code{mad (@var{x})} returns the mean absolute deviation of the values in
## @var{x}.  @code{mad} treats NaNs as missing values and removes them.
##
## @itemize
## @item
## If @var{x} is a vector, then @code{mad} returns the mean or median absolute
## deviation of the values in @var{X}.
## @item
## If @var{x} is a matrix, then @code{mad} returns the mean or median absolute
## deviation of each column of @var{X}.
## @item
## If @var{x} is an multidimensional array, then @code{mad (@var{x})} operates
## along the first non-singleton dimension of @var{x}.
## @end itemize
##
## @code{mad (@var{x}, @var{flag})} specifies whether to compute the mean
## absolute deviation (@qcode{flag = 0}, the default) or the median absolute
## deviation (@qcode{flag = 1}). Passing an empty variable, defaults to 0.
##
## @code{mad (@var{x}, @var{flag}, @qcode{"all"})} returns the MAD of all the
## elements in @var{x}.
##
## The optional variable @var{dim} forces @code{mad} to operate over the
## specified dimension, which must be a positive integer-valued number.
## Specifying any singleton dimension in @var{x}, including any dimension
## exceeding @code{ndims (@var{x})}, will result in a MAD equal to
## @code{zeros (size (@var{x}))}, while non-finite elements are returned as NaNs.
##
## @code{mad (@var{x}, @var{flag}, @var{vecdim})} returns the MAD over the
## dimensions specified in the vector @var{vecdim}.  For example, if @var{x}
## is a 2-by-3-by-4 array, then @code{mad (@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the median of the
## elements on the corresponding page of @var{x}.  If @var{vecdim} indexes all
## dimensions of @var{x}, then it is equivalent to @code{mad (@var{x}, "all")}.
## Any dimension in @var{vecdim} greater than @code{ndims (@var{x})} is ignored.
##
## @seealso{median, mean, mode}
## @end deftypefn

function m = mad (x, flag=0, varargin)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  if (! (isnumeric (x)))
    error ("mad: X must be numeric.");
  endif

  ## Set initial conditions
  all_flag    = false;
  perm_flag   = false;
  vecdim_flag = false;
  dim         = [];

  nvarg = numel (varargin);
  varg_chars = cellfun ("ischar", varargin);
  szx = sz_out = size (x);
  ndx = ndims (x);

  if (isempty (flag))
    flag = 0;
  endif

  if (! (flag == 0 || flag == 1))
    error ("mad: FLAG must be either 0 or 1.");
  endif

  ## Process optional char argument.
  if (any (varg_chars))
    for argin = varargin(varg_chars)
      switch (tolower (argin{:}))
        case "all"
          all_flag = true;

        otherwise
          print_usage ();
      endswitch
    endfor

    varargin(varg_chars) = [];
    nvarg = numel (varargin);
  endif

  ## Process special cases for in/out size
  if (nvarg > 0)
    ## dim or vecdim provided
    dim = varargin{1};
    vecdim_flag = ! isscalar (dim);

    if (! (isvector (dim) && dim > 0) || any (rem (dim, 1)))
      error ("mad: DIM must be a positive integer scalar or vector.");
    endif

    ## Adjust sz_out, account for possible dim > ndx by appending singletons
    sz_out(ndx + 1 : max (dim)) = 1;
    sz_out(dim(dim <= ndx)) = 1;
    szx(ndx + 1 : max (dim)) = 1;

    if (vecdim_flag)
      ## vecdim - try to simplify first
      dim = sort (dim);
      if (! all (diff (dim)))
         error ("mad: VECDIM must contain non-repeating positive integers.");
      endif

      ## dims > ndims(x) and dims only one element long don't affect mad
      sing_dim_x = find (szx != 1);
      dim(dim > ndx | szx(dim) == 1) = [];

      if (isempty (dim))
        ## No dims left to process, return input as output
        m = zeros (sz_out);
        return;
      elseif (numel (dim) == numel (sing_dim_x)
              && unique ([dim, sing_dim_x]) == dim)
        ## If DIMs cover all nonsingleton ndims(x) it's equivalent to "all"
        ##   (check lengths first to reduce unique overhead if not covered)
        all_flag = true;
      endif
    endif

  else
    ## Dim not provided.  Determine scalar dimension.
    if (all_flag)
      ## Special case 'all': Recast input as dim1 vector, process as normal.
      x = x(:);
      szx = [length(x), 1];
      dim = 1;
      sz_out = [1, 1];

    elseif (isrow (x))
      ## Special case row vector: Avoid setting dim to 1.
      dim = 2;
      sz_out = [1, 1];

    elseif (ndx == 2 && szx == [0, 0])
      ## Special case []: Do not apply sz_out(dim)=1 change.
      dim = 1;
      sz_out = [1, 1];

    else
      ## General case: Set dim to first non-singleton, contract sz_out along dim
      (dim = find (szx != 1, 1)) || (dim = 1);
      sz_out(dim) = 1;
    endif
  endif

  if (isempty (x))
    ## Empty input - output NaN
    m = NaN (sz_out);
    return;
  endif

  if (all (isnan (x)) || all (isinf (x)))
    ## NaN or Inf input - output NaN
    m = NaN (sz_out);
    return;
  endif

  if (szx(dim) == 1)
    ## Operation along singleton dimension - nothing to do
    m = zeros (sz_out);
    m(! isfinite (x)) = NaN;
    return;
  endif

  ## Permute dim to simplify all operations along dim1.  At func. end ipermute.
  if (numel (dim) > 1 || (dim != 1 && ! isvector (x)))
    perm = 1 : ndx;

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
  endif

  if (isvector (x))
    if (flag)         # Compute median absolute deviation
      ## Checks above ensure either dim1 or dim2 vector
      x = sort (x, dim);
      x = x(! isnan (x));
      n = length (x);
      k = floor ((n + 1) / 2);
      if (mod (n, 2))
        ## odd
        c = x(k);
        v = sort (abs (x - c));
        m = v(k);
      else
        ## even
        c = (x(k) + x(k + 1)) / 2;
        v = sort (abs (x - c));
        m = (v(k) + v(k + 1)) / 2;
      endif
      m(sum (isinf (x)) > 0) = Inf;

    else              # Compute mean absolute deviation
      x = x(! isnan (x));
      n = length (x);
      m = sum (abs (x - (sum (x) / n))) / n;
      m(sum (isinf (x)) > 0) = Inf;
    endif

  else
    if (flag)         # Compute median absolute deviation
      m = median (abs (x - median(x, dim, "omitnan")), dim, "omitnan");
      m(sum (isinf (x), 2) > 0) = Inf;

    else              # Compute mean absolute deviation
      idx = isnan (x);
      n = sum (! idx, dim);
      x1 = x;
      x1(idx) = 0;
      m1 = sum (x1, dim) ./ n;
      x2 = abs (x - m1);
      x2(idx) = 0;
      m = sum (x2, dim) ./ n;
      m(sum (isinf (x), 2) > 0) = Inf;
    endif
  endif

  if (perm_flag)
    ## Inverse permute back to correct dimensions
    m = ipermute (m, perm);
  endif

endfunction

%!assert (mad (1), 0)
%!assert (mad (1,1), 0)
%!assert (mad (1,0,3), 0)
%!assert (mad (1,1,3), 0)
%!assert (mad (1,[],5), 0)
%!assert (mad ([1,2,3]), 2/3)
%!assert (mad ([1,2,3],[]), 2/3)
%!assert (mad ([1,2,3],0), 2/3)
%!assert (mad ([1,2,3],1), 1)
%!assert (mad ([1,2,3],0,2), 2/3)
%!assert (mad ([1,2,3],[],2), 2/3)
%!assert (mad ([1,2,3],1,2), 1)
%!assert (mad ([1,2,3],0,1), zeros (1,3))
%!assert (mad ([1,2,3],1,1), zeros (1,3))
%!assert (mad ([1,2,3]',0,2), zeros (3,1))
%!assert (mad ([1,2,3]',1,2), zeros (3,1))
%!assert (mad ([1,2,3]',0,1), 2/3)
%!assert (mad ([1,2,3]',1,1), 1)

## Test vector or matrix input with scalar DIM
%!test
%! A = [57, 59, 60, 100, 59, 58, 57, 58, 300, 61, 62, 60, 62, 58, 57];
%! AA = [A;2*A;3*A];
%! m0 = [38.000, 39.333, 40.000, 66.667, 39.333, 38.667, 38.000, 38.667, ...
%!       200.000, 40.667, 41.333, 40.000, 41.333, 38.667, 38.000];
%! m1 = [32.569;65.138; 97.707];
%!
%! assert (mad (AA), m0, 1e-3);
%! assert (mad (AA,1), A);
%! assert (mad (AA,1,1), A);
%! assert (mad (AA,0,2), m1, 1e-3);
%! assert (mad (AA,1,2), [2;4;6]);
%! assert (mad (A,0,1), zeros (size (A)));
%! assert (mad (A,1,1), zeros (size (A)));

## Test n-dimensional input and optional arguments "all", VECDIM
%!test
%! x = repmat ([2 2.1 2.2 2 NaN; 3 1 2 NaN 5; 1 1.1 1.4 5 3], [1, 1, 4]);
%! m0 = repmat ([0.6667, 0.4667, 0.3111, 1.5, 1], [1, 1, 4]);
%! m1 = repmat ([1, 0.1, 0.2, 1.5, 1], [1, 1, 4]);
%! assert (mad (x), m0, 1e-4);
%! assert (mad (x, 1), m1, 1e-14);
%! assert (mad (x, [], [1, 2]), 1.0036 * ones(1,1,4), 1e-4)
%! assert (mad (x, 1, [1, 2]), 0.9 * ones(1,1,4), 1e-14)
%! assert (mad (x, 0, [1, 3]), m0(1,:,1), 1e-4)
%! assert (mad (x, 1, [1, 3]), m1(1,:,1), 1e-14)
%! assert (mad (x, 0, [2, 3]), [0.075; 1.25; 1.36], 1e-14)
%! assert (mad (x, 1, [2, 3]), [0.05; 1; 0.4], 1e-14)
%! assert (mad (x, 0, [1, 2, 3]) == mad (x, 0, "All"))
%! assert (mad (x, 1, [1, 2, 3]) == mad (x, 1, "All"))

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! assert (size (mad (x, [], [3 2])), [10 1 1 3]);
%! assert (size (mad (x, 0, [1 2])), [1 1 6 3]);
%! assert (size (mad (x, 1, [1 2 4])), [1 1 6]);
%! assert (size (mad (x, [], [1 4 3])), [1 40]);
%! assert (size (mad (x, 1, [1 2 3 4])), [1 1]);

## Test exceeding dimensions
%!assert (mad (ones (2,2), 0, 3), zeros (2,2))
%!assert (mad (ones (2,2,2), 1, 99), zeros (2,2,2))
%!assert (mad (magic (3), 1, 3), zeros (3))
%!assert (mad (magic (3), 1, [1 3]), [1, 4, 1])
%!assert (mad (magic (3), 1, [1 99]), [1, 4, 1])

## Test empty, NaN, Inf inputs
%!assert (mad ([]), NaN)
%!assert (mad ([], 1), NaN)
%!assert (mad (NaN), NaN)
%!assert (mad (NaN, 1), NaN)
%!assert (mad (Inf), NaN)
%!assert (mad (Inf, 1), NaN)
%!assert (mad (-Inf), NaN)
%!assert (mad (-Inf, 1), NaN)
%!assert (mad ([-Inf Inf]), NaN)
%!assert (mad ([-Inf Inf], 1), NaN)
%!assert (mad ([3 Inf]), Inf)
%!assert (mad ([3 4 Inf]), Inf)
%!assert (mad ([3 4 Inf], 1), Inf)
%!assert (mad ([Inf 3 4]), Inf)
%!assert (mad ([Inf 3 4], 1), Inf)
%!assert (mad ([Inf 3 Inf]), Inf)
%!assert (mad ([Inf 3 Inf], 1), Inf)
%!assert (mad ([1 2; 3 Inf]), [1 Inf])
%!assert (mad ([1 2; 3 Inf], 1), [1 Inf])

%!assert (mad ([]), NaN)
%!assert (mad (ones(1,0)), NaN)
%!assert (mad (ones(0,1)), NaN)
%!assert (mad ([], 0, 1), NaN(1,0))
%!assert (mad ([], 0, 2), NaN(0,1))
%!assert (mad ([], 0, 3), NaN(0,0))
%!assert (mad (ones(1,0), 0, 1), NaN(1,0))
%!assert (mad (ones(1,0), 0, 2), NaN(1,1))
%!assert (mad (ones(1,0), 0, 3), NaN(1,0))
%!assert (mad (ones(0,1), 0, 1), NaN(1,1))
%!assert (mad (ones(0,1), 0, 2), NaN(0,1))
%!assert (mad (ones(0,1), 0, 3), NaN(0,1))
%!assert (mad (ones(0,1,0,1), 0, 1), NaN(1,1,0))
%!assert (mad (ones(0,1,0,1), 0, 2), NaN(0,1,0))
%!assert (mad (ones(0,1,0,1), 0, 3), NaN(0,1,1))
%!assert (mad (ones(0,1,0,1), 0, 4), NaN(0,1,0))

## Test complex inputs (should sort by abs(a))
%!assert (mad([1 3 3i 2 1i]), 1.5297, 1e-4)
%!assert (mad([1 3 3i 2 1i], 1), 1)
%!assert (mad([1 2 4i; 3 2i 4]), [1, 1.4142, 2.8284], 1e-4)
%!assert (mad([1 2 4i; 3 2i 4], 1), [1, 1.4142, 2.8284], 1e-4)
%!assert (mad([1 2 4i; 3 2i 4], 1, 2), [1; 1])
%!assert (mad([1 2 4i; 3 2i 4], 0, 2), [1.9493; 1.8084], 1e-4)

## Test all-inf handling
%!assert <*65405> (mad ([-Inf Inf]), NaN)
%!assert <*65405> (mad ([-Inf Inf], 0), NaN)
%!assert <*65405> (mad ([-Inf Inf], 1), NaN)
%!assert <*65405> (mad ([-Inf Inf]', 0), NaN)
%!assert <*65405> (mad ([-Inf Inf]', 1), NaN)
%!assert <*65405> (mad ([-Inf Inf]', 0, 1), NaN)
%!assert <*65405> (mad ([-Inf Inf]', 0, 2), [NaN; NaN])
%!assert <*65405> (mad ([-Inf Inf]', 0, 3), [NaN; NaN])
%!assert <*65405> (mad ([-Inf Inf]', 1, 1), NaN)
%!assert <*65405> (mad ([-Inf Inf]', 1, 2), [NaN; NaN])
%!assert <*65405> (mad ([-Inf Inf]', 1, 3), [NaN; NaN])
%!assert <*65405> (mad (Inf(2), 0), [NaN, NaN])
%!assert <*65405> (mad (Inf(2), 1), [NaN, NaN])
%!assert <*65405> (mad (Inf(2), 0, 1), [NaN, NaN])
%!assert <*65405> (mad (Inf(2), 0, 2), [NaN; NaN])
%!assert <*65405> (mad (Inf(2), 0, 3), NaN(2))
%!assert <*65405> (mad (Inf(2), 1, 1), [NaN, NaN])
%!assert <*65405> (mad (Inf(2), 1, 2), [NaN; NaN])
%!assert <*65405> (mad (Inf(2), 1, 3), NaN(2))

## Test input case insensitivity
%!assert (mad ([1 2 3], 0, "aLL"), 2/3)
%!assert (mad ([1 2 3], 1, "aLL"), 1)

## Test input validation
%!error <Invalid call> mad ()
%!error <Invalid call> mad (1, 2, 3, 4)
%!error <mad: X must be numeric.> mad ("text")
%!error <mad: X must be numeric.> mad ({2 3 4})
%!error <mad: FLAG must be either 0 or 1.> mad (1, "all", 3)
%!error <mad: FLAG must be either 0 or 1.> mad (1, "b")
%!error <Invalid call> mad (1, 1, "foo")
%!error <DIM must be a positive integer> mad (1, [] ,ones (2,2))
%!error <DIM must be a positive integer> mad (1, [], 1.5)
%!error <DIM must be a positive integer> mad (1, [], 0)
%!error <DIM must be a positive integer> mad ([1 2 3], [], [-1 1])
%!error <VECDIM must contain non-repeating> mad(1, [], [1 2 2])
