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
## @deftypefn  {statistics} @var{m} = median (@var{x})
## @deftypefnx {statistics} @var{m} = median (@var{x}, "all")
## @deftypefnx {statistics} @var{m} = median (@var{x}, @var{dim})
## @deftypefnx {statistics} @var{m} = median (@var{x}, @var{vecdim})
## @deftypefnx {statistics} @var{m} = median (@dots{}, @var{outtype})
## @deftypefnx {statistics} @var{m} = median (@dots{}, @var{nanflag})
##
## Compute the median value of the elements of @var{x}.
##
## When the elements of @var{x} are sorted, say
## @code{@var{s} = sort (@var{x})}, the median is defined as
## @tex
## $$
## {\rm median} (x) =
##   \cases{s(\lceil N/2\rceil), & $N$ odd;\cr
##           (s(N/2)+s(N/2+1))/2, & $N$ even.}
## $$
## where $N$ is the number of elements of @var{x}.
##
## @end tex
## @ifnottex
##
## @example
## @group
##              |  @var{s}(ceil (N/2))          N odd
## median (@var{x}) = |
##              | (@var{s}(N/2) + @var{s}(N/2+1))/2   N even
## @end group
## @end example
##
## @end ifnottex
## @itemize
## @item
## If @var{x} is a matrix, then @code{median (@var{x})} returns a row vector
## with the mean of each column in @var{x}.
##
## @item
## If @var{x} is a multidimensional array, then @code{median (@var{x})}
## operates along the first non-singleton dimension of @var{x}.
## @end itemize
##
## @code{median (@var{x}, @var{dim})} returns the median along the operating
## dimension @var{dim} of @var{x}.  For @var{dim} greater than
## @code{ndims (@var{x})}, then @var{m} = @var{x}.
##
## @code{median (@var{x}, @var{vecdim})} returns the median over the
## dimensions specified in the vector @var{vecdim}.  For example, if @var{x}
## is a 2-by-3-by-4 array, then @code{median (@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the median of the
## elements on the corresponding page of @var{x}.  If @var{vecdim} indexes all
## dimensions of @var{x}, then it is equivalent to
## @code{median (@var{x}, "all")}.  Any dimension in @var{vecdim} greater than
## @code{ndims (@var{x})} is ignored.
##
## @code{median (@var{x}, "all")} returns the median of all the elements in
## @var{x}.  The optional flag "all" cannot be used together with @var{dim} or
## @var{vecdim} input arguments.
##
## @code{median (@dots{}, @var{outtype})} returns the median with a specified
## data type, using any of the input arguments in the previous syntaxes.
## @var{outtype} can take the following values:
## @table
## @item "default"
## Output is of type double, unless the input is single in which case the output
## is of type single.
##
## @item "double"
## Output is of type double.
##
## @item "native".
## Output is of the same type as the input (@code{class (@var{x})}), unless the
## input is logical in which case the output is of type double.
## @end table
##
## @code{median (@dots{}, @var{nanflag})} specifies whether to exclude NaN
## values from the calculation, using any of the input argument combinations in
## previous syntaxes.  By default, NaN values are included in the calculation
## (@var{nanflag} has the value "includenan").  To exclude NaN values, set the
## value of @var{nanflag} to "omitnan".
##
## @seealso{mean, mode}
## @end deftypefn

function m = median (x, varargin)

  if (nargin < 1 || nargin > 4)
    print_usage ();
  endif

  if (! (isnumeric (x) || islogical (x)))
    error ("median: X must be either numeric or logical");
  endif

  ## Check all char arguments.
  all_flag = false;
  omitnan = false;
  outtype = "default";
  nvarg = numel (varargin);
  varg_chars = cellfun ('ischar', varargin);

  if (nvarg == 2 && ! varg_chars(2))
    print_usage ();
  endif

if (any (varg_chars))
  for idx = varargin(varg_chars)
    switch (tolower (idx{:}))
      case "all"
        all_flag = true;
      case "omitnan"
        omitnan = true;
      case "includenan"
        omitnan = false;
      case {"default", "double", "native"}
        outtype = idx{:};
      otherwise
        print_usage ();
    endswitch
  endfor
  varargin(varg_chars) = [];
  nvarg = numel (varargin);
endif

  if (((nvarg == 1) && ! (isnumeric (varargin{1}))) || (nvarg > 1))
    print_usage ();
  endif

  if (nvarg == 0)

    ## Single numeric input argument, no dimensions given.
    if (all_flag)
      if (omitnan)
        x = x(! isnan (x));
      endif
      n = length (x(:));
      x = sort (x(:), 1);
      k = floor ((n + 1) / 2);
      if (mod (n, 2))
        m = x(k);
      else
        m = x(k) + x(k+1) / 2;
      endif
      ## Inject NaNs where needed, to be consistent with Matlab.
      if (! (omitnan || islogical (x)))
        m(any (isnan (x))) = NaN;
      endif
    else
      sz = size (x);
      (dim = find (sz != 1, 1)) || (dim = 1);
      x = sort (x, dim);
      if (omitnan)
        n = sum (! isnan (x), dim);
      else
        n = sum (ones (sz), dim);
      endif
      k = floor ((n + 1) ./ 2);
      for i = 1:numel (k)
        if (mod (n(i), 2))
          z(i) = {(nth_element (x, k(i), dim))};
        else
          z(i) = {(sum (nth_element (x, k(i):k(i)+1, dim), dim, "native") / 2)};
        endif
      endfor
      ## Collect correct elements
      szargs = cell (1, ndims (x));
      szz = size (z{1});
      for i = 1:numel (k)
        [szargs{:}] = ind2sub (szz, i);
        m(szargs{:}) = z{i}(szargs{:});
      endfor
      ## Inject NaNs where needed, to be consistent with Matlab.
      if (! (omitnan || islogical (x)))
        m(any (isnan (x), dim)) = NaN;
      endif
    endif

  else

    ## Two numeric input arguments, dimensions given.  Note scalar is vector!
    vecdim = varargin{1};
    if (! (isvector (vecdim) && all (vecdim)) || any (rem (vecdim, 1)))
      error ("median: DIM must be a positive integer scalar or vector");
    endif

    if (isscalar (vecdim))
      dim = vecdim;
      x = sort (x, dim);
      if (omitnan)
        n = sum (! isnan (x), dim);
      else
        n = sum (isnan (x) | ! isnan (x), dim);
      endif
      k = floor ((n + 1) ./ 2);
      for i = 1:numel (k)
        if (mod (n(i), 2) == 1)
          z(i) = {(nth_element (x, k(i), dim))};
        else
          z(i) = {(sum (nth_element (x, k(i):k(i)+1, dim), dim, "native") / 2)};
        endif
      endfor
      ## Collect correct elements
      szargs = cell (1, ndims (x));
      szz = size (z{1});
      for i = 1:numel (k)
        [szargs{:}] = ind2sub (szz, i);
        m(szargs{:}) = z{i}(szargs{:});
      endfor
      ## Inject NaNs where needed, to be consistent with Matlab.
      if (! omitnan && ! islogical (x))
        m(any (isnan (x), dim)) = NaN;
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
        if (omitnan)
          x = x(! isnan (x));
        endif
        n = length (x(:));
        x = sort (x(:), 1);
        k = floor ((n + 1) / 2);
        if (mod (n, 2) == 1)
          m = x(k);
        else
          m = x(k) + x(k+1) / 2;
        endif
        ## Inject NaNs where needed, to be consistent with Matlab.
        if (! omitnan && ! islogical (x))
          m(any (isnan (x))) = NaN;
        endif

      else
        ## Permute to bring remaining dims forward
        perm = [remdims, vecdim];
        m = permute (x, perm);

        ## Reshape to put all vecdims in final dimension
        szm = size (m);
        sznew = [szm(1:nremd), prod(szm(nremd+1:end))];
        m = reshape (m, sznew);

        ## Calculate median on single, squashed dimension
        dim = nremd + 1;
        m = sort (m, dim);
        if (omitnan)
          n = sum (! isnan (m), dim);
        else
          n = sum (isnan (m) | ! isnan (m), dim);
        endif
        k = floor ((n + 1) ./ 2);
        for i = 1:numel (k)
          if (mod (n(i), 2) == 1)
            z(i) = {(nth_element (m, k(i), dim))};
          else
            z(i) = {(sum (nth_element (m, k(i):k(i)+1, dim), dim, "native") ...
                     / 2)};
          endif
        endfor
        ## Collect correct elements
        szargs = cell (1, ndims (x));
        szz = size (z{1});
        for i = 1:numel (k)
          [szargs{:}] = ind2sub (szz, i);
          mm(szargs{:}) = z{i}(szargs{:});
        endfor
        ## Inject NaNs where needed, to be consistent with Matlab.
        if (! omitnan && ! islogical (x))
          mm(any (isnan (m), dim)) = NaN;
        endif

        ## Inverse permute back to correct dimensions
        m = ipermute (mm, perm);
      endif
    endif
  endif

  ## Convert output as requested
  switch (outtype)
    case "default"
      ## do nothing, the operators already do the right thing
    case "double"
      m = double (m);
    case "native"
      if (! islogical (x))
        m = cast (m, class (x));
      endif
    otherwise
      error ("mean: OUTTYPE '%s' not recognized", outtype);
  endswitch

endfunction

%!assert (median (1), 1)
%!assert (median ([1,2,3]), 2)
%!assert (median ([3,1,2]), 2)
%!assert (median ([2,4,6,8]), 5)
%!assert (median ([8,2,6,4]), 5)
%!assert (median (single ([1,2,3])), single (2))
%!assert (median ([1,2], 3), [1,2])

%!test
%! x = [1, 2, 3, 4, 5, 6];
%! x2 = x';
%! y = [1, 2, 3, 4, 5, 6, 7];
%! y2 = y';
%!
%! assert (median (x) == median (x2) && median (x) == 3.5);
%! assert (median (y) == median (y2) && median (y) == 4);
%! assert (median ([x2, 2*x2]), [3.5, 7]);
%! assert (median ([y2, 3*y2]), [4, 12]);

## Test outtype option
%!test
%! in = [1 2 3];
%! out = 2;
%! assert (median (in, "default"), median (in));
%! assert (median (in, "default"), out);
%!test
%! in = single ([1 2 3]);
%! out = 2;
%! assert (median (in, "default"), median (in));
%! assert (median (in, "default"), single (out));
%! assert (median (in, "double"), out);
%! assert (median (in, "native"), single (out));
%!test
%! in = uint8 ([1 2 3]);
%! out = 2;
%! assert (median (in, "default"), median (in));
%! assert (median (in, "default"), uint8 (out));
%! assert (median (in, "double"), out);
%! assert (median (in, "native"), uint8 (out));
%!test
%! in = logical ([1 0 1]);
%! out = 1;
%! assert (median (in, "default"), median (in));
%! assert (median (in, "default"), logical (out));
%! assert (median (in, "double"), out);
%! assert (median (in, "native"), logical (out));

## Test single input and optional arguments "all", DIM, "omitnan")
%!test
%! x = repmat ([2 2.1 2.2 2 NaN; 3 1 2 NaN 5; 1 1.1 1.4 5 3], [1, 1, 4]);
%! y = repmat ([2 1.1 2 NaN NaN], [1, 1, 4]);
%! assert (median (x), y);
%! assert (median (x, 1), y);
%! y = repmat ([2 1.1 2 3.5 4], [1, 1, 4]);
%! assert (median (x, "omitnan"), y);
%! assert (median (x, 1, "omitnan"), y);
%! y = repmat ([2.05; 2.5; 1.4], [1, 1, 4]);
%! assert (median (x, 2, "omitnan"), y);
%! y = repmat ([NaN; NaN; 1.4], [1, 1, 4]);
%! assert (median (x, 2), y);
%! assert (median (x, "all"), NaN);
%! assert (median (x, "all", "omitnan"), 3);

# Test boolean input
%!test
%! assert (median (true, "all"), logical (1));
%! assert (median (false), logical (0));
%! assert (median ([true false true]), true);
%! assert (median ([true false true], 2), true);
%! assert (median ([true false true], 1), logical ([1 0 1]));
%! assert (median ([true false NaN], 1), [1 0 NaN]);
%! assert (median ([true false NaN], 2), NaN);
%! assert (median ([true false NaN], 2, "omitnan"), 0.5);
%! assert (median ([true false NaN], 2, "omitnan", "native"), 0.5);

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! assert (size (median (x, [3 2])), [10 1 1 3]);
%! assert (size (median (x, [1 2])), [1 1 6 3]);
%! assert (size (median (x, [1 2 4])), [1 1 6]);
%! assert (size (median (x, [1 4 3])), [1 40]);
%! assert (size (median (x, [1 2 3 4])), [1 1]);

## Test exceeding dimensions
%!assert (median (ones (2,2), 3), ones (2,2));
%!assert (median (ones (2,2,2), 99), ones (2,2,2));
%!assert (median (magic (3), 3), magic (3));
%!assert (median (magic (3), [1 3]), [4, 5, 6]);
%!assert (median (magic (3), [1 99]), [4, 5, 6]);

## Test results with vecdim in n-dimensional arrays and "omitnan"
%!test
%! x = repmat ([2 2.1 2.2 2 NaN; 3 1 2 NaN 5; 1 1.1 1.4 5 3], [1, 1, 4]);
%! assert (median (x, [3 2]), [NaN NaN 1.4]');
%! assert (median (x, [3 2], "omitnan"), [2.05 2.5 1.4]');
%! assert (median (x, [1 3]), [2 1.1 2 NaN NaN]);
%! assert (median (x, [1 3], "omitnan"), [2 1.1 2 3.5 4]);


## Test empty, NaN, Inf inputs
%!assert (median (NaN), NaN)
%!assert <12345> (median (NaN, 'omitnan'), NaN)
%!assert (median (NaN (2)), [NaN NaN])
%!assert <12345> (median (NaN (2), 'omitnan'), [NaN NaN])
%!assert (median ([1 NaN 3]), NaN)
%!assert (median ([1 NaN 3], 1), [1 NaN 3])
%!assert (median ([1 NaN 3], 2), NaN)
%!assert (median ([1 NaN 3]'), NaN)
%!assert (median ([1 NaN 3]', 1), NaN)
%!assert (median ([1 NaN 3]', 2), [1; NaN; 3])
%!assert (median ([1 NaN 3], 'omitnan'), 2)
%!assert (median ([1 NaN 3]', 'omitnan'), 2)
%!assert <12345> (median ([1 NaN 3], 1, 'omitnan'), [1 NaN 3])
%!assert (median ([1 NaN 3], 2, 'omitnan'), 2)
%!assert (median ([1 NaN 3]', 1, 'omitnan'), 2)
%!assert <12345> (median ([1 NaN 3]', 2, 'omitnan'), [1; NaN; 3])
%!assert (median ([1 2 NaN 3]), NaN)
%!assert (median ([1 2 NaN 3], 'omitnan'), 2)
%!assert (median ([1,2,NaN;4,5,6;NaN,8,9]), [NaN, 5, NaN])
%!assert (median ([1 2 ; NaN 4]), [NaN 3])
%!assert (median ([1 2 ; NaN 4], 'omitnan'), [1 3])
%!assert (median ([1 2 ; NaN 4], 1, 'omitnan'), [1 3])
%!assert (median ([1 2 ; NaN 4], 2, 'omitnan'), [1.5; 4], eps)
%!assert <12345> (median ([1 2 ; NaN 4], 3, 'omitnan'), [1 2 ; NaN 4])
%!assert (median ([NaN 2 ; NaN 4]), [NaN 3])
%!assert <12345> (median ([NaN 2 ; NaN 4], 'omitnan'), [NaN 3])
%!assert <12345> (median (ones (1, 0, 3)), NaN (1, 1, 3))

%!assert (median (NaN('single')), NaN('single'));
%!assert <12345> (median (NaN('single'), 'omitnan'), NaN('single'));
%!assert (median (NaN('single'), 'double'), NaN('double'));
%!assert (median (single([1 2 ; NaN 4])), single([NaN 3]));
%!assert <12345> (median (single([1 2 ; NaN 4]), 'double'), double([1 3]));
%!assert (median (single([1 2 ; NaN 4]), 'omitnan'), single([1 3]));
%!assert (median (single([1 2 ; NaN 4]), 'omitnan', 'double'), double([1 3]));
%!assert (median (single([NaN 2 ; NaN 4]), 'double'), double([NaN 3]));
%!assert <12345> (median (single([NaN 2 ; NaN 4]), 'omitnan'), single([NaN 3]));
%!assert <12345> (median (single([NaN 2 ; NaN 4]), 'omitnan', 'double'), double([NaN 3]));

%!assert (median (Inf), Inf);
%!assert (median (-Inf), -Inf);
%!assert (median ([-Inf Inf]), NaN);
%!assert (median ([3 Inf]), Inf);
%!assert (median ([3 4 Inf]), 4);
%!assert (median ([Inf 3 4]), 4);
%!assert (median ([Inf 3 Inf]), Inf);

%!assert <12345> (median ([]), NaN);
%!assert <12345> (median (ones(1,0)), NaN);
%!assert <12345> (median (ones(0,1)), NaN);
%!assert <12345> (median ([], 1), NaN(1,0));
%!assert <12345> (median ([], 2), NaN(0,1));
%!assert <12345> (median ([], 3), NaN(0,1));
%!assert <12345> (median (ones(1,0), 1), NaN(1,0));
%!assert <12345> (median (ones(1,0), 2), NaN(1,1));
%!assert <12345> (median (ones(1,0), 3), NaN(1,0));
%!assert <12345> (median (ones(0,1), 1), NaN(1,1));
%!assert <12345> (median (ones(0,1), 2), NaN(0,1));
%!assert <12345> (median (ones(0,1), 3), NaN(0,1));
%!assert <12345> (median (ones(0,1,0,1), 1), NaN(1,1,0));
%!assert <12345> (median (ones(0,1,0,1), 2), NaN(0,1,0));
%!assert <12345> (median (ones(0,1,0,1), 3), NaN(0,1,1));
%!assert <12345> (median (ones(0,1,0,1), 4), NaN(0,1,0));

## Test complex inputs (should sort by abs(a))
%!assert (median([1 3 3i 2 1i]), 2)
%!assert (median([1 2 4i; 3 2i 4]), [2, 1+1i, 2+2i])

## Test multidimensional arrays
%!shared a, b, x, y
%! old_state = rand ("state");
%! restore_state = onCleanup (@() rand ("state", old_state));
%! rand ("state", 2);
%! a = rand (2,3,4,5);
%! b = rand (3,4,6,5);
%! x = sort (a, 4);
%! y = sort (b, 3);
%!assert <*35679> (median (a, 4), x(:, :, :, 3))
%!assert <*35679> (median (b, 3), (y(:, :, 3, :) + y(:, :, 4, :))/2)
%!shared   ## Clear shared to prevent variable echo for any later test failures

## Test non-floating point types
%!assert <12345> (median ([true, false]), true)
%!assert <12345> (median (logical ([])), logical ([NaN]))
%!assert (median (uint8 ([1, 3])), uint8 (2))
%!assert <12345> (median (uint8 ([])), uint8 ([NaN]))
%!assert (median (uint8 ([NaN 10])), uint8 (5))
%!assert (median (int8 ([1, 3, 4])), int8 (3))
%!assert <12345> (median (int8 ([])), int8 (NaN))
%!assert (median (single ([1, 3, 4])), single (3))
%!assert (median (single ([1, 3, NaN])), single (NaN))

## Test input validation
%!error <Invalid call> median ()
%!error <Invalid call> median (1, 2, 3)
%!error <Invalid call> median (1, 2, 3, 4)
%!error <Invalid call> median (1, "all", 3)
%!error <Invalid call> median (1, "b")
%!error <Invalid call> median (1, 1, "foo")
%!error <X must be either numeric or logical> median ({1:5})
%!error <X must be either numeric or logical> median ("char")
%!error <DIM must be a positive integer> median (1, ones (2,2))
%!error <DIM must be a positive integer> median (1, 1.5)
%!error <DIM must be a positive integer> median (1, 0)

