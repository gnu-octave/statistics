## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2023 Nick Jankowski
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
## @itemize
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
## @end itemize
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

### TODO: check relative speed of using nth_element

  if (nargin < 1 || nargin > 4)
    print_usage ();
  endif

  if (! (isnumeric (x) || islogical (x)))
    error ("median: X must be either numeric or logical.");
  endif

  ## Set initial conditions
  all_flag = 0;
  omitnan = 0;
  perm_flag = 0;
  out_flag = 0;
  vecdim_flag = 0;
  dim = [];

  nvarg = numel (varargin);
  varg_chars = cellfun ("ischar", varargin);
  szx = sz_out = size (x);
  ndx = ndims (x);
  outtype = class (x);

  if (nvarg > 1 && ! varg_chars(2:end))
    ## Only first varargin can be numeric
    print_usage ();
  endif

  ## Process any other char arguments.
  if (any (varg_chars))
    for idx = varargin(varg_chars)
      switch (tolower (idx{:}))
        case "all"
          all_flag = 1;

        case "omitnan"
          omitnan = 1;

        case "includenan"
          omitnan = 1;

        case "native"
          if (out_flag)
            error ("median: only one OUTTYPE can be specified.")
          endif
          if (strcmp (outtype, "logical"))
            outtype = "double";
          endif
          out_flag = 1;

        case "default"
          if (out_flag)
            error ("median: only one OUTTYPE can be specified.")
          endif
          if (! strcmp (outtype, "single"))
            outtype = "double";
          endif
          out_flag = 1;

        case "double"
          if (out_flag)
            error ("median: only one OUTTYPE can be specified.")
          endif
          outtype = "double";
          out_flag = 1;

        otherwise
          print_usage ();
      endswitch
    endfor

    varargin(varg_chars) = [];
    nvarg = numel (varargin);
  endif

  if (((nvarg == 1) && ! (isnumeric (varargin{1}))) || (nvarg > 1))
    ## After trimming char inputs can only be one varargin left, must be numeric
    print_usage ();
  endif

  ## Process special cases for in/out size
  if (nvarg > 0)
    ## dim or vecdim provided
    if (all_flag)
      error ("median: 'all' cannot be used with DIM or VECDIM options.");
    endif

    dim = varargin{1};
    vecdim_flag = ! isscalar (dim);

    if (! (isvector (dim) && (dim > 0)) || any (rem (dim, 1)))
      error ("median: DIM must be a positive integer scalar or vector.");
    endif

    ## Adjust sz_out, account for possible dim > ndx by appending singletons
    sz_out(ndx + 1 : max (dim)) = 1;
    sz_out(dim(dim <= ndx)) = 1;
    szx(ndx + 1 : max (dim)) = 1;

    if (vecdim_flag)
      ## vecdim - try to simplify first
      dim = sort (dim);
      if (! all (diff (dim)))
         error ("median: VECDIM must contain non-repeating positive integers.");
      endif

      ## dims > ndims(x) and dims only one element long don't affect median
      sing_dim_x = find (szx != 1);
      dim(dim > ndx | szx(dim) == 1) = [];

      if (isempty (dim))
        ## No dims left to process, return input as output
        if (! strcmp (class (x), outtype))
          m = outtype_convert (x, outtype);
        else
          m = x;
        endif
        return;
      elseif ((length(dim) == length(sing_dim_x)) ...
                 && unique ([dim, sing_dim_x]) == dim)
        ## If DIMs cover all nonsingleton ndims(x) it's equivalent to "all"
        ##   (check lengths first to reduce unique overhead if not covered)
        all_flag = true;
      endif
    endif

  else
    ## Dim not provided.  Determine scalar dimension
    if (all_flag)
      ## Special case 'all': Recast input as dim1 vector, process as normal
      x = x(:);
      szx = [length(x), 1];
      dim = 1;
      sz_out = [1 1];

    elseif (isrow (x))
      ## Special case row vector: Avoid setting dim to 1.
      dim = 2;
      sz_out = [1, 1];

    elseif (ndx == 2 && szx == [0, 0])
      ## Special case []: Do not apply sz_out(dim)=1 change
      dim = 1;
      sz_out = [1, 1];

    else
      ## General case: Set dim to first non-singleton, contract sz_out along dim
      (dim = find (szx != 1, 1)) || (dim = 1);
      sz_out(dim) = 1;
    endif
  endif

  if (isempty (x))
    ## Empty input - output NaN or class equivalent in pre-determined size
    switch (outtype)
      case {"double", "single"}
        m = NaN (sz_out, outtype);
      case ("logical")
        m = false (sz_out);
      otherwise
        m = cast (NaN (sz_out), outtype);
    endswitch
    return;
  endif

  if (szx(dim) == 1)
    ## Operation along singleton dimension - nothing to do
    if (! strcmp (class (x), outtype))
      m = outtype_convert (x, outtype);
    else
      m = x;
    endif
    return;
  endif

  ## Permute dim to simplify all operations along dim1.  At func. end ipermute.
  if ((length (dim) > 1) || (dim != 1 && ! isvector (x)))
    perm_vect = 1 : ndx;

    if (! vecdim_flag)
      ## Move dim to dim 1
      perm_vect([1, dim]) = [dim, 1];
      x = permute (x, perm_vect);
      szx([1, dim]) = szx([dim, 1]);
      dim = 1;

    else
      ## Move vecdims to front
      perm_vect(dim) = [];
      perm_vect = [dim, perm_vect];
      x = permute (x, perm_vect);

      ## Reshape all vecdims into dim1
      num_dim = prod (szx(dim));
      szx(dim) = [];
      szx = [ones(1, length(dim)), szx];
      szx(1) = num_dim;
      x = reshape (x, szx);
      dim = 1;
    endif

    perm_flag = true;
  endif

  ## Find column locations of NaNs
  hasnan = any (isnan (x), dim);
  if (! hasnan(:) && omitnan)
    ## Don't use omitnan path if no NaNs are present
    omitnan = 0;
  endif

  x = sort (x, dim); # Note- pushes any NaN's to end

  if (omitnan)
    ## Ignore any NaN's in data. Each operating vector might have a
    ## different number of non-NaN data points.

    if (isvector (x))
      ## Checks above ensure either dim1 or dim2 vector
      x = x(! isnan (x));
      n = length (x);
      k = floor ((n + 1) / 2);
      if (mod (n, 2))
        ## odd
        m = x(k);
      else
        ## even
        m = (x(k) + x(k + 1)) / 2;
      endif

    else
      n = sum (! isnan (x), 1);
      k = floor ((n + 1) ./ 2);
      m_idx_odd = mod (n, 2) & n;
      m_idx_even = ! m_idx_odd & n;

      m = NaN ([1, szx(2 : end)]);

      if (ndims (x) > 2)
        szx = [szx(1), prod(szx(2 : end))];
      endif

      ## Grab kth value, k possibly different for each column
      x_idx_odd = sub2ind (szx, (k(m_idx_odd))(:)', ...
                             (1 : szx(2))(m_idx_odd)(:)');
      x_idx_even = sub2ind (szx, ...
                             [(k(m_idx_even))(:)'; (k(m_idx_even) + 1)(:)'], ...
                               (1 : szx(2))(m_idx_even)([1 1], :));

      m(m_idx_odd) = x(x_idx_odd);
      m(m_idx_even) = sum (x(x_idx_even), 1) / 2;
    endif

  else
    ## No "omitnan". All 'vectors' uniform length.
    if (! all (hasnan))

      if (isvector (x))
        n = length (x);
        k = floor ((n + 1) / 2);
        if (mod (n, 2))
          ## Odd
          m = x(k);
        else
          ## Even
          m = (x(k) + x(k + 1)) / 2;
        endif

      else
        ## Nonvector, all operations permuted to be along dim 1
        n = szx(1);
        k = floor ((n + 1) / 2);
        m = NaN ([1, szx(2 : end)]);
        if (mod (n, 2))
          ## Odd
          m(1, :) = x(k, :);
        else
          ## Even
          m(1, :) = (x(k, :) + x(k + 1, :)) / 2;
        endif
      endif
      if (any (hasnan(:)))
        m(hasnan) = NaN;
      endif
    else
      m = NaN (sz_out);
    endif
  endif

  if (perm_flag)
    ## Inverse permute back to correct dimensions
    m = ipermute (m, perm_vect);
  endif

  ## Convert output type as requested
  if (! strcmp (class (m), outtype))
    m = outtype_convert (m, outtype);
  endif

endfunction

function m = outtype_convert (m, outtype)
  switch (outtype)
    case "single"
      m = single (m);
    case "double"
      m = double (m);
    otherwise
      m = cast (m, outtype);
  endswitch
endfunction

%!assert (median (1), 1)
%!assert (median ([1,2,3]), 2)
%!assert (median ([1,2,3]'), 2)
%!assert (median (cat(3,3,1,2)), 2)
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
%! assert (median ([x2, 2 * x2]), [3.5, 7]);
%! assert (median ([y2, 3 * y2]), [4, 12]);

## Test outtype option
%!test
%! in = [1 2 3];
%! out = 2;
%! assert (median (in, "default"), median (in));
%! assert (median (in, "default"), out);
%!test
%! in = single ([1 2 3]);
%! out = 2;
%! assert (median (in, "default"), single (median (in)));
%! assert (median (in, "default"), single (out));
%! assert (median (in, "double"), double (out));
%! assert (median (in, "native"), single (out));
%!test
%! in = uint8 ([1 2 3]);
%! out = 2;
%! assert (median (in, "default"), double (median (in)));
%! assert (median (in, "default"), double (out));
%! assert (median (in, "double"), out);
%! assert (median (in, "native"), uint8 (out));
%!test
%! in = logical ([1 0 1]);
%! out = 1;
%! assert (median (in, "default"), double (median (in)));
%! assert (median (in, "default"), double (out));
%! assert (median (in, "double"), double (out));
%! assert (median (in, "native"), double (out));

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
%! assert (median (x, "all", "omitnan"), 2);
%!assert (median (cat (3, 3, 1, NaN, 2), "omitnan"), 2)
%!assert (median (cat (3, 3, 1, NaN, 2), 3, "omitnan"), 2)

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
%! assert (median ([true false NaN], 2, "omitnan", "native"), double(0.5));

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
%!assert (median (NaN, "omitnan"), NaN)
%!assert (median (NaN (2)), [NaN NaN])
%!assert (median (NaN (2), "omitnan"), [NaN NaN])
%!assert (median ([1 NaN 3]), NaN)
%!assert (median ([1 NaN 3], 1), [1 NaN 3])
%!assert (median ([1 NaN 3], 2), NaN)
%!assert (median ([1 NaN 3]'), NaN)
%!assert (median ([1 NaN 3]', 1), NaN)
%!assert (median ([1 NaN 3]', 2), [1; NaN; 3])
%!assert (median ([1 NaN 3], "omitnan"), 2)
%!assert (median ([1 NaN 3]', "omitnan"), 2)
%!assert (median ([1 NaN 3], 1, "omitnan"), [1 NaN 3])
%!assert (median ([1 NaN 3], 2, "omitnan"), 2)
%!assert (median ([1 NaN 3]', 1, "omitnan"), 2)
%!assert (median ([1 NaN 3]', 2, "omitnan"), [1; NaN; 3])
%!assert (median ([1 2 NaN 3]), NaN)
%!assert (median ([1 2 NaN 3], "omitnan"), 2)
%!assert (median ([1,2,NaN;4,5,6;NaN,8,9]), [NaN, 5, NaN])
%!assert (median ([1 2 ; NaN 4]), [NaN 3])
%!assert (median ([1 2 ; NaN 4], "omitnan"), [1 3])
%!assert (median ([1 2 ; NaN 4], 1, "omitnan"), [1 3])
%!assert (median ([1 2 ; NaN 4], 2, "omitnan"), [1.5; 4], eps)
%!assert (median ([1 2 ; NaN 4], 3, "omitnan"), [1 2 ; NaN 4])
%!assert (median ([NaN 2 ; NaN 4]), [NaN 3])
%!assert (median ([NaN 2 ; NaN 4], "omitnan"), [NaN 3])
%!assert (median (ones (1, 0, 3)), NaN (1, 1, 3))

%!assert (median (NaN("single")), NaN("single"));
%!assert (median (NaN("single"), "omitnan"), NaN("single"));
%!assert (median (NaN("single"), "double"), NaN("double"));
%!assert (median (single([1 2 ; NaN 4])), single([NaN 3]));
%!assert (median (single([1 2 ; NaN 4]), "double"), double([NaN 3]));
%!assert (median (single([1 2 ; NaN 4]), "omitnan"), single([1 3]));
%!assert (median (single([1 2 ; NaN 4]), "omitnan", "double"), double([1 3]));
%!assert (median (single([NaN 2 ; NaN 4]), "double"), double([NaN 3]));
%!assert (median (single([NaN 2 ; NaN 4]), "omitnan"), single([NaN 3]));
%!assert (median (single([NaN 2 ; NaN 4]), "omitnan", "double"), double([NaN 3]));

%!assert (median (Inf), Inf);
%!assert (median (-Inf), -Inf);
%!assert (median ([-Inf Inf]), NaN);
%!assert (median ([3 Inf]), Inf);
%!assert (median ([3 4 Inf]), 4);
%!assert (median ([Inf 3 4]), 4);
%!assert (median ([Inf 3 Inf]), Inf);

%!assert (median ([]), NaN);
%!assert (median (ones(1,0)), NaN);
%!assert (median (ones(0,1)), NaN);
%!assert (median ([], 1), NaN(1,0));
%!assert (median ([], 2), NaN(0,1));
%!assert (median ([], 3), NaN(0,0));
%!assert (median (ones(1,0), 1), NaN(1,0));
%!assert (median (ones(1,0), 2), NaN(1,1));
%!assert (median (ones(1,0), 3), NaN(1,0));
%!assert (median (ones(0,1), 1), NaN(1,1));
%!assert (median (ones(0,1), 2), NaN(0,1));
%!assert (median (ones(0,1), 3), NaN(0,1));
%!assert (median (ones(0,1,0,1), 1), NaN(1,1,0));
%!assert (median (ones(0,1,0,1), 2), NaN(0,1,0));
%!assert (median (ones(0,1,0,1), 3), NaN(0,1,1));
%!assert (median (ones(0,1,0,1), 4), NaN(0,1,0));

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
%!assert (median ([true, false]), true)
%!assert (median (logical ([])), false)
%!assert (median (uint8 ([1, 3])), uint8 (2))
%!assert (median (uint8 ([])), uint8 (NaN))
%!assert (median (uint8 ([NaN 10])), uint8 (5))
%!assert <54567> (median (uint8 ([253, 255])), uint8 (254))
%!assert <54567> (median (uint8 ([253, 254])), uint8 (254))
%!assert (median (int8 ([1, 3, 4])), int8 (3))
%!assert (median (int8 ([])), int8 (NaN))
%!assert <54567> (median (int8 ([127, 126, 125, 124; 1 3 5 9])), int8 ([64 65 65 67]))
%!assert <54567> (median (int8 ([127, 126, 125, 124; 1 3 5 9]), 2), int8 ([126; 4]))
%!assert (median (single ([1, 3, 4])), single (3))
%!assert (median (single ([1, 3, NaN])), single (NaN))

## Test input case insensitivity
%!assert (median ([1 2 3], "aLL"), 2);
%!assert (median ([1 2 3], "OmitNan"), 2);
%!assert (median ([1 2 3], "DOUBle"), 2);

## Test input validation
%!error <Invalid call> median ()
%!error <Invalid call> median (1, 2, 3)
%!error <Invalid call> median (1, 2, 3, 4)
%!error <Invalid call> median (1, "all", 3)
%!error <Invalid call> median (1, "b")
%!error <Invalid call> median (1, 1, "foo")
%!error <'all' cannot be used with> median (1, 3, "all")
%!error <'all' cannot be used with> median (1, [2 3], "all")
%!error <X must be either numeric or logical> median ({1:5})
%!error <X must be either numeric or logical> median ("char")
%!error <only one OUTTYPE can be specified> median(1, "double", "native")
%!error <DIM must be a positive integer> median (1, ones (2,2))
%!error <DIM must be a positive integer> median (1, 1.5)
%!error <DIM must be a positive integer> median (1, 0)
%!error <DIM must be a positive integer> median ([1 2 3], [-1 1])
%!error <VECDIM must contain non-repeating> median(1, [1 2 2])
