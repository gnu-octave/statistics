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
##
## @end tex
## @ifnottex
##
## @example
## harmmean (@var{x}) = N / SUM_i @var{x}(i)^-1
## @end example
##
## @end ifnottex
## @noindent
## where @math{N} is the length of the @var{x} vector.
##
## @item If @var{x} is a matrix, then @code{harmmean(@var{x})} returns a row
## vector with the harmonic mean of each columns in @var{x}.
##
## @item If @var{x} is a multidimensional array, then @code{harmmean(@var{x})}
## operates along the first nonsingleton dimension of @var{x}.
##
## @item @var{x} must not contain any negative or complex values.
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
## a 2-by-3-by-4 array, then @code{harmmean(@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the harmonic mean of
## the elements on the corresponding page of @var{x}.  If @var{vecdim} indexes
## all dimensions of @var{x}, then it is equivalent to @code{harmmean (@var{x},
## "all")}.  Any dimension in @var{vecdim} greater than @code{ndims (@var{x})}
## is ignored.
##
## @code{harmmean(@dots{}, @var{nanflag})} specifies whether to exclude NaN
## values from the calculation, using any of the input argument combinations in
## previous syntaxes. By default, harmmean includes NaN values in the
## calculation (@var{nanflag} has the value "includenan").  To exclude NaN
## values, set the value of @var{nanflag} to "omitnan".
##
## @seealso{geomean, mean}
## @end deftypefn

function m = harmmean (x, varargin)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  if (! isnumeric (x) || ! isreal (x) || ! all (x(! isnan (x))(:) >= 0))
    error ("harmmean: X must contain real nonnegative values.");
  endif

  ## Set initial conditions
  all_flag = false;
  omitnan = false;

  nvarg = numel (varargin);
  varg_chars = cellfun ("ischar", varargin);
  szx = size (x);
  ndx = ndims (x);

  if (nvarg > 1 && ! varg_chars(2:end))
    ## Only first varargin can be numeric
    print_usage ();
  endif

  ## Process any other char arguments.
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

  ## Single numeric input argument, no dimensions given.
  if (nvarg == 0)

    if (all_flag)
      x = x(:);

      if (omitnan)
        x = x(! isnan (x));
      endif

      if (any (x == 0))
        m = 0;
        return;
      endif

      m = length (x) ./ sum (1 ./ x);
      m(m == Inf) = 0;  # handle zeros in X

    else
      ## Find the first non-singleton dimension.
      (dim = find (szx != 1, 1)) || (dim = 1);
      n = szx(dim);
      is_nan = 0;
      if (omitnan)
        idx = isnan (x);
        n = sum (! idx, dim);
        is_nan = sum (idx, dim);
        x(idx) = 1;     # remove NaNs by subtracting is_nan below
      endif

      m = n ./ (sum (1 ./ x, dim) - is_nan);
      m(m == Inf) = 0;  # handle zeros in X

    endif

  else

    ## Two numeric input arguments, dimensions given.  Note scalar is vector!
    vecdim = varargin{1};
    if (isempty (vecdim) || ! (isvector (vecdim) && all (vecdim > 0)) ...
          || any (rem (vecdim, 1)))
      error ("harmmean: DIM must be a positive integer scalar or vector.");
    endif

    if (ndx == 2 && isempty (x) && szx == [0,0])
      ## FIXME: this special case handling could be removed once sum
      ##        compatibly handles all sizes of empty inputs
      sz_out = szx;
      sz_out (vecdim(vecdim <= ndx)) = 1;
      m = NaN (sz_out);
    else

      if (isscalar (vecdim))
        if (vecdim > ndx)
          m = x;
        else
          n = szx(vecdim);
          is_nan = 0;
          if (omitnan)
            nanx = isnan (x);
            n = sum (! nanx, vecdim);
            is_nan = sum (nanx, vecdim);
            x(nanx) = 1;    # remove NaNs by subtracting is_nan below
          endif

          m = n ./ (sum (1 ./ x, vecdim) - is_nan);
          m(m == Inf) = 0;  # handle zeros in X

        endif

      else
        vecdim = sort (vecdim);
        if (! all (diff (vecdim)))
           error (strcat ("harmmean: VECDIM must contain non-repeating", ...
                          " positive integers."));
        endif
        ## Ignore exceeding dimensions in VECDIM
        vecdim(find (vecdim > ndims (x))) = [];

        if (isempty (vecdim))
          m = x;
        else
          ## Move vecdims to dim 1.

          ## Calculate permutation vector
          remdims = 1 : ndx;        # All dimensions
          remdims(vecdim) = [];     # Delete dimensions specified by vecdim
          nremd = numel (remdims);

          ## If all dimensions are given, it is similar to all flag
          if (nremd == 0)
            x = x(:);

            if (omitnan)
              x = x(! isnan (x));
            endif

            if (any (x == 0))
              m = 0;
              return;
            endif

            m = length (x) ./ sum (1 ./ x);
            m(m == Inf) = 0;  # handle zeros in X

          else
            ## Permute to bring vecdims to front
            perm = [vecdim, remdims];
            x = permute (x, perm);

            ## Reshape to squash all vecdims in dim1
            num_dim = prod (szx(vecdim));
            szx(vecdim) = [];
            szx = [ones(1, length(vecdim)), szx];
            szx(1) = num_dim;
            x = reshape (x, szx);

            ## Calculate mean on dim1
            if (omitnan)
              nanx = isnan (x);
              n = sum (! nanx, 1);
              is_nan = sum (nanx, 1);
              x(nanx) = 1;    # remove NaNs by subtracting is_nan below
            else
              n = szx(1);
              is_nan = 0;
            endif

            m = n ./ (sum (1 ./ x, 1) - is_nan);
            m(m == Inf) = 0;  # handle zeros in X

            ## Inverse permute back to correct dimensions
            m = ipermute (m, perm);
          endif
        endif
      endif
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
%! assert (harmmean (y, 2), [0 NaN m(3)]', 4e-14);
%! assert (harmmean (y', "omitnan"), m, 4e-14);
%! z = y + 20;
%! assert (harmmean (z, "all"), NaN);
%! assert (harmmean (z, "all", "includenan"), NaN);
%! assert (harmmean (z, "all", "omitnan"), 29.1108719858295, 4e-14);
%! m = [24.59488458841874 NaN 34.71244385944397];
%! assert (harmmean (z'), m, 4e-14);
%! assert (harmmean (z', "includenan"), m, 4e-14);
%! m(2) = 29.84104075528277;
%! assert (harmmean (z', "omitnan"), m, 4e-14);
%! assert (harmmean (z, 2, "omitnan"), m', 4e-14);

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! assert (size (harmmean (x, [3 2])), [10 1 1 3]);
%! assert (size (harmmean (x, [1 2])), [1 1 6 3]);
%! assert (size (harmmean (x, [1 2 4])), [1 1 6]);
%! assert (size (harmmean (x, [1 4 3])), [1 40]);
%! assert (size (harmmean (x, [1 2 3 4])), [1 1]);

## Test results with vecdim in n-dimensional arrays and "omitnan"
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! m = repmat ([5.559045930488016;13.04950789021461], [5 1 1 3]);
%! assert (harmmean (x, [3 2]), m, 4e-14);
%! x(2,5,6,3) = NaN;
%! m(2,3) = NaN;
%! assert (harmmean (x, [3 2]), m, 4e-14);
%! m(2,3) = 13.06617961315406;
%! assert (harmmean (x, [3 2], "omitnan"), m, 4e-14);

## Test errors
%!error <harmmean: X must contain real nonnegative values.> harmmean ("char")
%!error <harmmean: X must contain real nonnegative values.> harmmean ([1 -1 3])
%!error <harmmean: DIM must be a positive integer scalar or vector.> ...
%! harmmean (repmat ([1:20;6:25], [5 2 6 3 5]), -1)
%!error <harmmean: DIM must be a positive integer scalar or vector.> ...
%! harmmean (repmat ([1:20;6:25], [5 2 6 3 5]), 0)
%!error <harmmean: VECDIM must contain non-repeating positive integers.> ...
%! harmmean (repmat ([1:20;6:25], [5 2 6 3 5]), [1 1])
