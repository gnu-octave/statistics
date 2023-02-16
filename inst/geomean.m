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
## @deftypefn  {statistics} {@var{m} =} geomean (@var{x})
## @deftypefnx {statistics} {@var{m} =} geomean (@var{x}, "all")
## @deftypefnx {statistics} {@var{m} =} geomean (@var{x}, @var{dim})
## @deftypefnx {statistics} {@var{m} =} geomean (@var{x}, @var{vecdim})
## @deftypefnx {statistics} {@var{m} =} geomean (@dots{}, @var{nanflag})
##
## Compute the geometric mean of @var{x}.
##
## @itemize
## @item If @var{x} is a vector, then @code{geomean(@var{x})} returns the
## geometric mean of the elements in @var{x} defined as
## @tex
## $$ {\rm geomean}(x) = \left( \prod_{i=1}^N x_i \right)^\frac{1}{N}
## = exp \left({1\over N} \sum_{i=1}^N log x_i \right) $$
##
## @end tex
## @ifnottex
##
## @example
## geomean (@var{x}) = PROD_i @var{x}(i) ^ (1/N)
## @end example
##
## @end ifnottex
## @noindent
## where @math{N} is the length of the @var{x} vector.
##
## @item If @var{x} is a matrix, then @code{geomean(@var{x})} returns a row
## vector with the geometric mean of each columns in @var{x}.
##
## @item If @var{x} is a multidimensional array, then @code{geomean(@var{x})}
## operates along the first nonsingleton dimension of @var{x}.
##
## @item @var{x} must not contain any negative or complex values.
## @end itemize
##
## @code{geomean(@var{x}, "all")} returns the geometric mean of all the elements
## in @var{x}.  If @var{x} contains any 0, then the returned value is 0.
##
## @code{geomean(@var{x}, @var{dim})} returns the geometric mean along the
## operating dimension @var{dim} of @var{x}.  Calculating the harmonic mean of
## any subarray containing any 0 will return 0.
##
## @code{geomean(@var{x}, @var{vecdim})} returns the geometric mean over the
## dimensions specified in the vector @var{vecdim}.  For example, if @var{x} is
## a 2-by-3-by-4 array, then @code{geomean(@var{x}, [1 2])} returns a
## 1-by-1-by-4 array.  Each element of the output array is the geometric mean of
## the elements on the corresponding page of @var{x}.  If @var{vecdim} indexes
## all dimensions of @var{x}, then it is equivalent to @code{geomean (@var{x},
## "all")}.  Any dimension in @var{vecdim} greater than @code{ndims (@var{x})}
## is ignored.
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

  if (! isnumeric (x) || ! isreal (x) || ! all (x(! isnan (x))(:) >= 0))
    error ("geomean: X must contain real nonnegative values.");
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
        x = x(isnan (x));
      endif

      if (any (x == 0))
        m = 0;
        return;
      endif

      m = exp (sum (log (x(:)), 1) ./ length (x(:)));

    else
      ## Find the first non-singleton dimension.
      (dim = find (szx != 1, 1)) || (dim = 1);
      n = szx(dim);
      if (omitnan)
        idx = isnan (x);
        n = sum (! idx, dim);
        x(idx) = 1;     # log (1) = 0
      endif

      m = exp (sum (log (x), dim) ./ n);
      m(m == -Inf) = 0; # handle zeros in X

    endif

  else

    ## Two numeric input arguments, dimensions given.  Note scalar is vector!
    vecdim = varargin{1};
    if (isempty (vecdim) || ! (isvector (vecdim) && all (vecdim > 0)) ...
          || any (rem (vecdim, 1)))
      error ("geomean: DIM must be a positive integer scalar or vector.");
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
          if (omitnan)
            nanx = isnan (x);
            n = sum (! nanx, vecdim);
            x(nanx) = 1;    # log (1) = 0
          endif

          m = exp (sum (log (x), vecdim) ./ n);
          m(m == -Inf) = 0; # handle zeros in X

        endif

      else
        vecdim = sort (vecdim);
        if (! all (diff (vecdim)))
           error (strcat (["geomean: VECDIM must contain non-repeating"], ...
                          [" positive integers."]));
        endif
        ## Ignore exceeding dimensions in VECDIM
        vecdim(find (vecdim > ndims (x))) = [];

        if (isempty (vecdim))
          m = x;
        else
          ## Move vecdims to dim 1.

          ## Calculate permutation vector
          remdims = 1 : ndx;    # All dimensions
          remdims(vecdim) = [];     # Delete dimensions specified by vecdim
          nremd = numel (remdims);

          ## If all dimensions are given, it is similar to all flag
          if (nremd == 0)
            x = x(:);

            if (omitnan)
              x = x(isnan (x));
            endif

            if (any (x == 0))
              m = 0;
              return;
            endif

            m = exp (sum (log (x(:)), 1) ./ length (x(:)));

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
              x(nanx) = 1;    # log (1) = 0
            else
              n = szx(1);
            endif

            m = exp (sum (log (x), 1) ./ n);
            m(m == -Inf) = 0; # handle zeros in X

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
%! assert (geomean (x), 0);
%! assert (geomean (y, 2), [0, 9.462942809849169, 14.65658770861967]', 1e-14);
%! assert (geomean (y, "all"), 0);
%! y(2,4) = NaN;
%! assert (geomean (y, 2), [0 NaN 14.65658770861967]', 1e-14);
%! assert (geomean (y', "omitnan"), [0 9.623207231679554 14.65658770861967], 1e-14);
%! z = y + 20;
%! assert (geomean (z, "all"), NaN);
%! m = [24.79790781765634 NaN 34.85638839503932];
%! assert (geomean (z'), m, 4e-14);
%! assert (geomean (z', "includenan"), m, 4e-14);
%! m(2) = 30.02181156156319;
%! assert (geomean (z', "omitnan"), m, 4e-14);
%! assert (geomean (z, 2, "omitnan"), m', 4e-14);

## Test dimension indexing with vecdim in n-dimensional arrays
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! assert (size (geomean (x, [3 2])), [10 1 1 3]);
%! assert (size (geomean (x, [1 2])), [1 1 6 3]);
%! assert (size (geomean (x, [1 2 4])), [1 1 6]);
%! assert (size (geomean (x, [1 4 3])), [1 40]);
%! assert (size (geomean (x, [1 2 3 4])), [1 1]);

## Test results with vecdim in n-dimensional arrays and "omitnan"
%!test
%! x = repmat ([1:20;6:25], [5 2 6 3]);
%! m = repmat ([8.304361203739333;14.3078118884256], [5,1,1,3]);
%! assert (geomean (x, [3 2]), m, 4e-13);
%! x(2,5,6,3) = NaN;
%! m(2,3) = NaN;
%! assert (geomean (x, [3 2]), m, 4e-13);
%! m(2,3) = 14.3292729579901;
%! assert (geomean (x, [3 2], "omitnan"), m, 4e-13);

## Test errors
%!error <geomean: X must contain real nonnegative values.> geomean ("char")
%!error <geomean: X must contain real nonnegative values.> geomean ([1 -1 3])
%!error <geomean: DIM must be a positive integer scalar or vector.> ...
%! geomean (repmat ([1:20;6:25], [5 2 6 3 5]), -1)
%!error <geomean: DIM must be a positive integer scalar or vector.> ...
%! geomean (repmat ([1:20;6:25], [5 2 6 3 5]), 0)
%!error <geomean: VECDIM must contain non-repeating positive integers.> ...
%! geomean (repmat ([1:20;6:25], [5 2 6 3 5]), [1 1])
