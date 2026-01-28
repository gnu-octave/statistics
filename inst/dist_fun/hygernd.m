## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1997-2016 Kurt Hornik
## Copyright (C) 2022 Nicholas R. Jankowski
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
## @deftypefn  {statistics} {@var{r} =} hygernd (@var{m}, @var{k}, @var{n})
## @deftypefnx {statistics} {@var{r} =} hygernd (@var{m}, @var{k}, @var{n}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} hygernd (@var{m}, @var{k}, @var{n}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} hygernd (@var{m}, @var{k}, @var{n}, [@var{sz}])
##
## Random arrays from the hypergeometric distribution.
##
## @code{@var{r} = hygernd ((@var{m}, @var{k}, @var{n}} returns an array of
## random numbers chosen from the hypergeometric distribution with parameters
## @var{m}, @var{k}, and @var{n}.  The size of @var{r} is the common size of
## @var{m}, @var{k}, and @var{n}.  A scalar input functions as a constant matrix
## of the same size as the other inputs.
##
## The parameters @var{m}, @var{k}, and @var{n} must be positive integers
## with @var{k} and @var{n} not greater than @var{m}.
##
## When called with a single size argument, @code{hygernd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the hypergeometric distribution can be found at
## @url{https://en.wikipedia.org/wiki/Hypergeometric_distribution}
##
## @seealso{hygecdf, hygeinv, hygepdf, hygestat}
## @end deftypefn

function r = hygernd (m, k, n, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("hygernd: function called with too few input arguments.");
  endif

  ## Check for common size of T, M, and N
  if (! isscalar (m) || ! isscalar (k) || ! isscalar (n))
    [retval, m, k, n] = common_size (m, k, n);
    if (retval > 0)
      error ("hygernd: T, M, and N must be of common size or scalars.");
    endif
  endif

  ## Check for T, M, and N being reals
  if (iscomplex (m) || iscomplex (k) || iscomplex (n))
    error ("hygernd: T, M, and N must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 3)
    sz = size (m);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0)
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif (isempty (varargin{1}))
      r = [];
      return;
    else
      error (strcat ("hygernd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 4)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("hygernd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (m) && ! isequal (size (m), sz))
    error ("hygernd: T, M, and N must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (m, "single") || isa (k, "single") || isa (n, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ok = ((m >= 0) & (k >= 0) & (n > 0) & (k <= m) & (n <= m) &
        (m == fix (m)) & (k == fix (k)) & (n == fix (n)));

  ## Generate random sample from the hypergeometric distribution
  if (isscalar (m))
    if (ok)
      v = 0:n;
      p = hygepdf (v, m, k, n);
      r = v(lookup (cumsum (p(1:end-1)) / sum (p), rand (sz)) + 1);
      r = reshape (r, sz);
      if (strcmp (cls, "single"))
        r = single (r);
      endif
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);
    n = n(ok);
    num_n = numel (n);
    v = 0 : max (n(:));
    p = cumsum (hygepdf (v, m(ok), k(ok), n, "vectorexpand"), 2);

    ## Manual row-wise vectorization of lookup, which returns index of element
    ## less than or equal to test value, zero if test value is less than lowest
    ## number, and max index if greater than highest number.

    end_locs = sub2ind (size (p), [1 : num_n]', n(:) + 1);
    p = (p ./ p(end_locs)) - rand (num_n, 1);
    p(p>=0) = NaN;  # NaN values ignored by max
    [p_match, p_match_idx] = max (p, [], 2);
    p_match_idx(isnan(p_match)) = 0; # rand < min(p) gives NaN, reset to 0
    r(ok) = v(p_match_idx + 1);
  endif

endfunction

## Test output
%!assert (size (hygernd (4, 2, 2)), [1, 1])
%!assert (size (hygernd (4 * ones (2, 1), 2,2)), [2, 1])
%!assert (size (hygernd (4 * ones (2, 2), 2,2)), [2, 2])
%!assert (size (hygernd (4, 2 * ones (2, 1), 2)), [2, 1])
%!assert (size (hygernd (4, 2 * ones (2, 2), 2)), [2, 2])
%!assert (size (hygernd (4, 2, 2 * ones (2, 1))), [2, 1])
%!assert (size (hygernd (4, 2, 2 * ones (2, 2))), [2, 2])
%!assert (size (hygernd (4, 2, 2, 3)), [3, 3])
%!assert (size (hygernd (4, 2, 2, [4, 1])), [4, 1])
%!assert (size (hygernd (4, 2, 2, 4, 1)), [4, 1])
%!assert (size (hygernd (4, 2, 2, [])), [0, 0])
%!assert (size (hygernd (4, 2, 2, [2, 0, 2, 1])), [2, 0, 2])

## Test class of input preserved
%!assert (class (hygernd (4, 2, 2)), "double")
%!assert (class (hygernd (single (4), 2, 2)), "single")
%!assert (class (hygernd (single ([4, 4]), 2, 2)), "single")
%!assert (class (hygernd (4, single (2), 2)), "single")
%!assert (class (hygernd (4, single ([2, 2]),2)), "single")
%!assert (class (hygernd (4, 2, single (2))), "single")
%!assert (class (hygernd (4, 2, single ([2, 2]))), "single")

## Test input validation
%!error<hygernd: function called with too few input arguments.> hygernd ()
%!error<hygernd: function called with too few input arguments.> hygernd (1)
%!error<hygernd: function called with too few input arguments.> hygernd (1, 2)
%!error<hygernd: T, M, and N must be of common size or scalars.> ...
%! hygernd (ones (3), ones (2), ones (2))
%!error<hygernd: T, M, and N must be of common size or scalars.> ...
%! hygernd (ones (2), ones (3), ones (2))
%!error<hygernd: T, M, and N must be of common size or scalars.> ...
%! hygernd (ones (2), ones (2), ones (3))
%!error<hygernd: T, M, and N must not be complex.> hygernd (i, 2, 3)
%!error<hygernd: T, M, and N must not be complex.> hygernd (1, i, 3)
%!error<hygernd: T, M, and N must not be complex.> hygernd (1, 2, i)
%!error<hygernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hygernd (1, 2, 3, -1)
%!error<hygernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hygernd (1, 2, 3, 1.2)
%!error<hygernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hygernd (1, 2, 3, ones (2))
%!error<hygernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hygernd (1, 2, 3, [2 -1 2])
%!error<hygernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hygernd (1, 2, 3, [2 0 2.5])
%!error<hygernd: dimensions must be non-negative integers.> ...
%! hygernd (1, 2, 3, 2, -1, 5)
%!error<hygernd: dimensions must be non-negative integers.> ...
%! hygernd (1, 2, 3, 2, 1.5, 5)
%!error<hygernd: T, M, and N must be scalars or of size SZ.> ...
%! hygernd (2, ones (2), 2, 3)
%!error<hygernd: T, M, and N must be scalars or of size SZ.> ...
%! hygernd (2, ones (2), 2, [3, 2])
%!error<hygernd: T, M, and N must be scalars or of size SZ.> ...
%! hygernd (2, ones (2), 2, 3, 2)
