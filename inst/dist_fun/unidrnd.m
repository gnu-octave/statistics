## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 2005-2016 John W. Eaton
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
## @deftypefn  {statistics} {@var{r} =} unidrnd (@var{N})
## @deftypefnx {statistics} {@var{r} =} unidrnd (@var{N}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} unidrnd (@var{N}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} unidrnd (@var{N}, [@var{sz}])
##
## Random arrays from the discrete uniform distribution.
##
## @code{@var{r} = unidrnd (@var{N})} returns an array of random numbers chosen
## from the discrete uniform distribution with parameter @var{N}, which
## corresponds to the maximum observable value.  @code{unidrnd} assumes the
## integer values in the range @math{[1,N]} with equal probability.  The size of
## @var{r} is the size of @var{N}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## The maximum observable values in @var{N} must be positive integers, otherwise
## @qcode{NaN} is returned.
##
## When called with a single size argument, @code{unidrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Warning: The underlying implementation uses the double class and will only
## be accurate for @var{N} < @code{flintmax} (@w{@math{2^{53}}} on
## IEEE 754 compatible systems).
##
## Further information about the discrete uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Discrete_uniform_distribution}
##
## @seealso{unidcdf, unidinv, unidrnd, unidfit, unidstat}
## @end deftypefn

function r = unidrnd (N, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("unidrnd: function called with too few input arguments.");
  endif

  ## Check for N being real
  if (iscomplex (N))
    error ("unidrnd: N must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 1)
    sz = size (N);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif ((isrow (varargin{1}) || isempty (varargin{1})) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("unidrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 2)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("unidrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (N) && ! isequal (size (N), sz))
    error ("unidrnd: N must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (N, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (N))
    if (N > 0 && N == fix (N))
      r = ceil (rand (sz, cls) * N);
    else
      r = NaN (sz, cls);
    endif
  else
    r = ceil (rand (sz, cls) .* N);

    k = ! (N > 0 & N == fix (N));
    r(k) = NaN;
  endif

endfunction

## Test output
%!assert (size (unidrnd (2)), [1, 1])
%!assert (size (unidrnd (ones (2,1))), [2, 1])
%!assert (size (unidrnd (ones (2,2))), [2, 2])
%!assert (size (unidrnd (1, 3)), [3, 3])
%!assert (size (unidrnd (1, [4 1])), [4, 1])
%!assert (size (unidrnd (1, 4, 1)), [4, 1])
%!assert (size (unidrnd (1, 4, 1)), [4, 1])
%!assert (size (unidrnd (1, 4, 1, 5)), [4, 1, 5])
%!assert (size (unidrnd (1, 0, 1)), [0, 1])
%!assert (size (unidrnd (1, 1, 0)), [1, 0])
%!assert (size (unidrnd (1, 1, 2, 0, 5)), [1, 2, 0, 5])
%!assert (unidrnd (0, 1, 1), NaN)
%!assert (unidrnd ([0, 0, 0], [1, 3]), [NaN, NaN, NaN])

## Test class of input preserved
%!assert (class (unidrnd (2)), "double")
%!assert (class (unidrnd (single (2))), "single")
%!assert (class (unidrnd (single ([2 2]))), "single")

## Test input validation
%!error<unidrnd: function called with too few input arguments.> unidrnd ()
%!error<unidrnd: N must not be complex.> unidrnd (i)
%!error<unidrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! unidrnd (1, -1)
%!error<unidrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! unidrnd (1, 1.2)
%!error<unidrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! unidrnd (1, ones (2))
%!error<unidrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! unidrnd (1, [2 -1 2])
%!error<unidrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! unidrnd (1, [2 0 2.5])
%!error<unidrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! unidrnd (ones (2), ones (2))
%!error<unidrnd: dimensions must be non-negative integers.> ...
%! unidrnd (1, 2, -1, 5)
%!error<unidrnd: dimensions must be non-negative integers.> ...
%! unidrnd (1, 2, 1.5, 5)
%!error<unidrnd: N must be scalar or of size SZ.> unidrnd (ones (2,2), 3)
%!error<unidrnd: N must be scalar or of size SZ.> unidrnd (ones (2,2), [3, 2])
%!error<unidrnd: N must be scalar or of size SZ.> unidrnd (ones (2,2), 2, 3)
