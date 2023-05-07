## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR LAMBDA PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} burrrnd (@var{lambda}, @var{c}, @var{k})
## @deftypefnx {statistics} {@var{r} =} burrrnd (@var{lambda}, @var{c}, @var{k}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} burrrnd (@var{lambda}, @var{c}, @var{k}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} burrrnd (@var{lambda}, @var{c}, @var{k}, [@var{sz}])
##
## Random arrays from the Burr type XII distribution.
##
## @code{@var{r} = burrrnd (@var{lambda}, @var{c}, @var{k})} returns an array of
## random numbers chosen from the Burr type XII distribution with scale
## parameter @var{lambda}, first shape parameter @var{c}, and second shape
## parameter@var{c}, and @var{k}.  The size of @var{r} is the common size of
## @var{lambda}, @var{c}, and @var{k}.  LAMBDA scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## When called with a single size argument, @code{burrrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the Burr distribution can be found at
## @url{https://en.wikipedia.org/wiki/Burr_distribution}
##
## @seealso{burrcdf, burrinv, burrpdf, burrfit, burrlike}
## @end deftypefn

function r = burrrnd (lambda, c, k, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("burrrnd: function called with too few input arguments.");
  endif

  ## Check for common size of LAMBDA, C, and K
  if (! isscalar (lambda) || ! isscalar (c) || ! isscalar (k))
    [retval, lambda, c, k] = common_size (lambda, c, k);
    if (retval > 0)
      error ("burrrnd: LAMBDA, C, and K must be of common size or scalars.");
    endif
  endif

  ## Check for LAMBDA, C, and K being reals
  if (iscomplex (lambda) || iscomplex (c) || iscomplex (k))
    error ("burrrnd: LAMBDA, C, and K must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 3)
    sz = size (lambda);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["burrrnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 4)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("burrrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (lambda) && ! isequal (size (lambda), sz))
    error ("burrrnd: LAMBDA, C, and K must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (lambda, "single") || isa (c, "single") || isa (k, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Burr type XII distribution
  lambda(lambda <= 0) = NaN;
  c(c <= 0) = NaN;
  k(k <= 0) = NaN;
  r = lambda .* (((1 - rand (sz, cls)) .^ (-(1./k))) - 1) .^ (1./c);

endfunction

## Test output
%!assert (size (burrrnd (1, 1, 1)), [1 1])
%!assert (size (burrrnd (ones (2,1), 1, 1)), [2, 1])
%!assert (size (burrrnd (ones (2,2), 1, 1)), [2, 2])
%!assert (size (burrrnd (1, ones (2,1), 1)), [2, 1])
%!assert (size (burrrnd (1, ones (2,2), 1)), [2, 2])
%!assert (size (burrrnd (1, 1, ones (2,1))), [2, 1])
%!assert (size (burrrnd (1, 1, ones (2,2))), [2, 2])
%!assert (size (burrrnd (1, 1, 1, 3)), [3, 3])
%!assert (size (burrrnd (1, 1, 1, [4 1])), [4, 1])
%!assert (size (burrrnd (1, 1, 1, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (burrrnd (1,1,1)), "double")
%!assert (class (burrrnd (single (1),1,1)), "single")
%!assert (class (burrrnd (single ([1 1]),1,1)), "single")
%!assert (class (burrrnd (1,single (1),1)), "single")
%!assert (class (burrrnd (1,single ([1 1]),1)), "single")
%!assert (class (burrrnd (1,1,single (1))), "single")
%!assert (class (burrrnd (1,1,single ([1 1]))), "single")

## Test input validation
%!error<burrrnd: function called with too few input arguments.> burrrnd ()
%!error<burrrnd: function called with too few input arguments.> burrrnd (1)
%!error<burrrnd: function called with too few input arguments.> burrrnd (1, 2)
%!error<burrrnd: LAMBDA, C, and K must be of common size or scalars.> ...
%! burrrnd (ones (3), ones (2), ones (2))
%!error<burrrnd: LAMBDA, C, and K must be of common size or scalars.> ...
%! burrrnd (ones (2), ones (3), ones (2))
%!error<burrrnd: LAMBDA, C, and K must be of common size or scalars.> ...
%! burrrnd (ones (2), ones (2), ones (3))
%!error<burrrnd: LAMBDA, C, and K must not be complex.> burrrnd (i, 2, 3)
%!error<burrrnd: LAMBDA, C, and K must not be complex.> burrrnd (1, i, 3)
%!error<burrrnd: LAMBDA, C, and K must not be complex.> burrrnd (1, 2, i)
%!error<burrrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! burrrnd (1, 2, 3, -1)
%!error<burrrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! burrrnd (1, 2, 3, 1.2)
%!error<burrrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! burrrnd (1, 2, 3, ones (2))
%!error<burrrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! burrrnd (1, 2, 3, [2 -1 2])
%!error<burrrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! burrrnd (1, 2, 3, [2 0 2.5])
%!error<burrrnd: dimensions must be non-negative integers.> ...
%! burrrnd (1, 2, 3, 2, -1, 5)
%!error<burrrnd: dimensions must be non-negative integers.> ...
%! burrrnd (1, 2, 3, 2, 1.5, 5)
%!error<burrrnd: LAMBDA, C, and K must be scalar or of size SZ.> ...
%! burrrnd (2, ones (2), 2, 3)
%!error<burrrnd: LAMBDA, C, and K must be scalar or of size SZ.> ...
%! burrrnd (2, ones (2), 2, [3, 2])
%!error<burrrnd: LAMBDA, C, and K must be scalar or of size SZ.> ...
%! burrrnd (2, ones (2), 2, 3, 2)
