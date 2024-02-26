## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{r} =} gamrnd (@var{a}, @var{b})
## @deftypefnx {statistics} {@var{r} =} gamrnd (@var{a}, @var{b}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} gamrnd (@var{a}, @var{b}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} gamrnd (@var{a}, @var{b}, [@var{sz}])
##
## Random arrays from the Gamma distribution.
##
## @code{@var{r} = gamrnd (@var{a}, @var{b})} returns an array of random
## numbers chosen from the Gamma distribution with shape parameter @var{a} and
## scale parameter @var{b}.  The size of @var{r} is the common size of
## @var{a} and @var{b}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## When called with a single size argument, @code{gamrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## OCTAVE/MATLAB use the alternative parameterization given by the pair
## @math{α, β}, i.e. shape @var{a} and scale @var{b}.  In Wikipedia, the two
## common parameterizations use the pairs @math{k, θ}, as shape and scale, and
## @math{α, β}, as shape and rate, respectively.  The parameter names @var{a}
## and @var{b} used here (for MATLAB compatibility) correspond to the parameter
## notation @math{k, θ} instead of the @math{α, β} as reported in Wikipedia.
##
## Further information about the Gamma distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gamma_distribution}
##
## @seealso{gamcdf, gaminv, gampdf, gamfit, gamlike, gamstat}
## @end deftypefn

function r = gamrnd (a, b, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("gamrnd: function called with too few input arguments.");
  endif

  ## Check for common size of A and B
  if (! isscalar (a) || ! isscalar (b))
    [retval, a, b] = common_size (a, b);
    if (retval > 0)
      error ("gamrnd: A and B must be of common size or scalars.");
    endif
  endif

  ## Check for A and B being reals
  if (iscomplex (a) || iscomplex (b))
    error ("gamrnd: A and B must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (a);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["gamrnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("gamrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (a) && ! isequal (size (a), sz))
    error ("gamrnd: A and B must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (a, "single") || isa (b, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Gamma distribution
  if (isscalar (a) && isscalar (b))
    if ((a > 0) && (a < Inf) && (b > 0) && (b < Inf))
      r = b * randg (a, sz, cls);
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    valid = (a > 0) & (a < Inf) & (b > 0) & (b < Inf);
    r(valid) = b(valid) .* randg (a(valid), cls);
  endif

endfunction

## Test output
%!assert (size (gamrnd (1, 1)), [1 1])
%!assert (size (gamrnd (1, ones (2,1))), [2, 1])
%!assert (size (gamrnd (1, ones (2,2))), [2, 2])
%!assert (size (gamrnd (ones (2,1), 1)), [2, 1])
%!assert (size (gamrnd (ones (2,2), 1)), [2, 2])
%!assert (size (gamrnd (1, 1, 3)), [3, 3])
%!assert (size (gamrnd (1, 1, [4, 1])), [4, 1])
%!assert (size (gamrnd (1, 1, 4, 1)), [4, 1])
%!assert (size (gamrnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (gamrnd (1, 1, 0, 1)), [0, 1])
%!assert (size (gamrnd (1, 1, 1, 0)), [1, 0])
%!assert (size (gamrnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (gamrnd (1, 1)), "double")
%!assert (class (gamrnd (1, single (1))), "single")
%!assert (class (gamrnd (1, single ([1, 1]))), "single")
%!assert (class (gamrnd (single (1), 1)), "single")
%!assert (class (gamrnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<gamrnd: function called with too few input arguments.> gamrnd ()
%!error<gamrnd: function called with too few input arguments.> gamrnd (1)
%!error<gamrnd: A and B must be of common size or scalars.> ...
%! gamrnd (ones (3), ones (2))
%!error<gamrnd: A and B must be of common size or scalars.> ...
%! gamrnd (ones (2), ones (3))
%!error<gamrnd: A and B must not be complex.> gamrnd (i, 2, 3)
%!error<gamrnd: A and B must not be complex.> gamrnd (1, i, 3)
%!error<gamrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gamrnd (1, 2, -1)
%!error<gamrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gamrnd (1, 2, 1.2)
%!error<gamrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gamrnd (1, 2, ones (2))
%!error<gamrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gamrnd (1, 2, [2 -1 2])
%!error<gamrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gamrnd (1, 2, [2 0 2.5])
%!error<gamrnd: dimensions must be non-negative integers.> ...
%! gamrnd (1, 2, 2, -1, 5)
%!error<gamrnd: dimensions must be non-negative integers.> ...
%! gamrnd (1, 2, 2, 1.5, 5)
%!error<gamrnd: A and B must be scalars or of size SZ.> ...
%! gamrnd (2, ones (2), 3)
%!error<gamrnd: A and B must be scalars or of size SZ.> ...
%! gamrnd (2, ones (2), [3, 2])
%!error<gamrnd: A and B must be scalars or of size SZ.> ...
%! gamrnd (2, ones (2), 3, 2)
