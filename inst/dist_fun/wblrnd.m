## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{r} =} wblrnd (@var{lambda}, @var{k})
## @deftypefnx {statistics} {@var{r} =} wblrnd (@var{lambda}, @var{k}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} wblrnd (@var{lambda}, @var{k}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} wblrnd (@var{lambda}, @var{k}, [@var{sz}])
##
## Random arrays from the Weibull distribution.
##
## @code{@var{r} = wblrnd (@var{lambda}, @var{k})} returns an array of random
## numbers chosen from the Weibull distribution with scale parameter
## @var{lambda} and shape parameter @var{k}.  The size of @var{r} is the common
## size of @var{lambda} and @var{k}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.  Both parameters must be
## positive reals.
##
## When called with a single size argument, @code{wblrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the Weibull distribution can be found at
## @url{https://en.wikipedia.org/wiki/Weibull_distribution}
##
## @seealso{wblcdf, wblinv, wblpdf, wblfit, wbllike, wblstat, wblplot}
## @end deftypefn

function r = wblrnd (lambda, k, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("wblrnd: function called with too few input arguments.");
  endif

  ## Check for common size of LAMBDA and K
  if (! isscalar (lambda) || ! isscalar (k))
    [retval, lambda, k] = common_size (lambda, k);
    if (retval > 0)
      error ("wblrnd: LAMBDA and K must be of common size or scalars.");
    endif
  endif

  ## Check for LAMBDA and K being reals
  if (iscomplex (lambda) || iscomplex (k))
    error ("wblrnd: LAMBDA and K must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (lambda);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif ((isrow (varargin{1}) || isempty (varargin{1})) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("wblrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("wblrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (lambda) && ! isequal (size (lambda), sz))
    error ("wblrnd: LAMBDA and K must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (lambda, "single") || isa (k, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Weibull distribution
  if (isscalar (lambda) && isscalar (k))
    if ((lambda > 0) && (lambda < Inf) && (k > 0) && (k < Inf))
      r = lambda * rande (sz, cls) .^ (1/k);
    else
      r = NaN (sz, cls);
    endif
  else
    r = lambda .* rande (sz, cls) .^ (1./k);

    is_nan = (lambda <= 0) | (lambda == Inf) | (k <= 0) | (k == Inf);
    r(is_nan) = NaN;
  endif

endfunction

## Test output
%!assert (size (wblrnd (1, 1)), [1 1])
%!assert (size (wblrnd (1, ones (2,1))), [2, 1])
%!assert (size (wblrnd (1, ones (2,2))), [2, 2])
%!assert (size (wblrnd (ones (2,1), 1)), [2, 1])
%!assert (size (wblrnd (ones (2,2), 1)), [2, 2])
%!assert (size (wblrnd (1, 1, 3)), [3, 3])
%!assert (size (wblrnd (1, 1, [4, 1])), [4, 1])
%!assert (size (wblrnd (1, 1, 4, 1)), [4, 1])
%!assert (size (wblrnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (wblrnd (1, 1, 0, 1)), [0, 1])
%!assert (size (wblrnd (1, 1, 1, 0)), [1, 0])
%!assert (size (wblrnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (wblrnd (1, 1)), "double")
%!assert (class (wblrnd (1, single (1))), "single")
%!assert (class (wblrnd (1, single ([1, 1]))), "single")
%!assert (class (wblrnd (single (1), 1)), "single")
%!assert (class (wblrnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<wblrnd: function called with too few input arguments.> wblrnd ()
%!error<wblrnd: function called with too few input arguments.> wblrnd (1)
%!error<wblrnd: LAMBDA and K must be of common size or scalars.> ...
%! wblrnd (ones (3), ones (2))
%!error<wblrnd: LAMBDA and K must be of common size or scalars.> ...
%! wblrnd (ones (2), ones (3))
%!error<wblrnd: LAMBDA and K must not be complex.> wblrnd (i, 2, 3)
%!error<wblrnd: LAMBDA and K must not be complex.> wblrnd (1, i, 3)
%!error<wblrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! wblrnd (1, 2, -1)
%!error<wblrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! wblrnd (1, 2, 1.2)
%!error<wblrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! wblrnd (1, 2, ones (2))
%!error<wblrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! wblrnd (1, 2, [2 -1 2])
%!error<wblrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! wblrnd (1, 2, [2 0 2.5])
%!error<wblrnd: dimensions must be non-negative integers.> ...
%! wblrnd (1, 2, 2, -1, 5)
%!error<wblrnd: dimensions must be non-negative integers.> ...
%! wblrnd (1, 2, 2, 1.5, 5)
%!error<wblrnd: LAMBDA and K must be scalar or of size SZ.> ...
%! wblrnd (2, ones (2), 3)
%!error<wblrnd: LAMBDA and K must be scalar or of size SZ.> ...
%! wblrnd (2, ones (2), [3, 2])
%!error<wblrnd: LAMBDA and K must be scalar or of size SZ.> ...
%! wblrnd (2, ones (2), 3, 2)
