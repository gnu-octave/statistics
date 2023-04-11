## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} wblrnd (@var{lambda}, @var{k})
## @deftypefnx {statistics} {@var{r} =} wblrnd (@var{lambda}, @var{k}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} wblrnd (@var{lambda}, @var{k}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} wblrnd (@var{lambda}, @var{k}, [@var{sz}])
##
## Random arrays from the Weibull distribution.
##
## @code{@var{r} = wblrnd (@var{lambda}, @var{k})} returns an array of random
## numbers chosen from the Weibull distribution with parameters @var{lambda} and
## @var{k}.  The size of @var{r} is the common size of @var{lambda} and @var{k}.
##  A scalar input functions as a constant matrix of the same size as the other
## inputs.  Both parameters must be positive reals.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{wblcdf, wblinv, wblpdf, wblstat, wblplot}
## @end deftypefn

function r = wblrnd (lambda, k, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  if (! isscalar (lambda) || ! isscalar (k))
    [retval, lambda, k] = common_size (lambda, k);
    if (retval > 0)
      error ("wblrnd: SCALE and SHAPE must be of common size or scalars");
    endif
  endif

  if (iscomplex (lambda) || iscomplex (k))
    error ("wblrnd: SCALE and SHAPE must not be complex");
  endif

  if (nargin == 2)
    sz = size (lambda);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error ("wblrnd: dimension vector must be row vector of non-negative integers");
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("wblrnd: dimensions must be non-negative integers");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (lambda) && ! isequal (size (lambda), sz))
    error ("wblrnd: SCALE and SHAPE must be scalar or of size SZ");
  endif

  if (isa (lambda, "single") || isa (k, "single"))
    cls = "single";
  else
    cls = "double";
  endif

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


%!assert (size (wblrnd (1,2)), [1, 1])
%!assert (size (wblrnd (ones (2,1), 2)), [2, 1])
%!assert (size (wblrnd (ones (2,2), 2)), [2, 2])
%!assert (size (wblrnd (1, 2*ones (2,1))), [2, 1])
%!assert (size (wblrnd (1, 2*ones (2,2))), [2, 2])
%!assert (size (wblrnd (1, 2, 3)), [3, 3])
%!assert (size (wblrnd (1, 2, [4 1])), [4, 1])
%!assert (size (wblrnd (1, 2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (wblrnd (1, 2)), "double")
%!assert (class (wblrnd (single (1), 2)), "single")
%!assert (class (wblrnd (single ([1 1]), 2)), "single")
%!assert (class (wblrnd (1, single (2))), "single")
%!assert (class (wblrnd (1, single ([2 2]))), "single")

## Test input validation
%!error wblrnd ()
%!error wblrnd (1)
%!error wblrnd (ones (3), ones (2))
%!error wblrnd (ones (2), ones (3))
%!error wblrnd (i, 2)
%!error wblrnd (2, i)
%!error wblrnd (1,2, -1)
%!error wblrnd (1,2, ones (2))
%!error wblrnd (1, 2, [2 -1 2])
%!error wblrnd (1,2, 1, ones (2))
%!error wblrnd (1,2, 1, -1)
%!error wblrnd (ones (2,2), 2, 3)
%!error wblrnd (ones (2,2), 2, [3, 2])
%!error wblrnd (ones (2,2), 2, 2, 3)
