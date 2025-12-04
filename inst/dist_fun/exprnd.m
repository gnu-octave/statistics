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
## @deftypefn  {statistics} {@var{r} =} exprnd (@var{mu})
## @deftypefnx {statistics} {@var{r} =} exprnd (@var{mu}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} exprnd (@var{mu}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} exprnd (@var{mu}, [@var{sz}])
##
## Random arrays from the exponential distribution.
##
## @code{@var{r} = exprnd (@var{mu})} returns an array of random numbers chosen
## from the exponential distribution with mean parameter @var{mu}.  The size of
## @var{r} is the size of @var{mu}.
##
## When called with a single size argument, @code{exprnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## A common alternative parameterization of the exponential distribution is to
## use the parameter @math{λ} defined as the mean number of events in an
## interval as opposed to the parameter @math{μ}, which is the mean wait time
## for an event to occur. @math{λ} and @math{μ} are reciprocals,
## i.e. @math{μ = 1 / λ}.
##
## Further information about the exponential distribution can be found at
## @url{https://en.wikipedia.org/wiki/Exponential_distribution}
##
## @seealso{expcdf, expinv, exppdf, expfit, explike, expstat}
## @end deftypefn

function r = exprnd (mu, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("exprnd: function called with too few input arguments.");
  endif

  ## Check for MU being real
  if (iscomplex (mu))
    error ("exprnd: MU must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 1)
    sz = size (mu);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("exprnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 2)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("exprnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("exprnd: MU must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from exponential distribution
  if (isscalar (mu))
    if ((mu > 0) && (mu < Inf))
      r = rande (sz, cls) * mu;
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    k = (mu > 0) & (mu < Inf);
    r(k) = rande (sum (k(:)), 1, cls) .* mu(k)(:);
  endif

endfunction

## Test output
%!assert (size (exprnd (2)), [1, 1])
%!assert (size (exprnd (ones (2,1))), [2, 1])
%!assert (size (exprnd (ones (2,2))), [2, 2])
%!assert (size (exprnd (1, 3)), [3, 3])
%!assert (size (exprnd (1, [4 1])), [4, 1])
%!assert (size (exprnd (1, 4, 1)), [4, 1])
%!assert (size (exprnd (1, 4, 1)), [4, 1])
%!assert (size (exprnd (1, 4, 1, 5)), [4, 1, 5])
%!assert (size (exprnd (1, 0, 1)), [0, 1])
%!assert (size (exprnd (1, 1, 0)), [1, 0])
%!assert (size (exprnd (1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (exprnd (2)), "double")
%!assert (class (exprnd (single (2))), "single")
%!assert (class (exprnd (single ([2 2]))), "single")

## Test input validation
%!error<exprnd: function called with too few input arguments.> exprnd ()
%!error<exprnd: MU must not be complex.> exprnd (i)
%!error<exprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! exprnd (1, -1)
%!error<exprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! exprnd (1, 1.2)
%!error<exprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! exprnd (1, ones (2))
%!error<exprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! exprnd (1, [2 -1 2])
%!error<exprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! exprnd (1, [2 0 2.5])
%!error<exprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! exprnd (ones (2), ones (2))
%!error<exprnd: dimensions must be non-negative integers.> ...
%! exprnd (1, 2, -1, 5)
%!error<exprnd: dimensions must be non-negative integers.> ...
%! exprnd (1, 2, 1.5, 5)
%!error<exprnd: MU must be scalar or of size SZ.> exprnd (ones (2,2), 3)
%!error<exprnd: MU must be scalar or of size SZ.> exprnd (ones (2,2), [3, 2])
%!error<exprnd: MU must be scalar or of size SZ.> exprnd (ones (2,2), 2, 3)
