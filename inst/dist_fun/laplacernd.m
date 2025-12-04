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
## @deftypefn  {statistics} {@var{r} =} laplacernd (@var{mu}, @var{beta})
## @deftypefnx {statistics} {@var{r} =} laplacernd (@var{mu}, @var{beta}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} laplacernd (@var{mu}, @var{beta}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} laplacernd (@var{mu}, @var{beta}, [@var{sz}])
##
## Random arrays from the Laplace distribution.
##
## @code{@var{r} = laplacernd (@var{mu}, @var{beta})} returns an array of
## random numbers chosen from the Laplace distribution with location parameter
## @var{mu} and scale parameter @var{beta}.  The size of @var{r} is the common
## size of @var{mu} and @var{beta}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Both parameters must be reals and @qcode{@var{beta} > 0}.
## For @qcode{@var{beta} <= 0}, @qcode{NaN} is returned.
##
## When called with a single size argument, @code{laplacernd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the Laplace distribution can be found at
## @url{https://en.wikipedia.org/wiki/Laplace_distribution}
##
## @seealso{laplacecdf, laplaceinv, laplacernd}
## @end deftypefn

function r = laplacernd (mu, beta, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("laplacernd: function called with too few input arguments.");
  endif

  ## Check for common size of MU, and BETA
  if (! isscalar (mu) || ! isscalar (beta))
    [retval, mu, beta] = common_size (mu, beta);
    if (retval > 0)
      error ("laplacernd: MU and BETA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and BETA being reals
  if (iscomplex (mu) || iscomplex (beta))
    error ("laplacernd: MU and BETA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (mu);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("laplacernd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("laplacernd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("laplacernd: MU and BETA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single") || isa (beta, "single"))
    is_type = "single";
  else
    is_type = "double";
  endif

  ## Generate random sample from Laplace distribution
  tmp = rand (sz, is_type);
  r = ((tmp < 1/2) .* log (2 * tmp) - ...
       (tmp > 1/2) .* log (2 * (1 - tmp))) .* beta + mu;

  ## Force output to NaN for invalid parameter BETA <= 0
  k = (beta <= 0);
  r(k) = NaN;

endfunction

## Test output
%!assert (size (laplacernd (1, 1)), [1 1])
%!assert (size (laplacernd (1, ones (2,1))), [2, 1])
%!assert (size (laplacernd (1, ones (2,2))), [2, 2])
%!assert (size (laplacernd (ones (2,1), 1)), [2, 1])
%!assert (size (laplacernd (ones (2,2), 1)), [2, 2])
%!assert (size (laplacernd (1, 1, 3)), [3, 3])
%!assert (size (laplacernd (1, 1, [4, 1])), [4, 1])
%!assert (size (laplacernd (1, 1, 4, 1)), [4, 1])
%!assert (size (laplacernd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (laplacernd (1, 1, 0, 1)), [0, 1])
%!assert (size (laplacernd (1, 1, 1, 0)), [1, 0])
%!assert (size (laplacernd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (laplacernd (1, 1)), "double")
%!assert (class (laplacernd (1, single (1))), "single")
%!assert (class (laplacernd (1, single ([1, 1]))), "single")
%!assert (class (laplacernd (single (1), 1)), "single")
%!assert (class (laplacernd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<laplacernd: function called with too few input arguments.> laplacernd ()
%!error<laplacernd: function called with too few input arguments.> laplacernd (1)
%!error<laplacernd: MU and BETA must be of common size or scalars.> ...
%! laplacernd (ones (3), ones (2))
%!error<laplacernd: MU and BETA must be of common size or scalars.> ...
%! laplacernd (ones (2), ones (3))
%!error<laplacernd: MU and BETA must not be complex.> laplacernd (i, 2, 3)
%!error<laplacernd: MU and BETA must not be complex.> laplacernd (1, i, 3)
%!error<laplacernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! laplacernd (1, 2, -1)
%!error<laplacernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! laplacernd (1, 2, 1.2)
%!error<laplacernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! laplacernd (1, 2, ones (2))
%!error<laplacernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! laplacernd (1, 2, [2 -1 2])
%!error<laplacernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! laplacernd (1, 2, [2 0 2.5])
%!error<laplacernd: dimensions must be non-negative integers.> ...
%! laplacernd (1, 2, 2, -1, 5)
%!error<laplacernd: dimensions must be non-negative integers.> ...
%! laplacernd (1, 2, 2, 1.5, 5)
%!error<laplacernd: MU and BETA must be scalars or of size SZ.> ...
%! laplacernd (2, ones (2), 3)
%!error<laplacernd: MU and BETA must be scalars or of size SZ.> ...
%! laplacernd (2, ones (2), [3, 2])
%!error<laplacernd: MU and BETA must be scalars or of size SZ.> ...
%! laplacernd (2, ones (2), 3, 2)
