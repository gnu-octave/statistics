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
## @deftypefn  {statistics} {@var{r} =} gumbelrnd (@var{mu}, @var{beta})
## @deftypefnx {statistics} {@var{r} =} gumbelrnd (@var{mu}, @var{beta}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} gumbelrnd (@var{mu}, @var{beta}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} gumbelrnd (@var{mu}, @var{beta}, [@var{sz}])
##
## Random arrays from the Gumbel distribution.
##
## @code{@var{r} = gumbelrnd (@var{mu}, @var{beta})} returns an array of random
## numbers chosen from the Gumbel distribution (also known as the extreme value
## or the type I generalized extreme value distribution) with location
## parameter @var{mu} and scale parameter @var{beta}.  The size of @var{r} is
## the common size of @var{mu} and @var{beta}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## When called with a single size argument, @code{gumbelrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## The Gumbel distribution is used to model the distribution of the maximum (or
## the minimum) of a number of samples of various distributions.  This version
## is suitable for modeling maxima.  For modeling minima, use the alternative
## extreme value iCDF, @code{evinv}.
##
## Further information about the Gumbel distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gumbel_distribution}
##
## @seealso{gumbelcdf, gumbelinv, gumbelpdf, gumbelfit, gumbellike, gumbelstat,
## evrnd}
## @end deftypefn

function r = gumbelrnd (mu, beta, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("gumbelrnd: function called with too few input arguments.");
  endif

  ## Check for common size of MU and BETA
  if (! isscalar (mu) || ! isscalar (beta))
    [retval, mu, beta] = common_size (mu, beta);
    if (retval > 0)
      error ("gumbelrnd: MU and BETA must be of common size or scalars.");
    endif
  endif

  ## Check for MU and BETA being reals
  if (iscomplex (mu) || iscomplex (beta))
    error ("gumbelrnd: MU and BETA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (mu);
  elseif (nargin == 3)
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
      error (strcat ("gumbelrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("gumbelrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  ## Use 'size (ones (sz))' to ignore any trailing singleton dimensions in SZ
  if (! isscalar (mu) && ! isequal (size (mu), size (ones (sz))))
    error ("gumbelrnd: MU and BETA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single") || isa (beta, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Return NaNs for out of range values of BETA
  beta(beta < 0) = NaN;

  ## Generate uniform random values, and apply the extreme value inverse CDF.
  r = -log (-log (rand (sz, cls))) .* beta + mu;

endfunction

## Test output
%!assert (size (gumbelrnd (1, 1)), [1 1])
%!assert (size (gumbelrnd (1, ones (2,1))), [2, 1])
%!assert (size (gumbelrnd (1, ones (2,2))), [2, 2])
%!assert (size (gumbelrnd (ones (2,1), 1)), [2, 1])
%!assert (size (gumbelrnd (ones (2,2), 1)), [2, 2])
%!assert (size (gumbelrnd (1, 1, 3)), [3, 3])
%!assert (size (gumbelrnd (1, 1, [4, 1])), [4, 1])
%!assert (size (gumbelrnd (1, 1, 4, 1)), [4, 1])
%!assert (size (gumbelrnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (gumbelrnd (1, 1, 0, 1)), [0, 1])
%!assert (size (gumbelrnd (1, 1, 1, 0)), [1, 0])
%!assert (size (gumbelrnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])
%!assert (size (gumbelrnd (1, 1, [])), [0, 0])
%!assert (size (gumbelrnd (1, 1, [2, 0, 2, 1])), [2, 0, 2])

## Test class of input preserved
%!assert (class (gumbelrnd (1, 1)), "double")
%!assert (class (gumbelrnd (1, single (1))), "single")
%!assert (class (gumbelrnd (1, single ([1, 1]))), "single")
%!assert (class (gumbelrnd (single (1), 1)), "single")
%!assert (class (gumbelrnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<gumbelrnd: function called with too few input arguments.> gumbelrnd ()
%!error<gumbelrnd: function called with too few input arguments.> gumbelrnd (1)
%!error<gumbelrnd: MU and BETA must be of common size or scalars.> ...
%! gumbelrnd (ones (3), ones (2))
%!error<gumbelrnd: MU and BETA must be of common size or scalars.> ...
%! gumbelrnd (ones (2), ones (3))
%!error<gumbelrnd: MU and BETA must not be complex.> gumbelrnd (i, 2, 3)
%!error<gumbelrnd: MU and BETA must not be complex.> gumbelrnd (1, i, 3)
%!error<gumbelrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gumbelrnd (1, 2, -1)
%!error<gumbelrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gumbelrnd (1, 2, 1.2)
%!error<gumbelrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gumbelrnd (1, 2, ones (2))
%!error<gumbelrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gumbelrnd (1, 2, [2 -1 2])
%!error<gumbelrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gumbelrnd (1, 2, [2 0 2.5])
%!error<gumbelrnd: dimensions must be non-negative integers.> ...
%! gumbelrnd (1, 2, 2, -1, 5)
%!error<gumbelrnd: dimensions must be non-negative integers.> ...
%! gumbelrnd (1, 2, 2, 1.5, 5)
%!error<gumbelrnd: MU and BETA must be scalars or of size SZ.> ...
%! gumbelrnd (2, ones (2), 3)
%!error<gumbelrnd: MU and BETA must be scalars or of size SZ.> ...
%! gumbelrnd (2, ones (2), [3, 2])
%!error<gumbelrnd: MU and BETA must be scalars or of size SZ.> ...
%! gumbelrnd (2, ones (2), 3, 2)
