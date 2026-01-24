## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{r} =} ricernd (@var{s}, @var{sigma})
## @deftypefnx {statistics} {@var{r} =} ricernd (@var{s}, @var{sigma}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} ricernd (@var{s}, @var{sigma}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} ricernd (@var{s}, @var{sigma}, [@var{sz}])
##
## Random arrays from the Rician distribution.
##
## @code{@var{r} = ricernd (@var{s}, @var{sigma})} returns an array of random
## numbers chosen from the Rician distribution with noncentrality parameter
## @var{s} and scale parameter @var{sigma}.  The size of @var{r} is the common
## size of @var{s} and @var{sigma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## When called with a single size argument, @code{ricernd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the Rician distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rice_distribution}
##
## @seealso{ricecdf, riceinv, ricepdf, ricefit, ricelike, ricestat}
## @end deftypefn

function r = ricernd (s, sigma, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("ricernd: function called with too few input arguments.");
  endif

  ## Check for common size of S and SIGMA
  if (! isscalar (s) || ! isscalar (sigma))
    [retval, s, sigma] = common_size (s, sigma);
    if (retval > 0)
      error ("ricernd: S and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for S and SIGMA being reals
  if (iscomplex (s) || iscomplex (sigma))
    error ("ricernd: S and SIGMA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (s);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif ((isrow (varargin{1}) || isempty (varargin{1})) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("ricernd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("ricernd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (s) && ! isequal (size (s), sz))
    error ("ricernd: S and SIGMA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (s, "single") || isa (sigma, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Return NaNs for out of range values of S and SIGMA
  s(s < 0) = NaN;
  sigma(sigma <= 0) = NaN;

  ## Force S and SIGMA into the same size as SZ (if necessary)
  if (isscalar (s))
    s = repmat (s, sz);
  endif
  if (isscalar (sigma))
    sigma = repmat (sigma, sz);
  endif

  ## Generate random sample from the Rician distribution
  r = sigma .* sqrt (ncx2rnd (2, (s ./ sigma) .^ 2, sz));
  ## Cast to appropriate class
  r = cast (r, cls);

endfunction

## Test output
%!assert (size (ricernd (2, 1/2)), [1, 1])
%!assert (size (ricernd (2 * ones (2, 1), 1/2)), [2, 1])
%!assert (size (ricernd (2 * ones (2, 2), 1/2)), [2, 2])
%!assert (size (ricernd (2, 1/2 * ones (2, 1))), [2, 1])
%!assert (size (ricernd (1, 1/2 * ones (2, 2))), [2, 2])
%!assert (size (ricernd (ones (2, 1), 1)), [2, 1])
%!assert (size (ricernd (ones (2, 2), 1)), [2, 2])
%!assert (size (ricernd (2, 1/2, 3)), [3, 3])
%!assert (size (ricernd (1, 1, [4, 1])), [4, 1])
%!assert (size (ricernd (1, 1, 4, 1)), [4, 1])
%!assert (size (ricernd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (ricernd (1, 1, 0, 1)), [0, 1])
%!assert (size (ricernd (1, 1, 1, 0)), [1, 0])
%!assert (size (ricernd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (ricernd (1, 1)), "double")
%!assert (class (ricernd (1, single (0))), "single")
%!assert (class (ricernd (1, single ([0, 0]))), "single")
%!assert (class (ricernd (1, single (1), 2)), "single")
%!assert (class (ricernd (1, single ([1, 1]), 1, 2)), "single")
%!assert (class (ricernd (single (1), 1, 2)), "single")
%!assert (class (ricernd (single ([1, 1]), 1, 1, 2)), "single")

## Test input validation
%!error<ricernd: function called with too few input arguments.> ricernd ()
%!error<ricernd: function called with too few input arguments.> ricernd (1)
%!error<ricernd: S and SIGMA must be of common size or scalars.> ...
%! ricernd (ones (3), ones (2))
%!error<ricernd: S and SIGMA must be of common size or scalars.> ...
%! ricernd (ones (2), ones (3))
%!error<ricernd: S and SIGMA must not be complex.> ricernd (i, 2)
%!error<ricernd: S and SIGMA must not be complex.> ricernd (1, i)
%!error<ricernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ricernd (1, 1/2, -1)
%!error<ricernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ricernd (1, 1/2, 1.2)
%!error<ricernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ricernd (1, 1/2, ones (2))
%!error<ricernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ricernd (1, 1/2, [2 -1 2])
%!error<ricernd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ricernd (1, 1/2, [2 0 2.5])
%!error<ricernd: dimensions must be non-negative integers.> ...
%! ricernd (1, 1/2, 2, -1, 5)
%!error<ricernd: dimensions must be non-negative integers.> ...
%! ricernd (1, 1/2, 2, 1.5, 5)
%!error<ricernd: S and SIGMA must be scalars or of size SZ.> ...
%! ricernd (2, 1/2 * ones (2), 3)
%!error<ricernd: S and SIGMA must be scalars or of size SZ.> ...
%! ricernd (2, 1/2 * ones (2), [3, 2])
%!error<ricernd: S and SIGMA must be scalars or of size SZ.> ...
%! ricernd (2, 1/2 * ones (2), 3, 2)
%!error<ricernd: S and SIGMA must be scalars or of size SZ.> ...
%! ricernd (2 * ones (2), 1/2, 3)
%!error<ricernd: S and SIGMA must be scalars or of size SZ.> ...
%! ricernd (2 * ones (2), 1/2, [3, 2])
%!error<ricernd: S and SIGMA must be scalars or of size SZ.> ...
%! ricernd (2 * ones (2), 1/2, 3, 2)
