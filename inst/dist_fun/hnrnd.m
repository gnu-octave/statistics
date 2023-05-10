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
## @deftypefn  {statistics} {@var{r} =} hnrnd (@var{mu}, @var{sigma})
## @deftypefnx {statistics} {@var{r} =} hnrnd (@var{mu}, @var{sigma}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} hnrnd (@var{mu}, @var{sigma}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} hnrnd (@var{mu}, @var{sigma}, [@var{sz}])
##
## Random arrays from the half-normal distribution.
##
## @code{@var{r} = hnrnd (@var{mu}, @var{sigma})} returns an array of random
## numbers chosen from the half-normal distribution with location parameter
## @var{mu} and scale parameter @var{sigma}.  The size of @var{r} is the common
## size of @var{mu} and @var{sigma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## When called with a single size argument, @code{hnrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the half-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Half-normal_distribution}
##
## @seealso{hncdf, hninv, hnpdf, hnfit, hnlike}
## @end deftypefn

function r = hnrnd (mu, sigma, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("hnrnd: function called with too few input arguments.");
  endif

  ## Check for common size of MU, and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("hnrnd: MU and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (mu) || iscomplex (sigma))
    error ("hnrnd: MU and SIGMA must not be complex.");
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
      error (strcat (["hnrnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("hnrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("hnrnd: MU and SIGMA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single") || isa (sigma, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from half-normal distribution
  r = abs (randn (sz, cls)) .* sigma + mu;

  ## Force output to NaN for invalid parameter SIGMA <= 0
  k = (sigma <= 0);
  r(k) = NaN;

endfunction

## Test output
%!assert (size (hnrnd (1, 1, 1)), [1, 1])
%!assert (size (hnrnd (1, 1, 2)), [2, 2])
%!assert (size (hnrnd (1, 1, [2, 1])), [2, 1])
%!assert (size (hnrnd (1, zeros (2, 2))), [2, 2])
%!assert (size (hnrnd (1, ones (2, 1))), [2, 1])
%!assert (size (hnrnd (1, ones (2, 2))), [2, 2])
%!assert (size (hnrnd (ones (2, 1), 1)), [2, 1])
%!assert (size (hnrnd (ones (2, 2), 1)), [2, 2])
%!assert (size (hnrnd (1, 1, 3)), [3, 3])
%!assert (size (hnrnd (1, 1, [4 1])), [4, 1])
%!assert (size (hnrnd (1, 1, 4, 1)), [4, 1])
%!test
%! r =  hnrnd (1, [1, 0, -1]);
%! assert (r([2:3]), [NaN, NaN])

## Test class of input preserved
%!assert (class (hnrnd (1, 0)), "double")
%!assert (class (hnrnd (1, single (0))), "single")
%!assert (class (hnrnd (1, single ([0 0]))), "single")
%!assert (class (hnrnd (1, single (1))), "single")
%!assert (class (hnrnd (1, single ([1 1]))), "single")
%!assert (class (hnrnd (single (1), 1)), "single")
%!assert (class (hnrnd (single ([1 1]), 1)), "single")

## Test input validation
%!error<hnrnd: function called with too few input arguments.> hnrnd ()
%!error<hnrnd: function called with too few input arguments.> hnrnd (1)
%!error<hnrnd: MU and SIGMA must be of common size or scalars.> ...
%! hnrnd (ones (3), ones (2))
%!error<hnrnd: MU and SIGMA must be of common size or scalars.> ...
%! hnrnd (ones (2), ones (3))
%!error<hnrnd: MU and SIGMA must not be complex.> hnrnd (i, 2, 3)
%!error<hnrnd: MU and SIGMA must not be complex.> hnrnd (1, i, 3)
%!error<hnrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hnrnd (1, 2, -1)
%!error<hnrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hnrnd (1, 2, 1.2)
%!error<hnrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hnrnd (1, 2, ones (2))
%!error<hnrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hnrnd (1, 2, [2 -1 2])
%!error<hnrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! hnrnd (1, 2, [2 0 2.5])
%!error<hnrnd: dimensions must be non-negative integers.> ...
%! hnrnd (1, 2, 2, -1, 5)
%!error<hnrnd: dimensions must be non-negative integers.> ...
%! hnrnd (1, 2, 2, 1.5, 5)
%!error<hnrnd: MU and SIGMA must be scalars or of size SZ.> ...
%! hnrnd (2, ones (2), 3)
%!error<hnrnd: MU and SIGMA must be scalars or of size SZ.> ...
%! hnrnd (2, ones (2), [3, 2])
%!error<hnrnd: MU and SIGMA must be scalars or of size SZ.> ...
%! hnrnd (2, ones (2), 3, 2)
