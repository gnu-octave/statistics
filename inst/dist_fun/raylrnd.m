## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
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
## @deftypefn  {statistics} {@var{r} =} raylrnd (@var{sigma})
## @deftypefnx {statistics} {@var{r} =} raylrnd (@var{sigma}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} raylrnd (@var{sigma}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} raylrnd (@var{sigma}, [@var{sz}])
##
## Random arrays from the Rayleigh distribution.
##
## @code{@var{r} = raylrnd (@var{sigma})} returns an array of random numbers
## chosen from the Rayleigh distribution with scale parameter @var{sigma}.  The
## size of @var{r} is the size of @var{sigma}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.  @var{sigma} must be a
## finite real number greater than 0, otherwise @qcode{NaN} is returned.
##
## When called with a single size argument, @code{raylrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the Rayleigh distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rayleigh_distribution}
##
## @seealso{raylcdf, raylinv, raylpdf, raylfit, rayllike, raylstat}
## @end deftypefn

function r = raylrnd (sigma, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("raylrnd: function called with too few input arguments.");
  endif

  ## Check for SIGMA being real
  if (iscomplex (sigma))
    error ("raylrnd: SIGMA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 1)
    sz = size (sigma);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["raylrnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 2)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("raylrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (sigma) && ! isequal (size (sigma), sz))
    error ("raylrnd: SIGMA must be scalar or of size SZ.");
  endif

  ## Generate random sample from Rayleigh distribution
  r = sqrt (-2 .* log (1 - rand (sz)) .* sigma .^ 2);

  ## Check for valid parameter
  k = find (! (sigma > 0));
  if (any (k))
    r(k) = NaN;
  endif

  ## Cast into appropriate class
  if (isa (sigma, "single"))
    r = cast (r, "single");
  endif

endfunction

## Test output
%!assert (size (raylrnd (2)), [1, 1])
%!assert (size (raylrnd (ones (2,1))), [2, 1])
%!assert (size (raylrnd (ones (2,2))), [2, 2])
%!assert (size (raylrnd (1, 3)), [3, 3])
%!assert (size (raylrnd (1, [4 1])), [4, 1])
%!assert (size (raylrnd (1, 4, 1)), [4, 1])
%!assert (size (raylrnd (1, 4, 1)), [4, 1])
%!assert (size (raylrnd (1, 4, 1, 5)), [4, 1, 5])
%!assert (size (raylrnd (1, 0, 1)), [0, 1])
%!assert (size (raylrnd (1, 1, 0)), [1, 0])
%!assert (size (raylrnd (1, 1, 2, 0, 5)), [1, 2, 0, 5])
%!assert (raylrnd (0, 1, 1), NaN)
%!assert (raylrnd ([0, 0, 0], [1, 3]), [NaN, NaN, NaN])

## Test class of input preserved
%!assert (class (raylrnd (2)), "double")
%!assert (class (raylrnd (single (2))), "single")
%!assert (class (raylrnd (single ([2 2]))), "single")

## Test input validation
%!error<raylrnd: function called with too few input arguments.> raylrnd ()
%!error<raylrnd: SIGMA must not be complex.> raylrnd (i)
%!error<raylrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! raylrnd (1, -1)
%!error<raylrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! raylrnd (1, 1.2)
%!error<raylrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! raylrnd (1, ones (2))
%!error<raylrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! raylrnd (1, [2 -1 2])
%!error<raylrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! raylrnd (1, [2 0 2.5])
%!error<raylrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! raylrnd (ones (2), ones (2))
%!error<raylrnd: dimensions must be non-negative integers.> ...
%! raylrnd (1, 2, -1, 5)
%!error<raylrnd: dimensions must be non-negative integers.> ...
%! raylrnd (1, 2, 1.5, 5)
%!error<raylrnd: SIGMA must be scalar or of size SZ.> raylrnd (ones (2,2), 3)
%!error<raylrnd: SIGMA must be scalar or of size SZ.> raylrnd (ones (2,2), [3, 2])
%!error<raylrnd: SIGMA must be scalar or of size SZ.> raylrnd (ones (2,2), 2, 3)
