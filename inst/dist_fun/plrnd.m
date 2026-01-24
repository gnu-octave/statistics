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
## @deftypefn  {statistics} {@var{r} =} plrnd (@var{x}, @var{Fx})
## @deftypefnx {statistics} {@var{r} =} plrnd (@var{x}, @var{Fx}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} plrnd (@var{x}, @var{Fx}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} plrnd (@var{x}, @var{Fx}, [@var{sz}])
##
## Random arrays from the piecewise linear distribution.
##
## @code{@var{r} = plrnd (@var{x}, @var{Fx})} returns a random number chosen
## from the piecewise linear distribution with a vector of @var{x} values at
## which the CDF changes slope and a vector of CDF values @var{Fx} that
## correspond to each value in @var{x}.  Both @var{x} and @var{Fx} must be
## vectors of the same size and at least 2-elements long.
##
## When called with a single size argument, @code{plrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the piecewise linear distribution can be found at
## @url{https://en.wikipedia.org/wiki/Piecewise_linear_function}
##
## @seealso{plcdf, plinv, plpdf, plstat}
## @end deftypefn

function r = plrnd (x, Fx, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("plrnd: function called with too few input arguments.");
  endif

  ## Check for common size of X and FX
  if (! isvector (x) || ! isvector (Fx) || ! isequal (size (x), size (Fx)))
    error ("plrnd: X and FX must be vectors of equal size.");
  endif

  ## Check for X and FX being at least 2-elements long
  if (length (x) < 2 || length (Fx) < 2)
    error ("plrnd: X and FX must be at least two-elements long.");
  endif

  ## Check for Fx being bounded in [0, 1]
  if (any (Fx < 0) || any (Fx > 1))
    error ("plrnd: FX must be bounded in the range [0, 1].");
  endif

  ## Check for X and FX being reals
  if (iscomplex (x) || iscomplex (Fx))
    error ("plrnd: X and FX must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = 1;
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif ((isrow (varargin{1}) || isempty (varargin{1})) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("plrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("plrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Force X and FX into row vectors
  x = x(:)';
  Fx = Fx(:)';

  ## Check for class type
  if (isa (x, "single") || isa (Fx, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from the piecewise linear distribution
  u = rand (sz);
  r = zeros (sz);
  [~, bin] = histc (u(:)', Fx);
  r0 = x(bin);
  dx = diff (x);
  dF = diff (Fx);
  dr = (u(:)' - Fx(bin)) .* dx(bin) ./ dF(bin);
  r(:) = r0 + dr;

  ## Cast to appropriate class
  r = cast (r, cls);

endfunction

## Test output
%!shared x, Fx
%! x = [0, 1, 3, 4, 7, 10];
%! Fx = [0, 0.2, 0.5, 0.6, 0.7, 1];
%!assert (size (plrnd (x, Fx)), [1, 1])
%!assert (size (plrnd (x, Fx, 3)), [3, 3])
%!assert (size (plrnd (x, Fx, [4, 1])), [4, 1])
%!assert (size (plrnd (x, Fx, 4, 1)), [4, 1])
%!assert (size (plrnd (x, Fx, 4, 1, 5)), [4, 1, 5])
%!assert (size (plrnd (x, Fx, 0, 1)), [0, 1])
%!assert (size (plrnd (x, Fx, 1, 0)), [1, 0])
%!assert (size (plrnd (x, Fx, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (plrnd (x, Fx)), "double")
%!assert (class (plrnd (x, single (Fx))), "single")
%!assert (class (plrnd (single (x), Fx)), "single")

## Test input validation
%!error<plrnd: function called with too few input arguments.> plrnd ()
%!error<plrnd: function called with too few input arguments.> plrnd (1)
%!error<plrnd: X and FX must be vectors of equal size.> ...
%! plrnd ([0, 1, 2], [0, 1])
%!error<plrnd: X and FX must be at least two-elements long.> ...
%! plrnd ([0], [1])
%!error<plrnd: FX must be bounded in the range> ...
%! plrnd ([0, 1, 2], [0, 1, 1.5])
%!error<plrnd: FX must be bounded in the range> ...
%! plrnd ([0, 1, 2], [0, i, 1])
%!error<plrnd: X and FX must not be complex.> ...
%! plrnd ([0, i, 2], [0, 0.5, 1])
%!error<plrnd: X and FX must not be complex.> ...
%! plrnd ([0, i, 2], [0, 0.5i, 1])
%!error<plrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! plrnd (x, Fx, -1)
%!error<plrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! plrnd (x, Fx, 1.2)
%!error<plrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! plrnd (x, Fx, ones (2))
%!error<plrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! plrnd (x, Fx, [2 -1 2])
%!error<plrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! plrnd (x, Fx, [2 0 2.5])
%!error<plrnd: dimensions must be non-negative integers.> ...
%! plrnd (x, Fx, 2, -1, 5)
%!error<plrnd: dimensions must be non-negative integers.> ...
%! plrnd (x, Fx, 2, 1.5, 5)
