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
## @deftypefn  {statistics} {@var{r} =} cauchyrnd (@var{x0}, @var{gamma})
## @deftypefnx {statistics} {@var{r} =} cauchyrnd (@var{x0}, @var{gamma}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} cauchyrnd (@var{x0}, @var{gamma}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} cauchyrnd (@var{x0}, @var{gamma}, [@var{sz}])
##
## Random arrays from the Cauchy distribution.
##
## @code{@var{r} = cauchyrnd (@var{x0}, @var{gamma})} returns an array of
## random numbers chosen from the Cauchy distribution with location parameter
## @var{x0} and scale parameter @var{gamma}.  The size of @var{r} is the common
## size of @var{x0} and @var{gamma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## When called with a single size argument, @code{cauchyrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the Cauchy distribution can be found at
## @url{https://en.wikipedia.org/wiki/Cauchy_distribution}
##
## @seealso{cauchycdf, cauchyinv, cauchypdf}
## @end deftypefn

function r = cauchyrnd (x0, gamma, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("cauchyrnd: function called with too few input arguments.");
  endif

  ## Check for common size of X0 and GAMMA
  if (! isscalar (x0) || ! isscalar (gamma))
    [retval, x0, gamma] = common_size (x0, gamma);
    if (retval > 0)
      error (strcat (["cauchyrnd: X0 and GAMMA must be of common"], ...
                     [" size or scalars."]));
    endif
  endif

  ## Check for X0 and GAMMA being reals
  if (iscomplex (x0) || iscomplex (gamma))
    error ("cauchyrnd: X0 and GAMMA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (x0);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["cauchyrnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("cauchyrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (x0) && ! isequal (size (x0), sz))
    error ("cauchyrnd: X0 and GAMMA must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (x0, "single") || isa (gamma, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Cauchy distribution
  if (isscalar (x0) && isscalar (gamma))
    if (! isinf (x0) && (gamma > 0) && (gamma < Inf))
      r = x0 - cot (pi * rand (sz, cls)) * gamma;
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    k = ! isinf (x0) & (gamma > 0) & (gamma < Inf);
    r(k) = x0(k)(:) - cot (pi * rand (sum (k(:)), 1, cls)) .* gamma(k)(:);
  endif

endfunction

## Test output
%!assert (size (cauchyrnd (1, 1)), [1 1])
%!assert (size (cauchyrnd (1, ones (2,1))), [2, 1])
%!assert (size (cauchyrnd (1, ones (2,2))), [2, 2])
%!assert (size (cauchyrnd (ones (2,1), 1)), [2, 1])
%!assert (size (cauchyrnd (ones (2,2), 1)), [2, 2])
%!assert (size (cauchyrnd (1, 1, 3)), [3, 3])
%!assert (size (cauchyrnd (1, 1, [4, 1])), [4, 1])
%!assert (size (cauchyrnd (1, 1, 4, 1)), [4, 1])
%!assert (size (cauchyrnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (cauchyrnd (1, 1, 0, 1)), [0, 1])
%!assert (size (cauchyrnd (1, 1, 1, 0)), [1, 0])
%!assert (size (cauchyrnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (cauchyrnd (1, 1)), "double")
%!assert (class (cauchyrnd (1, single (1))), "single")
%!assert (class (cauchyrnd (1, single ([1, 1]))), "single")
%!assert (class (cauchyrnd (single (1), 1)), "single")
%!assert (class (cauchyrnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<cauchyrnd: function called with too few input arguments.> cauchyrnd ()
%!error<cauchyrnd: function called with too few input arguments.> cauchyrnd (1)
%!error<cauchyrnd: X0 and GAMMA must be of common size or scalars.> ...
%! cauchyrnd (ones (3), ones (2))
%!error<cauchyrnd: X0 and GAMMA must be of common size or scalars.> ...
%! cauchyrnd (ones (2), ones (3))
%!error<cauchyrnd: X0 and GAMMA must not be complex.> cauchyrnd (i, 2, 3)
%!error<cauchyrnd: X0 and GAMMA must not be complex.> cauchyrnd (1, i, 3)
%!error<cauchyrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! cauchyrnd (1, 2, -1)
%!error<cauchyrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! cauchyrnd (1, 2, 1.2)
%!error<cauchyrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! cauchyrnd (1, 2, ones (2))
%!error<cauchyrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! cauchyrnd (1, 2, [2 -1 2])
%!error<cauchyrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! cauchyrnd (1, 2, [2 0 2.5])
%!error<cauchyrnd: dimensions must be non-negative integers.> ...
%! cauchyrnd (1, 2, 2, -1, 5)
%!error<cauchyrnd: dimensions must be non-negative integers.> ...
%! cauchyrnd (1, 2, 2, 1.5, 5)
%!error<cauchyrnd: X0 and GAMMA must be scalar or of size SZ.> ...
%! cauchyrnd (2, ones (2), 3)
%!error<cauchyrnd: X0 and GAMMA must be scalar or of size SZ.> ...
%! cauchyrnd (2, ones (2), [3, 2])
%!error<cauchyrnd: X0 and GAMMA must be scalar or of size SZ.> ...
%! cauchyrnd (2, ones (2), 3, 2)
