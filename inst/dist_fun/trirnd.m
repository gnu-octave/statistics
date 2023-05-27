## Copyright (C) 1997-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} trirnd (@var{a}, @var{b}, @var{c})
## @deftypefnx {statistics} {@var{r} =} trirnd (@var{a}, @var{b}, @var{c}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} trirnd (@var{a}, @var{b}, @var{c}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} trirnd (@var{a}, @var{b}, @var{c}, [@var{sz}])
##
## Random arrays from the triangular distribution.
##
## @code{@var{r} = trirnd (@var{sigma})} returns an array of random numbers
## chosen from the triangular distribution with parameters @var{a}, @var{b}, and
## @var{c} on the interval [@var{a}, @var{b}].  The size of @var{r} is the
## common size of @var{a}, @var{b}, and @var{c}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## When called with a single size argument, @code{trirnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the triangular distribution can be found at
## @url{https://en.wikipedia.org/wiki/Triangular_distribution}
##
## @seealso{tricdf, triinv, tripdf}
## @end deftypefn

function rnd = trirnd (a, b, c, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("trirnd: function called with too few input arguments.");
  endif

  ## Check for common size of A, B, and C
  if (! isscalar (a) || ! isscalar (b) || ! isscalar (c))
    [retval, a, b, c] = common_size (a, b, c);
    scalarABC = false;
    if (retval > 0)
      error ("trirnd: A, B, and C must be of common size or scalars.");
    endif
  else
    scalarABC = true;
  endif

  ## Check for A, B, and C being reals
  if (iscomplex (a) || iscomplex (b) || iscomplex (c))
    error ("trirnd: A, B, and C must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 3)
    sz = size (a);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["trirnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 4)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("trirnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (a) && ! isequal (size (a), sz))
    error ("trirnd: A, B, and C must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (a, "single") || isa (b, "single") || isa (c, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from triangular distribution
  if (scalarABC)
    if ((-Inf < a) && (a < b) && (a <= c) && (c <= b) && (b < Inf))
      w = b-a;
      left_width = c-a;
      right_width = b-c;
      h = 2 / w;
      left_area = h * left_width / 2;
      rnd = rand (sz, cls);
      idx = rnd < left_area;
      rnd(idx) = a + (rnd(idx) * w * left_width).^0.5;
      rnd(~idx) = b - ((1-rnd(~idx)) * w * right_width).^0.5;
    else
      rnd = NaN (sz, cls);
    endif
  else
    w = b-a;
    left_width = c-a;
    right_width = b-c;
    h = 2 ./ w;
    left_area = h .* left_width / 2;
    rnd = rand (sz, cls);
    k = rnd < left_area;
    rnd(k) = a(k) + (rnd(k) .* w(k) .* left_width(k)).^0.5;
    rnd(~k) = b(~k) - ((1-rnd(~k)) .* w(~k) .* right_width(~k)).^0.5;

    k = ! (-Inf < a) | ! (a < b) | ! (a <= c) | ! (c <= b) | ! (b < Inf);
    rnd(k) = NaN;
  endif

endfunction

## Test results
%!assert (size (trirnd (1,2,1.5)), [1, 1])
%!assert (size (trirnd (1*ones (2,1), 2,1.5)), [2, 1])
%!assert (size (trirnd (1*ones (2,2), 2,1.5)), [2, 2])
%!assert (size (trirnd (1, 2*ones (2,1), 1.5)), [2, 1])
%!assert (size (trirnd (1, 2*ones (2,2), 1.5)), [2, 2])
%!assert (size (trirnd (1, 2, 1.5*ones (2,1))), [2, 1])
%!assert (size (trirnd (1, 2, 1.5*ones (2,2))), [2, 2])
%!assert (size (trirnd (1, 2, 1.5, 3)), [3, 3])
%!assert (size (trirnd (1, 2, 1.5, [4 1])), [4, 1])
%!assert (size (trirnd (1, 2, 1.5, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (trirnd (1,2,1.5)), "double")
%!assert (class (trirnd (single (1),2,1.5)), "single")
%!assert (class (trirnd (single ([1 1]),2,1.5)), "single")
%!assert (class (trirnd (1,single (2),1.5)), "single")
%!assert (class (trirnd (1,single ([2 2]),1.5)), "single")
%!assert (class (trirnd (1,2,single (1.5))), "single")
%!assert (class (trirnd (1,2,single ([1.5 1.5]))), "single")

## Test input validation
%!error<trirnd: function called with too few input arguments.> trirnd ()
%!error<trirnd: function called with too few input arguments.> trirnd (1)
%!error<trirnd: function called with too few input arguments.> trirnd (1, 2)
%!error<trirnd: A, B, and C must be of common size or scalars.> ...
%! trirnd (ones (3), 5 * ones (2), ones (2))
%!error<trirnd: A, B, and C must be of common size or scalars.> ...
%! trirnd (ones (2), 5 * ones (3), ones (2))
%!error<trirnd: A, B, and C must be of common size or scalars.> ...
%! trirnd (ones (2), 5 * ones (2), ones (3))
%!error<trirnd: A, B, and C must not be complex.> trirnd (i, 5, 3)
%!error<trirnd: A, B, and C must not be complex.> trirnd (1, 5+i, 3)
%!error<trirnd: A, B, and C must not be complex.> trirnd (1, 5, i)
%!error<trirnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! trirnd (1, 5, 3, -1)
%!error<trirnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! trirnd (1, 5, 3, 1.2)
%!error<trirnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! trirnd (1, 5, 3, ones (2))
%!error<trirnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! trirnd (1, 5, 3, [2 -1 2])
%!error<trirnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! trirnd (1, 5, 3, [2 0 2.5])
%!error<trirnd: dimensions must be non-negative integers.> ...
%! trirnd (1, 5, 3, 2, -1, 5)
%!error<trirnd: dimensions must be non-negative integers.> ...
%! trirnd (1, 5, 3, 2, 1.5, 5)
%!error<trirnd: A, B, and C must be scalar or of size SZ.> ...
%! trirnd (2, 5 * ones (2), 2, 3)
%!error<trirnd: A, B, and C must be scalar or of size SZ.> ...
%! trirnd (2, 5 * ones (2), 2, [3, 2])
%!error<trirnd: A, B, and C must be scalar or of size SZ.> ...
%! trirnd (2, 5 * ones (2), 2, 3, 2)
