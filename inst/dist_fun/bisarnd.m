## Copyright (C) 2018 John Donoghue
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
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
## @deftypefn  {statistics} {@var{r} =} bisarnd (@var{a}, @var{b}, @var{mu})
## @deftypefnx {statistics} {@var{r} =} bisarnd (@var{a}, @var{b}, @var{mu}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} bisarnd (@var{a}, @var{b}, @var{mu}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} bisarnd (@var{a}, @var{b}, @var{mu}, [@var{sz}])
##
## Random arrays from the Birnbaum-Saunders distribution.
##
## @code{@var{r} = bisarnd (@var{a}, @var{b}, @var{mu})} returns an array of
## random numbers chosen from the Birnbaum-Saunders distribution with shape
## parameter @var{a}, scale parameter @var{b}, and location parameter @var{mu}.
## The size of @var{r} is the common size of @var{a}, @var{b}, and @var{mu}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## When called with a single size argument, it returns a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a row vector of dimensions @var{sz}.
##
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bisacdf, bisainv, bisapdf, bisafit, bisalike, bisastat}
## @end deftypefn

function r = bisarnd (a, b, mu, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("bisarnd: function called with too few input arguments.");
  endif

  ## Check for common size of A, B, and MU
  if (! isscalar (a) || ! isscalar (b) || ! isscalar (mu))
    [retval, a, b, mu] = common_size (a, b, mu);
    if (retval > 0)
      error (strcat (["bisarnd: A, B, and MU must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for A, B, and MU being reals
  if (iscomplex (a) || iscomplex (b) || iscomplex (mu))
    error ("bisarnd: A, B, and MU must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 3)
    sz = size (mu);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["bisarnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 4)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("bisarnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("bisarnd: A, B, and MU must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single") || isa (b, "single") || isa (a, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Birnbaum-Saunders distribution
  if (isscalar (a) && isscalar (b) && isscalar (mu))
    if ((mu > -Inf) && (mu < Inf) && (b > 0) && (b < Inf) ...
                    && (a > 0) && (a < Inf))
      r = rand (sz, cls);
      y = a * norminv (r);
      r = mu + b * (y + sqrt (4 + y .^ 2)) .^ 2 / 4;
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);
    k = (mu > -Inf) & (mu < Inf) & (b > 0) & (b < Inf) & (a > 0) & (a < Inf);
    r(k) = rand (sum (k(:)),1);
    y = a(k) .* norminv (r(k));
    r(k) = mu(k) + b(k) .* (y + sqrt (4 + y.^2)).^2 / 4;
  endif
endfunction

## Test results
%!assert (size (bisarnd (1, 1, 0)), [1 1])
%!assert (size (bisarnd (1, 1, zeros (2,1))), [2, 1])
%!assert (size (bisarnd (1, 1, zeros (2,2))), [2, 2])
%!assert (size (bisarnd (1, ones (2,1), 0)), [2, 1])
%!assert (size (bisarnd (1, ones (2,2), 0)), [2, 2])
%!assert (size (bisarnd (ones (2,1), 1, 0)), [2, 1])
%!assert (size (bisarnd (ones (2,2), 1, 0)), [2, 2])
%!assert (size (bisarnd (1, 1, 0, 3)), [3, 3])
%!assert (size (bisarnd (1, 1, 0, [4, 1])), [4, 1])
%!assert (size (bisarnd (1, 1, 0, 4, 1)), [4, 1])
%!assert (size (bisarnd (1, 1, 0, 4, 1, 5)), [4, 1, 5])
%!assert (size (bisarnd (1, 1, 0, 0, 1)), [0, 1])
%!assert (size (bisarnd (1, 1, 0, 1, 0)), [1, 0])
%!assert (size (bisarnd (1, 1, 0, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (bisarnd (1, 1, 0)), "double")
%!assert (class (bisarnd (1, 1, single (0))), "single")
%!assert (class (bisarnd (1, 1, single ([0, 0]))), "single")
%!assert (class (bisarnd (1, single (1), 0)), "single")
%!assert (class (bisarnd (1, single ([1, 1]), 0)), "single")
%!assert (class (bisarnd (single (1), 1, 0)), "single")
%!assert (class (bisarnd (single ([1, 1]), 1, 0)), "single")

## Test input validation
%!error<bisarnd: function called with too few input arguments.> bisarnd ()
%!error<bisarnd: function called with too few input arguments.> bisarnd (1)
%!error<bisarnd: function called with too few input arguments.> bisarnd (1, 2)
%!error<bisarnd: A, B, and MU must be of common size or scalars.> ...
%! bisarnd (ones (3), ones (2), ones (2), 2)
%!error<bisarnd: A, B, and MU must be of common size or scalars.> ...
%! bisarnd (ones (2), ones (3), ones (2), 2)
%!error<bisarnd: A, B, and MU must be of common size or scalars.> ...
%! bisarnd (ones (2), ones (2), ones (3), 2)
%!error<bisarnd: A, B, and MU must not be complex.> bisarnd (i, 2, 3)
%!error<bisarnd: A, B, and MU must not be complex.> bisarnd (1, i, 3)
%!error<bisarnd: A, B, and MU must not be complex.> bisarnd (1, 2, i)
%!error<bisarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! bisarnd (1, 2, 3, -1)
%!error<bisarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! bisarnd (1, 2, 3, 1.2)
%!error<bisarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! bisarnd (1, 2, 3, ones (2))
%!error<bisarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! bisarnd (1, 2, 3, [2 -1 2])
%!error<bisarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! bisarnd (1, 2, 3, [2 0 2.5])
%!error<bisarnd: dimensions must be non-negative integers.> ...
%! bisarnd (1, 2, 3, 2, -1, 5)
%!error<bisarnd: dimensions must be non-negative integers.> ...
%! bisarnd (1, 2, 3, 2, 1.5, 5)
%!error<bisarnd: A, B, and MU must be scalar or of size SZ.> ...
%! bisarnd (2, 1, ones (2), 3)
%!error<bisarnd: A, B, and MU must be scalar or of size SZ.> ...
%! bisarnd (2, 1, ones (2), [3, 2])
%!error<bisarnd: A, B, and MU must be scalar or of size SZ.> ...
%! bisarnd (2, 1, ones (2), 3, 2)
