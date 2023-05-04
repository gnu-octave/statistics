## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{r} =} betarnd (@var{a}, @var{b})
## @deftypefnx {statistics} {@var{r} =} betarnd (@var{a}, @var{b}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} betarnd (@var{a}, @var{b}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} betarnd (@var{a}, @var{b}, [@var{sz}])
##
## Random arrays from the Beta distribution.
##
## @code{@var{r} = betarnd (@var{a}, @var{b})} returns an array of random
## numbers chosen from the Beta distribution with shape parameters @var{a} and
## @var{b}.  The size of @var{r} is the common size of @var{a} and @var{b}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## Further information about the Beta distribution can be found at
## @url{https://en.wikipedia.org/wiki/Beta_distribution}
##
## @seealso{betacdf, betainv, betapdf, betafit, betalike, betastat}
## @end deftypefn

function r = betarnd (a, b, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("betarnd: function called with too few input arguments.");
  endif

  ## Check for common size of A and B
  if (! isscalar (a) || ! isscalar (b))
    [retval, a, b] = common_size (a, b);
    if (retval > 0)
      error ("betarnd: A and B must be of common size or scalars.");
    endif
  endif

  ## Check for A and B being reals
  if (iscomplex (a) || iscomplex (b))
    error ("betarnd: A and B must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (a);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["betarnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("betarnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (a) && ! isequal (size (a), sz))
    error ("betarnd: A and B must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (a, "single") || isa (b, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Beta distribution
  if (isscalar (a) && isscalar (b))
    if ((a > 0) && (a < Inf) && (b > 0) && (b < Inf))
      tmpr = randg (a, sz, cls);
      r = tmpr ./ (tmpr + randg (b, sz, cls));
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    k = (a > 0) & (a < Inf) & (b > 0) & (b < Inf);
    tmpr = randg (a(k), cls);
    r(k) = tmpr ./ (tmpr + randg (b(k), cls));
  endif

endfunction

## Test output
%!assert (size (betarnd (2, 1/2)), [1 1])
%!assert (size (betarnd (2 * ones (2, 1), 1/2)), [2, 1])
%!assert (size (betarnd (2 * ones (2, 2), 1/2)), [2, 2])
%!assert (size (betarnd (2, 1/2 * ones (2, 1))), [2, 1])
%!assert (size (betarnd (1, 1/2 * ones (2, 2))), [2, 2])
%!assert (size (betarnd (ones (2, 1), 1)), [2, 1])
%!assert (size (betarnd (ones (2, 2), 1)), [2, 2])
%!assert (size (betarnd (2, 1/2, 3)), [3, 3])
%!assert (size (betarnd (1, 1, [4, 1])), [4, 1])
%!assert (size (betarnd (1, 1, 4, 1)), [4, 1])
%!assert (size (betarnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (betarnd (1, 1, 0, 1)), [0, 1])
%!assert (size (betarnd (1, 1, 1, 0)), [1, 0])
%!assert (size (betarnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (betarnd (1, 1)), "double")
%!assert (class (betarnd (1, single (0))), "single")
%!assert (class (betarnd (1, single ([0, 0]))), "single")
%!assert (class (betarnd (1, single (1), 2)), "single")
%!assert (class (betarnd (1, single ([1, 1]), 1, 2)), "single")
%!assert (class (betarnd (single (1), 1, 2)), "single")
%!assert (class (betarnd (single ([1, 1]), 1, 1, 2)), "single")

## Test input validation
%!error<betarnd: function called with too few input arguments.> betarnd ()
%!error<betarnd: function called with too few input arguments.> betarnd (1)
%!error<betarnd: A and B must be of common size or scalars.> ...
%! betarnd (ones (3), ones (2))
%!error<betarnd: A and B must be of common size or scalars.> ...
%! betarnd (ones (2), ones (3))
%!error<betarnd: A and B must not be complex.> betarnd (i, 2)
%!error<betarnd: A and B must not be complex.> betarnd (1, i)
%!error<betarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! betarnd (1, 1/2, -1)
%!error<betarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! betarnd (1, 1/2, 1.2)
%!error<betarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! betarnd (1, 1/2, ones (2))
%!error<betarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! betarnd (1, 1/2, [2 -1 2])
%!error<betarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! betarnd (1, 1/2, [2 0 2.5])
%!error<betarnd: dimensions must be non-negative integers.> ...
%! betarnd (1, 1/2, 2, -1, 5)
%!error<betarnd: dimensions must be non-negative integers.> ...
%! betarnd (1, 1/2, 2, 1.5, 5)
%!error<betarnd: A and B must be scalar or of size SZ.> ...
%! betarnd (2, 1/2 * ones (2), 3)
%!error<betarnd: A and B must be scalar or of size SZ.> ...
%! betarnd (2, 1/2 * ones (2), [3, 2])
%!error<betarnd: A and B must be scalar or of size SZ.> ...
%! betarnd (2, 1/2 * ones (2), 3, 2)
%!error<betarnd: A and B must be scalar or of size SZ.> ...
%! betarnd (2 * ones (2), 1/2, 3)
%!error<betarnd: A and B must be scalar or of size SZ.> ...
%! betarnd (2 * ones (2), 1/2, [3, 2])
%!error<betarnd: A and B must be scalar or of size SZ.> ...
%! betarnd (2 * ones (2), 1/2, 3, 2)
