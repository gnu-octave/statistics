## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
##
## This file is part of Octave.
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
## @deftypefn  {} {} burrrnd (@var{alpha}, @var{c}, @var{k})
## @deftypefnx {} {} burrrnd (@var{alpha}, @var{c}, @var{k}, @var{r})
## @deftypefnx {} {} burrrnd (@var{alpha}, @var{c}, @var{k}, @var{r}, @var{c}, @dots{})
## @deftypefnx {} {} burrrnd (@var{alpha}, @var{c}, @var{k}, [@var{sz}])
## Return a matrix of random samples from the generalized Pareto distribution
## with scale parameter @var{alpha} and shape parameters @var{c} and @var{k}.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## If no size arguments are given then the result matrix is the common size of
## @var{alpha}, @var{c} and @var{k}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: Random deviates from the generalized extreme value (GEV) distribution

function rnd = burrrnd (alpha, c, k, varargin)

  if (nargin < 3)
    print_usage ();
  endif

  if (! isscalar (alpha) || ! isscalar (c) || ! isscalar (k))
    [retval, alpha, c, k] = common_size (alpha, c, k);
    if (retval > 0)
      error ("burrrnd: ALPHA, C and K must be of common size or scalars");
    endif
  endif

  if (iscomplex (alpha) || iscomplex (c) || iscomplex (k))
    error ("burrrnd: ALPHA, C and K must not be complex");
  endif

  if (nargin == 3)
    sz = size (alpha);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error ("burrrnd: dimension vector must be row vector of non-negative integers");
    endif
  elseif (nargin > 4)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("burrrnd: dimensions must be non-negative integers");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (alpha) && ! isequal (size (c), sz) && ! isequal (size (k), sz))
    error ("burrrnd: ALPHA, C and K must be scalar or of size SZ");
  endif

  if (isa (alpha, "single") || isa (c, "single") || isa (k, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (alpha) && isscalar (c) && isscalar(k))
    if ((0 < alpha) && (alpha < Inf) && (0 < c) && (c < Inf) ...
        && (0 < k) && (k < Inf))
      rnd = rand (sz, cls);
      rnd(:) = ((1 - rnd(:) / alpha).^(-1 / k) - 1).^(1 / c);
    else
      rnd = NaN (sz, cls);
    endif
  else
    rnd = NaN (sz, cls);

    j = (0 < alpha) && (alpha < Inf) && (0 < c) && (c < Inf) ...
        && (0 < k) && (k < Inf);
    rnd(k) = rand(sum(j(:)),1);
    rnd(k) = ((1 - rnd(j) / alpha(j)).^(-1 ./ k(j)) - 1).^(1 ./ c(j));
  endif

endfunction


%!assert (size (burrrnd (1, 1, 1)), [1 1])
%!assert (size (burrrnd (ones (2,1), 1, 1)), [2, 1])
%!assert (size (burrrnd (ones (2,2), 1, 1)), [2, 2])
%!assert (size (burrrnd (1, ones (2,1), 1)), [2, 1])
%!assert (size (burrrnd (1, ones (2,2), 1)), [2, 2])
%!assert (size (burrrnd (1, 1, ones (2,1))), [2, 1])
%!assert (size (burrrnd (1, 1, ones (2,2))), [2, 2])
%!assert (size (burrrnd (1, 1, 1, 3)), [3, 3])
%!assert (size (burrrnd (1, 1, 1, [4 1])), [4, 1])
%!assert (size (burrrnd (1, 1, 1, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (burrrnd (1,1,1)), "double")
%!assert (class (burrrnd (single (1),1,1)), "single")
%!assert (class (burrrnd (single ([1 1]),1,1)), "single")
%!assert (class (burrrnd (1,single (1),1)), "single")
%!assert (class (burrrnd (1,single ([1 1]),1)), "single")
%!assert (class (burrrnd (1,1,single (1))), "single")
%!assert (class (burrrnd (1,1,single ([1 1]))), "single")

## Test input validation
%!error burrrnd ()
%!error burrrnd (1)
%!error burrrnd (1,2)
%!error burrrnd (ones (3), ones (2), ones (2), 2)
%!error burrrnd (ones (2), ones (3), ones (2), 2)
%!error burrrnd (ones (2), ones (2), ones (3), 2)
%!error burrrnd (i, 2, 2)
%!error burrrnd (2, i, 2)
%!error burrrnd (2, 2, i)
%!error burrrnd (4,2,2, -1)
%!error burrrnd (4,2,2, ones (2))
%!error burrrnd (4,2,2, [2 -1 2])
%!error burrrnd (4*ones (2),2,2, 3)
%!error burrrnd (4*ones (2),2,2, [3, 2])
%!error burrrnd (4*ones (2),2,2, 3, 2)

