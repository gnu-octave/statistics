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
## @deftypefn  {} {} bbsrnd (@var{location}, @var{scale}, @var{shape})
## @deftypefnx {} {} bbsrnd (@var{location}, @var{scale}, @var{shape}, @var{r})
## @deftypefnx {} {} bbsrnd (@var{location}, @var{scale}, @var{shape}, @var{r}, @var{c}, @dots{})
## @deftypefnx {} {} bbsrnd (@var{location}, @var{scale}, @var{shape}, [@var{sz}])
## Return a matrix of random samples from the Birnbaum-Saunders
##  distribution with parameters @var{location}, @var{scale} and @var{shape}.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## If no size arguments are given then the result matrix is the common size of
## @var{location}, @var{scale} and @var{shape}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: Random deviates from the Birnbaum-Saunders distribution

function rnd = bbsrnd (location, scale, shape, varargin)

  if (nargin < 3)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale) || ! isscalar (shape))
    [retval, location, scale, shape] = common_size (location, scale, shape);
    if (retval > 0)
      error ("bbsrnd: LOCATION, SCALE and SHAPE must be of common size or scalars");
    endif
  endif

  if (iscomplex (location) || iscomplex (scale) || iscomplex (shape))
    error ("bbsrnd: LOCATION, SCALE and SHAPE must not be complex");
  endif

  if (nargin == 3)
    sz = size (location);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error ("bbsrnd: dimension vector must be row vector of non-negative integers");
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("bbsrnd: dimensions must be non-negative integers");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (location) && ! isequal (size (location), sz))
    error ("bbsrnd: LOCATION, SCALE and SHAPE must be scalar or of size SZ");
  endif

  if (isa (location, "single") || isa (scale, "single") || isa (shape, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (location) && isscalar (scale) && isscalar (shape))
    if ((-Inf < location) && (location < Inf) ...
        && (0 < scale) && (scale < Inf) ...
        && (0 < shape) && (shape < Inf))
      rnd = rand(sz,cls);
      y = shape * norminv (rnd);
      rnd = location + scale * (y + sqrt (4 + y.^2)).^2 / 4;
    else
      rnd = NaN (sz, cls);
    endif
  else
    rnd = NaN (sz, cls);

    k = (-Inf < location) & (location < Inf) ...
        & (0 < scale) & (scale < Inf) ...
        & (0 < shape) & (shape < Inf);
    rnd(k) = rand(sum(k(:)),1);
    y = shape(k) .* norminv (rnd(k));
    rnd(k) = location(k) + scale(k) .* (y + sqrt (4 + y.^2)).^2 / 4;
  endif
endfunction


%!assert (size (bbsrnd (0, 1, 1)), [1 1])
%!assert (size (bbsrnd (zeros (2,1), 1, 1)), [2, 1])
%!assert (size (bbsrnd (zeros (2,2), 1, 1)), [2, 2])
%!assert (size (bbsrnd (0, ones (2,1), 1)), [2, 1])
%!assert (size (bbsrnd (0, ones (2,2), 1)), [2, 2])
%!assert (size (bbsrnd (0, 1, ones (2,1))), [2, 1])
%!assert (size (bbsrnd (0, 1, ones (2,2))), [2, 2])
%!assert (size (bbsrnd (0, 1, 1, 3)), [3, 3])
%!assert (size (bbsrnd (0, 1, 1, [4 1])), [4, 1])
%!assert (size (bbsrnd (0, 1, 1, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (bbsrnd (0,1,1)), "double")
%!assert (class (bbsrnd (single (0),1,1)), "single")
%!assert (class (bbsrnd (single ([0 0]),1,1)), "single")
%!assert (class (bbsrnd (0,single (1),1)), "single")
%!assert (class (bbsrnd (0,single ([1 1]),1)), "single")
%!assert (class (bbsrnd (0,1,single (1))), "single")
%!assert (class (bbsrnd (0,1,single ([1 1]))), "single")

## Test input validation
%!error bbsrnd ()
%!error bbsrnd (1)
%!error bbsrnd (1,2)
%!error bbsrnd (ones (3), ones (2), ones (2), 2)
%!error bbsrnd (ones (2), ones (3), ones (2), 2)
%!error bbsrnd (ones (2), ones (2), ones (3), 2)
%!error bbsrnd (i, 2, 3)
%!error bbsrnd (1, i, 3)
%!error bbsrnd (1, 2, i)
%!error bbsrnd (1,2,3, -1)
%!error bbsrnd (1,2,3, ones (2))
%!error bbsrnd (1,2,3, [2 -1 2])
%!error bbsrnd (ones (2),1,2, 3)
%!error bbsrnd (ones (2),1,2, [3, 2])
%!error bbsrnd (ones (2),1,2, 3, 2)

