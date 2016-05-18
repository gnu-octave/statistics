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
## @deftypefn  {} {} gprnd (@var{location}, @var{scale}, @var{shape})
## @deftypefnx {} {} gprnd (@var{location}, @var{scale}, @var{shape}, @var{r})
## @deftypefnx {} {} gprnd (@var{location}, @var{scale}, @var{shape}, @var{r}, @var{c}, @dots{})
## @deftypefnx {} {} gprnd (@var{location}, @var{scale}, @var{shape}, [@var{sz}])
## Return a matrix of random samples from the generalized Pareto distribution
## with parameters @var{location}, @var{scale} and @var{shape}.
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
## Description: Random deviates from the generalized Pareto distribution

function rnd = gprnd (location, scale, shape, varargin)

  if (nargin < 3)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale) || ! isscalar (shape))
    [retval, location, scale, shape] = common_size (location, scale, shape);
    if (retval > 0)
      error ("gpgrnd: LOCATION, SCALE and SHAPE must be of common size or scalars");
    endif
  endif

  if (iscomplex (location) || iscomplex (scale) || iscomplex (shape))
    error ("gprnd: LOCATION, SCALE and SHAPE must not be complex");
  endif

  if (nargin == 3)
    sz = size (location);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error ("gprnd: dimension vector must be row vector of non-negative integers");
    endif
  elseif (nargin > 4)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("gprnd: dimensions must be non-negative integers");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (location) && ! isequal (size (location), sz))
    error ("gprnd: LOCATION, SCALE and SHAPE must be scalar or of size SZ");
  endif

  if (isa (location, "single") || isa (scale, "single") || isa (shape, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (location) && isscalar (scale) && isscalar (shape))
    if ((-Inf < location) && (location < Inf) && (0 < scale) && (scale < Inf) ...
          && (-Inf < shape) && (shape < Inf))
      rnd = rand(sz,cls);
      if (shape == 0)
        rnd = -log(1 - rnd);
        rnd = scale * rnd + location;
      elseif ((shape < 0) || (shape > 0))
        rnd = (1 - rnd).^(-shape) - 1;
        rnd = (scale / shape) * rnd + location;
      end
    else
      rnd = NaN (sz, cls);
    endif
  else
    rnd = NaN (sz, cls);
    
    k = (-Inf < location) & (location < Inf) & (scale > 0) ...
        & (-Inf < shape) & (shape < Inf);
    rnd(k(:)) = rand (1, sum(k(:)), cls);
    if (any (shape == 0))
        rnd(k) = -log(1 - rnd(k));
        rnd(k) = scale(k) .* rnd(k) + location(k);
    elseif (any (shape < 0 | shape > 0))
      rnd(k) = (1 - rnd(k)).^(-shape(k)) - 1;
      rnd(k) = (scale(k) ./ shape(k)) .* rnd(k) + location(k);
    end
  endif
endfunction


%!assert (size (gprnd (0,1,0)), [1, 1])
%!assert (size (gprnd (zeros (2,1), 1,0)), [2, 1])
%!assert (size (gprnd (zeros (2,2), 1,0)), [2, 2])
%!assert (size (gprnd (0, ones (2,1), 0)), [2, 1])
%!assert (size (gprnd (0, ones (2,2), 0)), [2, 2])
%!assert (size (gprnd (0,1, zeros (2,1))), [2, 1])
%!assert (size (gprnd (0,1, zeros (2,2))), [2, 2])
%!assert (size (gprnd (0,1, 0, 3)), [3, 3])
%!assert (size (gprnd (0,1, 0, [4 1])), [4, 1])
%!assert (size (gprnd (0,1, 0, 4, 1)), [4, 1])

%!assert (size (gprnd (0,1,1)), [1, 1])
%!assert (size (gprnd (zeros (2,1), 1,1)), [2, 1])
%!assert (size (gprnd (zeros (2,2), 1,1)), [2, 2])
%!assert (size (gprnd (0, ones (2,1), 1)), [2, 1])
%!assert (size (gprnd (0, ones (2,2), 1)), [2, 2])
%!assert (size (gprnd (0,1, ones (2,1))), [2, 1])
%!assert (size (gprnd (0,1, ones (2,2))), [2, 2])
%!assert (size (gprnd (0,1, 1, 3)), [3, 3])
%!assert (size (gprnd (0,1, 1, [4 1])), [4, 1])
%!assert (size (gprnd (0,1, 1, 4, 1)), [4, 1])

%!assert (size (gprnd (0,1,-1)), [1, 1])
%!assert (size (gprnd (zeros (2,1), -1,1)), [2, 1])
%!assert (size (gprnd (zeros (2,2), -1,1)), [2, 2])
%!assert (size (gprnd (0, ones (2,1), -1)), [2, 1])
%!assert (size (gprnd (0, ones (2,2), -1)), [2, 2])
%!assert (size (gprnd (0,1, -ones (2,1))), [2, 1])
%!assert (size (gprnd (0,1, -ones (2,2))), [2, 2])
%!assert (size (gprnd (0,1, -1, 3)), [3, 3])
%!assert (size (gprnd (0,1, -1, [4 1])), [4, 1])
%!assert (size (gprnd (0,1, -1, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (gprnd (0,1,0)), "double")
%!assert (class (gprnd (single (0),1,0)), "single")
%!assert (class (gprnd (single ([0 0]),1,0)), "single")
%!assert (class (gprnd (0,single (1),0)), "single")
%!assert (class (gprnd (0,single ([1 1]),0)), "single")
%!assert (class (gprnd (0,1,single (0))), "single")
%!assert (class (gprnd (0,1,single ([0 0]))), "single")

## Test input validation
%!error gprnd ()
%!error gprnd (1)
%!error gprnd (1,2)
%!error gprnd (zeros (3), ones (2), zeros (2))
%!error gprnd (zeros (2), ones (3), zeros (2))
%!error gprnd (zeros (2), ones (2), zeros (3))
%!error gprnd (i, 1, 0)
%!error gprnd (0, i, 0)
%!error gprnd (0, 1, i)
%!error gprnd (0,1,0, -1)
%!error gprnd (0,1,0, ones (2))
%!error gprnd (0,1,0, [2 -1 2])
%!error gprnd (zeros (2),1,0, 3)
%!error gprnd (zeros (2),1,0, [3, 2])
%!error gprnd (zeros (2),1,0, 3, 2)

