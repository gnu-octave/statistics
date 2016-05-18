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
## @deftypefn  {} {} trirnd (@var{a}, @var{b}, @var{c})
## @deftypefnx {} {} trirnd (@var{a}, @var{b}, @var{c}, @var{r})
## @deftypefnx {} {} trirnd (@var{a}, @var{b}, @var{c}, @var{r}, @var{c}, @dots{})
## @deftypefnx {} {} trirnd (@var{a}, @var{b}, @var{c}, [@var{sz}])
## Return a matrix of random samples from the rectangular distribution with
## parameters @var{a}, @var{b}, and @var{c} on the interval [@var{a}, @var{b}].
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## If no size arguments are given then the result matrix is the common size of
## @var{a}, @var{b} and @var{c}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: Random deviates from the triangular distribution

function rnd = trirnd (a, b, c, varargin)

  if (nargin < 3)
    print_usage ();
  endif

  if (! isscalar (a) || ! isscalar (b) || ! isscalar (c))
    [retval, a, b, c] = common_size (a, b, c);
    if (retval > 0)
      error ("trirnd: A, B, and C must be of common size or scalars");
    endif
  endif

  if (iscomplex (a) || iscomplex (b) || iscomplex (c))
    error ("trirnd: A, B, and C must not be complex");
  endif

  if (nargin == 3)
    sz = size (a);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error ("trirnd: dimension vector must be row vector of non-negative integers");
    endif
  elseif (nargin > 4)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("trirnd: dimensions must be non-negative integers");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (a) && ! isequal (size (b), sz))
    error ("trirnd: A, B, and C must be scalar or of size SZ");
  endif

  if (isa (a, "single") || isa (b, "single") || isa (c, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (a) && isscalar (b) && isscalar (c))
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
%!error trirnd ()
%!error trirnd (1)
%!error trirnd (1,2)
%!error trirnd (ones (3), 2*ones (2), 1.5*ones (2), 2)
%!error trirnd (ones (2), 2*ones (3), 1.5*ones (2), 2)
%!error trirnd (ones (2), 2*ones (2), 1.5*ones (3), 2)
%!error trirnd (i, 2, 1.5)
%!error trirnd (1, i, 1.5)
%!error trirnd (1, 2, i)
%!error trirnd (1,2,1.5, -1)
%!error trirnd (1,2,1.5, ones (2))
%!error trirnd (1,2,1.5, [2 -1 2])
%!error trirnd (1*ones (2),2,1.5, 3)
%!error trirnd (1*ones (2),2,1.5, [3, 2])
%!error trirnd (1*ones (2),2,1.5, 3, 2)

