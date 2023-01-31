## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2019 Anthony Morast
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
## @deftypefn  {statistics} @var{r} = unifrnd (@var{a}, @var{b})
## @deftypefnx {statistics} @var{r} = unifrnd (@var{a}, @var{b}, @var{rows})
## @deftypefnx {statistics} @var{r} = unifrnd (@var{a}, @var{b}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} @var{r} = unifrnd (@var{a}, @var{b}, [@var{sz}])
##
## Random arrays from the uniform distribution.
##
## @code{@var{r} = unifrnd (@var{a}, @var{b})} returns an array of random
## numbers chosen from a uniform distribution on the interval [@var{a},
## @var{b}].  The size of @var{y} is the common size of the input arguments.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{unifcdf, unifinv, unifpdf, unifstat}
## @end deftypefn

function r = unifrnd (a, b, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  if (! isscalar (a) || ! isscalar (b))
    [retval, a, b] = common_size (a, b);
    if (retval > 0)
      error ("unifrnd: A and B must be of common size or scalars.");
    endif
  endif

  if (iscomplex (a) || iscomplex (b))
    error ("unifrnd: A and B must not be complex.");
  endif

  if (nargin == 2)
    sz = size (a);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["unifrnd: dimension vector must be row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("unifrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (a) && ! isequal (size (a), sz))
    error ("unifrnd: A and B must be scalar or of size SZ.");
  endif

  if (isa (a, "single") || isa (b, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (a) && isscalar (b))
    if ((-Inf < a) && (a <= b) && (b < Inf))
      r = a + (b - a) * rand (sz, cls);
    else
      r = NaN (sz, cls);
    endif
  else
    r = a + (b - a) .* rand (sz, cls);

    k = !(-Inf < a) | !(a <= b) | !(b < Inf);
    r(k) = NaN;
  endif

endfunction


%!assert (size (unifrnd (1,2)), [1, 1])
%!assert (size (unifrnd (ones (2,1), 2)), [2, 1])
%!assert (size (unifrnd (ones (2,2), 2)), [2, 2])
%!assert (size (unifrnd (1, 2*ones (2,1))), [2, 1])
%!assert (size (unifrnd (1, 2*ones (2,2))), [2, 2])
%!assert (size (unifrnd (1, 2, 3)), [3, 3])
%!assert (size (unifrnd (1, 2, [4 1])), [4, 1])
%!assert (size (unifrnd (1, 2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (unifrnd (1, 2)), "double")
%!assert (class (unifrnd (single (1), 2)), "single")
%!assert (class (unifrnd (single ([1 1]), 2)), "single")
%!assert (class (unifrnd (1, single (2))), "single")
%!assert (class (unifrnd (1, single ([2 2]))), "single")

## Test input validation
%!error unifrnd ()
%!error unifrnd (1)
%!error unifrnd (ones (3), ones (2))
%!error unifrnd (ones (2), ones (3))
%!error unifrnd (i, 2)
%!error unifrnd (2, i)
%!error unifrnd (1,2, -1)
%!error unifrnd (1,2, ones (2))
%!error unifrnd (1, 2, [2 -1 2])
%!error unifrnd (1,2, 1, ones (2))
%!error unifrnd (1,2, 1, -1)
%!error unifrnd (ones (2,2), 2, 3)
%!error unifrnd (ones (2,2), 2, [3, 2])
%!error unifrnd (ones (2,2), 2, 2, 3)

%!assert (unifrnd (0,0), 0)
%!assert (unifrnd (1,1), 1)
%!assert (unifrnd (1,0), NaN)

