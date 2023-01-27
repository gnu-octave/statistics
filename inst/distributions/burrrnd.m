## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{r} = burrrnd (@var{a}, @var{c}, @var{k})
## @deftypefnx {statistics} @var{r} = burrrnd (@var{a}, @var{c}, @var{k}, @var{m})
## @deftypefnx {statistics} @var{r} = burrrnd (@var{a}, @var{c}, @var{k}, @var{m}, @var{n}, @dots{})
## @deftypefnx {statistics} @var{r} = burrrnd (@var{a}, @var{c}, @var{k}, [@var{sz}])
##
## Random arrays from the Burr type XII distribution.
##
## @code{@var{r} = burrrnd (@var{a}, @var{c}, @var{k})} returns an array of
## random numbers chosen from the Burr type XII distribution with parameters
## @var{a}, @var{c}, and @var{k}.  The size of @var{r} is the common size of
## @var{a}, @var{c}, and @var{k}.  A scalar input functions as a constant matrix
## of the same size as the other inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## If no size arguments are given then the result matrix is the common size of
## @var{a}, @var{c} and @var{k}.
##
## @seealso{burrcdf, burrinv, burrpdf}
## @end deftypefn

function r = burrrnd (a, c, k, varargin)

  if (nargin < 3)
    print_usage ();
  endif

  if (! isscalar (a) || ! isscalar (c) || ! isscalar (k))
    [retval, a, c, k] = common_size (a, c, k);
    if (retval > 0)
      error ("burrrnd: A, C, and K must be of common size or scalars.");
    endif
  endif

  if (iscomplex (a) || iscomplex (c) || iscomplex (k))
    error ("burrrnd: A, C, and K must not be complex.");
  endif

  if (nargin == 3)
    sz = size (a);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["burrrnd: dimension vector must be row vector of"], ...
                     [" non-negative integers."]));
    endif
  elseif (nargin > 4)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("burrrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (a) && ! isequal (size (c), sz) && ! isequal (size (k), sz))
    error ("burrrnd: A, C, and K must be scalar or of size SZ.");
  endif

  if (isa (a, "single") || isa (c, "single") || isa (k, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (a) && isscalar (c) && isscalar(k))
    if ((0 < a) && (a < Inf) && (0 < c) && (c < Inf) ...
        && (0 < k) && (k < Inf))
      r = rand (sz, cls);
      r(:) = ((1 - r(:) / a).^(-1 / k) - 1).^(1 / c);
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    j = (0 < a) && (a < Inf) && (0 < c) && (c < Inf) ...
        && (0 < k) && (k < Inf);
    r(k) = rand(sum(j(:)),1);
    r(k) = ((1 - r(j) / a(j)).^(-1 ./ k(j)) - 1).^(1 ./ c(j));
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

