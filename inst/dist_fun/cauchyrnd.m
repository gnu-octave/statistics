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
## @deftypefn  {statistics} {@var{r} =} cauchy_rnd (@var{location}, @var{scale})
## @deftypefnx {statistics} {@var{r} =} cauchy_rnd (@var{location}, @var{scale}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} cauchy_rnd (@var{location}, @var{scale}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} cauchy_rnd (@var{location}, @var{scale}, [@var{sz}])
##
## Random arrays from the Cauchy distribution.
##
## @code{@var{r} = cauchy_rnd (@var{location}, @var{scale})} returns an array of
## random numbers chosen from the Cauchy distribution with parameters
## @var{location} and @var{scale}.  The size of @var{r} is the common size of
## @var{location} and @var{scale}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{cauchy_cdf, cauchy_inv, cauchy_pdf}
## @end deftypefn

function r = cauchy_rnd (location, scale, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale))
    [retval, location, scale] = common_size (location, scale);
    if (retval > 0)
      error (strcat (["cauchy_rnd: LOCATION and SCALE must be of common"], ...
                     [" size or scalars."]));
    endif
  endif

  if (iscomplex (location) || iscomplex (scale))
    error ("cauchy_rnd: LOCATION and SCALE must not be complex.");
  endif

  if (nargin == 2)
    sz = size (location);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["cauchy_rnd: dimension vector must be a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("cauchy_rnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (location) && ! isequal (size (location), sz))
    error ("cauchy_rnd: LOCATION and SCALE must be scalar or of size SZ.");
  endif

  if (isa (location, "single") || isa (scale, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (location) && isscalar (scale))
    if (! isinf (location) && (scale > 0) && (scale < Inf))
      r = location - cot (pi * rand (sz, cls)) * scale;
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    k = ! isinf (location) & (scale > 0) & (scale < Inf);
    r(k) = location(k)(:) ...
             - cot (pi * rand (sum (k(:)), 1, cls)) .* scale(k)(:);
  endif

endfunction


%!assert (size (cauchy_rnd (1,2)), [1, 1])
%!assert (size (cauchy_rnd (ones (2,1), 2)), [2, 1])
%!assert (size (cauchy_rnd (ones (2,2), 2)), [2, 2])
%!assert (size (cauchy_rnd (1, 2*ones (2,1))), [2, 1])
%!assert (size (cauchy_rnd (1, 2*ones (2,2))), [2, 2])
%!assert (size (cauchy_rnd (1, 2, 3)), [3, 3])
%!assert (size (cauchy_rnd (1, 2, [4 1])), [4, 1])
%!assert (size (cauchy_rnd (1, 2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (cauchy_rnd (1, 2)), "double")
%!assert (class (cauchy_rnd (single (1), 2)), "single")
%!assert (class (cauchy_rnd (single ([1 1]), 2)), "single")
%!assert (class (cauchy_rnd (1, single (2))), "single")
%!assert (class (cauchy_rnd (1, single ([2 2]))), "single")

## Test input validation
%!error cauchy_rnd ()
%!error cauchy_rnd (1)
%!error cauchy_rnd (ones (3), ones (2))
%!error cauchy_rnd (ones (2), ones (3))
%!error cauchy_rnd (i, 2)
%!error cauchy_rnd (2, i)
%!error cauchy_rnd (1,2, -1)
%!error cauchy_rnd (1,2, ones (2))
%!error cauchy_rnd (1,2, [2 -1 2])
%!error cauchy_rnd (1,2, 1, ones (2))
%!error cauchy_rnd (1,2, 1, -1)
%!error cauchy_rnd (ones (2,2), 2, 3)
%!error cauchy_rnd (ones (2,2), 2, [3, 2])
%!error cauchy_rnd (ones (2,2), 2, 2, 3)
