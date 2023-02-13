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
## @deftypefn  {statistics} {@var{r} =} gamrnd (@var{a}, @var{b})
## @deftypefnx {statistics} {@var{r} =} gamrnd (@var{a}, @var{b}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} gamrnd (@var{a}, @var{b}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} gamrnd (@var{a}, @var{b}, [@var{sz}])
##
## Random arrays from the Gamma distribution.
##
## @code{@var{r} = gamrnd (@var{a}, @var{b})} returns an array of random numbers
## chosen from the Gamma distribution with shape parameter @var{a} and scale
## parameter @var{b}.  The size of @var{r} is the common size of @var{a} and
## @var{b}.  A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{gamcdf, gaminv, gampdf, gamfit, gamlike, gamstat}
## @end deftypefn

function r = gamrnd (a, b, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  if (! isscalar (a) || ! isscalar (b))
    [retval, a, b] = common_size (a, b);
    if (retval > 0)
      error ("gamrnd: A and B must be of common size or scalars.");
    endif
  endif

  if (iscomplex (a) || iscomplex (b))
    error ("gamrnd: A and B must not be complex.");
  endif

  if (nargin == 2)
    sz = size (a);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["gamrnd: dimension vector must be row vector of"], ...
                     [" non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("gamrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (a) && ! isequal (size (a), sz))
    error ("gamrnd: A and B must be scalar or of size SZ.");
  endif

  if (isa (a, "single") || isa (b, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (a) && isscalar (b))
    if ((a > 0) && (a < Inf) && (b > 0) && (b < Inf))
      r = b * randg (a, sz, cls);
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    k = (a > 0) & (a < Inf) & (b > 0) & (b < Inf);
    r(k) = b(k) .* randg (a(k), cls);
  endif

endfunction


%!assert (size (gamrnd (1,2)), [1, 1])
%!assert (size (gamrnd (ones (2,1), 2)), [2, 1])
%!assert (size (gamrnd (ones (2,2), 2)), [2, 2])
%!assert (size (gamrnd (1, 2*ones (2,1))), [2, 1])
%!assert (size (gamrnd (1, 2*ones (2,2))), [2, 2])
%!assert (size (gamrnd (1, 2, 3)), [3, 3])
%!assert (size (gamrnd (1, 2, [4 1])), [4, 1])
%!assert (size (gamrnd (1, 2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (gamrnd (1, 2)), "double")
%!assert (class (gamrnd (single (1), 2)), "single")
%!assert (class (gamrnd (single ([1 1]), 2)), "single")
%!assert (class (gamrnd (1, single (2))), "single")
%!assert (class (gamrnd (1, single ([2 2]))), "single")

## Test input validation
%!error gamrnd ()
%!error gamrnd (1)
%!error gamrnd (ones (3), ones (2))
%!error gamrnd (ones (2), ones (3))
%!error gamrnd (i, 2)
%!error gamrnd (2, i)
%!error gamrnd (1,2, -1)
%!error gamrnd (1,2, ones (2))
%!error gamrnd (1, 2, [2 -1 2])
%!error gamrnd (1,2, 1, ones (2))
%!error gamrnd (1,2, 1, -1)
%!error gamrnd (ones (2,2), 2, 3)
%!error gamrnd (ones (2,2), 2, [3, 2])
%!error gamrnd (ones (2,2), 2, 2, 3)
