## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 2005-2016 John W. Eaton
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
## @deftypefn  {statistics} {@var{r} =} unidrnd (@var{df})
## @deftypefnx {statistics} {@var{r} =} unidrnd (@var{df}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} unidrnd (@var{df}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} unidrnd (@var{df}, [@var{sz}])
##
## Random arrays from the discrete uniform distribution.
##
## @code{@var{r} = unidrnd (@var{df})} returns an array of random numbers chosen
## from the discrete uniform distribution with @var{df} degrees of freedom.  The
## size of @var{r} is the size of @var{df}.  @var{sigma} must be a finite
## integer greater than 0, otherwise NaN is returned.
##
## @var{df} may be a scalar or a multi-dimensional array.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{unidcdf, unidinv, unidrnd, unidstat}
## @end deftypefn

function r = unidrnd (df, varargin)

  if (nargin < 1)
    print_usage ();
  endif

  if (nargin == 1)
    sz = size (df);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["unidrnd: dimension vector must be a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 2)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("unidrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (df) && ! isequal (size (df), sz))
    error ("unidrnd: DF must be scalar or of size SZ.");
  endif

  if (iscomplex (df))
    error ("unidrnd: DF must not be complex.");
  endif

  if (isa (df, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (df))
    if (df > 0 && df == fix (df))
      r = ceil (rand (sz, cls) * df);
    else
      r = NaN (sz, cls);
    endif
  else
    r = ceil (rand (sz, cls) .* df);

    k = ! (df > 0 & df == fix (df));
    r(k) = NaN;
  endif

endfunction


%!assert (size (unidrnd (2)), [1, 1])
%!assert (size (unidrnd (ones (2,1))), [2, 1])
%!assert (size (unidrnd (ones (2,2))), [2, 2])
%!assert (size (unidrnd (10, [4 1])), [4, 1])
%!assert (size (unidrnd (10, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (unidrnd (2)), "double")
%!assert (class (unidrnd (single (2))), "single")
%!assert (class (unidrnd (single ([2 2]))), "single")

## Test input validation
%!error unidrnd ()
%!error unidrnd (10, [1;2;3])
%!error unidrnd (10, 2, ones (2))
%!error unidrnd (10*ones (2), 2, 1)
%!error unidrnd (i)
