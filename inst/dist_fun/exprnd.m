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
## @deftypefn  {statistics} {@var{r} =} exprnd (@var{mu})
## @deftypefnx {statistics} {@var{r} =} exprnd (@var{mu}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} exprnd (@var{mu}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} exprnd (@var{mu}, [@var{sz}])
##
## Random arrays from the exponential distribution.
##
## @code{@var{r} = exprnd (@var{mu})} returns an array of random numbers chosen
## from the exponential distribution with mean parameter @var{mu}.  The size of
## @var{r} is the size of @var{mu}.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{expcdf, expinv, exppdf, expfit, explike, expstat}
## @end deftypefn

function r = exprnd (mu, varargin)

  if (nargin < 1)
    print_usage ();
  endif

  if (nargin == 1)
    sz = size (mu);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["exprnd: dimension vector must be row vector of"], ...
                     [" non-negative integers."]));
    endif
  elseif (nargin > 2)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("exprnd: dimensions must be non-negative integers");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("exprnd: MU must be scalar or of size SZ.");
  endif

  if (iscomplex (mu))
    error ("exprnd: MU must not be complex.");
  endif

  if (isa (mu, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  if (isscalar (mu))
    if ((mu > 0) && (mu < Inf))
      r = rande (sz, cls) * mu;
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);

    k = (mu > 0) & (mu < Inf);
    r(k) = rande (sum (k(:)), 1, cls) .* mu(k)(:);
  endif

endfunction


%!assert (size (exprnd (2)), [1, 1])
%!assert (size (exprnd (ones (2,1))), [2, 1])
%!assert (size (exprnd (ones (2,2))), [2, 2])
%!assert (size (exprnd (1, 3)), [3, 3])
%!assert (size (exprnd (1, [4 1])), [4, 1])
%!assert (size (exprnd (1, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (exprnd (1)), "double")
%!assert (class (exprnd (single (1))), "single")
%!assert (class (exprnd (single ([1 1]))), "single")

## Test input validation
%!error exprnd ()
%!error exprnd (1, -1)
%!error exprnd (1, ones (2))
%!error exprnd (i)
%!error exprnd (1, [2 -1 2])
%!error exprnd (1, 2, -1)
%!error exprnd (1, 2, ones (2))
%!error exprnd (ones (2,2), 3)
%!error exprnd (ones (2,2), [3, 2])
%!error exprnd (ones (2,2), 2, 3)
