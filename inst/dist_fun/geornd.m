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
## @deftypefn  {statistics} {@var{r} =} geornd (@var{ps})
## @deftypefnx {statistics} {@var{r} =} geornd (@var{ps}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} geornd (@var{ps}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} geornd (@var{ps}, [@var{sz}])
##
## Random arrays from the geometric distribution.
##
## @code{@var{r} = geornd (@var{ps})} returns an array of random numbers chosen
## from the Birnbaum-Saunders distribution with probability of success parameter
## @var{ps}.  The size of @var{r} is the size of @var{ps}.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## The geometric distribution models the number of failures (@var{x}) of a
## Bernoulli trial with probability @var{ps} before the first success.
##
## Further information about the geometric distribution can be found at
## @url{https://en.wikipedia.org/wiki/Geometric_distribution}
##
## @seealso{geocdf, geoinv, geopdf, geostat}
## @end deftypefn

function r = geornd (ps, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("geornd: function called with too few input arguments.");
  endif

  ## Check for PS being reals
  if (iscomplex (ps))
    error ("geornd: PS must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 1)
    sz = size (ps);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["geornd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 2)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("geornd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameter match requested dimensions in size
  if (! isscalar (ps) && ! isequal (size (ps), sz))
    error ("geornd: PS must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (ps, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from geometric distribution
  if (isscalar (ps))
    if (ps > 0 && ps < 1);
      r = floor (- rande (sz, cls) ./ log (1 - ps));
    elseif (ps == 0)
      r = Inf (sz, cls);
    elseif (ps == 1)
      r = zeros (sz, cls);
    elseif (ps < 0 || ps > 1)
      r = NaN (sz, cls);
    endif
  else
    r = floor (- rande (sz, cls) ./ log (1 - ps));

    k = ! (ps >= 0) | ! (ps <= 1);
    r(k) = NaN;

    k = (ps == 0);
    r(k) = Inf;
  endif

endfunction

## Test output
%!assert (size (geornd (0.5)), [1, 1])
%!assert (size (geornd (0.5*ones (2,1))), [2, 1])
%!assert (size (geornd (0.5*ones (2,2))), [2, 2])
%!assert (size (geornd (0.5, 3)), [3, 3])
%!assert (size (geornd (0.5, [4 1])), [4, 1])
%!assert (size (geornd (0.5, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (geornd (0.5)), "double")
%!assert (class (geornd (single (0.5))), "single")
%!assert (class (geornd (single ([0.5 0.5]))), "single")
%!assert (class (geornd (single (0))), "single")
%!assert (class (geornd (single (1))), "single")

## Test input validation
%!error<geornd: function called with too few input arguments.> geornd ()
%!error<geornd: DF must not be complex.> geornd (i)
%!error<geornd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! geornd (1, -1)
%!error<geornd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! geornd (1, 1.2)
%!error<geornd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! geornd (1, ones (2))
%!error<geornd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! geornd (1, [2 -1 2])
%!error<geornd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! geornd (1, [2 0 2.5])
%!error<geornd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! geornd (ones (2), ones (2))
%!error<geornd: dimensions must be non-negative integers.> ...
%! geornd (1, 2, -1, 5)
%!error<geornd: dimensions must be non-negative integers.> ...
%! geornd (1, 2, 1.5, 5)
%!error<geornd: DF must be scalar or of size SZ.> geornd (ones (2,2), 3)
%!error<geornd: DF must be scalar or of size SZ.> geornd (ones (2,2), [3, 2])
%!error<geornd: DF must be scalar or of size SZ.> geornd (ones (2,2), 2, 3)
