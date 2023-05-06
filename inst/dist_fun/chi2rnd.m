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
## @deftypefn  {statistics} {@var{r} =} chi2rnd (@var{df})
## @deftypefnx {statistics} {@var{r} =} chi2rnd (@var{df}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} chi2rnd (@var{df}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} chi2rnd (@var{df}, [@var{sz}])
##
## Random arrays from the chi-squared distribution.
##
## @code{@var{r} = chi2rnd (@var{df})} returns an array of random numbers chosen
## from the chi-squared distribution with @var{df} degrees of freedom.  The size
## of @var{r} is the size of @var{df}.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## Further information about the chi-squared distribution can be found at
## @url{https://en.wikipedia.org/wiki/Chi-squared_distribution}
##
## @seealso{chi2cdf, chi2inv, chi2pdf, chi2stat}
## @end deftypefn

function r = chi2rnd (df, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("chi2rnd: function called with too few input arguments.");
  endif

  ## Check for DF being reals
  if (iscomplex (df))
    error ("chi2rnd: DF must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 1)
    sz = size (df);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["chi2rnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 2)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("chi2rnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameter match requested dimensions in size
  if (! isscalar (df) && ! isequal (size (df), sz))
    error ("chi2rnd: DF must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (df, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from chi-squared distribution
  if (isscalar (df))
    if ((df > 0) && (df < Inf))
      r = 2 * randg (df/2, sz, cls);
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);
    k = (df > 0) | (df < Inf);
    r(k) = 2 * randg (df(k)/2, cls);
  endif

endfunction

## Test output
%!assert (size (chi2rnd (2)), [1, 1])
%!assert (size (chi2rnd (ones (2,1))), [2, 1])
%!assert (size (chi2rnd (ones (2,2))), [2, 2])
%!assert (size (chi2rnd (1, 3)), [3, 3])
%!assert (size (chi2rnd (1, [4 1])), [4, 1])
%!assert (size (chi2rnd (1, 4, 1)), [4, 1])
%!assert (size (chi2rnd (1, 4, 1)), [4, 1])
%!assert (size (chi2rnd (1, 4, 1, 5)), [4, 1, 5])
%!assert (size (chi2rnd (1, 0, 1)), [0, 1])
%!assert (size (chi2rnd (1, 1, 0)), [1, 0])
%!assert (size (chi2rnd (1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (chi2rnd (2)), "double")
%!assert (class (chi2rnd (single (2))), "single")
%!assert (class (chi2rnd (single ([2 2]))), "single")

## Test input validation
%!error<chi2rnd: function called with too few input arguments.> chi2rnd ()
%!error<chi2rnd: DF must not be complex.> chi2rnd (i)
%!error<chi2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! chi2rnd (1, -1)
%!error<chi2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! chi2rnd (1, 1.2)
%!error<chi2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! chi2rnd (1, ones (2))
%!error<chi2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! chi2rnd (1, [2 -1 2])
%!error<chi2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! chi2rnd (1, [2 0 2.5])
%!error<chi2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! chi2rnd (ones (2), ones (2))
%!error<chi2rnd: dimensions must be non-negative integers.> ...
%! chi2rnd (1, 2, -1, 5)
%!error<chi2rnd: dimensions must be non-negative integers.> ...
%! chi2rnd (1, 2, 1.5, 5)
%!error<chi2rnd: DF must be scalar or of size SZ.> chi2rnd (ones (2,2), 3)
%!error<chi2rnd: DF must be scalar or of size SZ.> chi2rnd (ones (2,2), [3, 2])
%!error<chi2rnd: DF must be scalar or of size SZ.> chi2rnd (ones (2,2), 2, 3)
