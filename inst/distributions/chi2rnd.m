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
## Random arrays from the Chi-square distribution.
##
## @code{@var{r} = chi2rnd (@var{df})} returns an array of random numbers chosen
## from the Chi-square distribution with @var{df} degrees of freedom.  The size
## of @var{r} is the size of @var{df}.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{chi2cdf, chi2inv, chi2pdf, chi2stat}
## @end deftypefn

function r = chi2rnd (df, varargin)

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
      error (strcat (["chi2rnd: dimension vector must be row vector of"], ...
                     [" non-negative integers."]));
    endif
  elseif (nargin > 2)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("chi2rnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (df) && ! isequal (size (df), sz))
    error ("chi2rnd: DF must be scalar or of size SZ.");
  endif

  if (iscomplex (df))
    error ("chi2rnd: DF must not be complex.");
  endif

  if (isa (df, "single"))
    cls = "single";
  else
    cls = "double";
  endif

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


%!assert (size (chi2rnd (2)), [1, 1])
%!assert (size (chi2rnd (ones (2,1))), [2, 1])
%!assert (size (chi2rnd (ones (2,2))), [2, 2])
%!assert (size (chi2rnd (1, 3)), [3, 3])
%!assert (size (chi2rnd (1, [4 1])), [4, 1])
%!assert (size (chi2rnd (1, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (chi2rnd (2)), "double")
%!assert (class (chi2rnd (single (2))), "single")
%!assert (class (chi2rnd (single ([2 2]))), "single")

## Test input validation
%!error chi2rnd ()
%!error chi2rnd (ones (3), ones (2))
%!error chi2rnd (ones (2), ones (3))
%!error chi2rnd (i)
%!error chi2rnd (1, -1)
%!error chi2rnd (1, ones (2))
%!error chi2rnd (1, [2 -1 2])
%!error chi2rnd (ones (2,2), 3)
%!error chi2rnd (ones (2,2), [3, 2])
%!error chi2rnd (ones (2,2), 2, 3)
