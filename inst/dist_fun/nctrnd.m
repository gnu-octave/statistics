## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} nctrnd (@var{df}, @var{mu})
## @deftypefnx {statistics} {@var{r} =} nctrnd (@var{df}, @var{mu}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} nctrnd (@var{df}, @var{mu}, [@var{sz}])
##
## Random arrays from the noncentral T distribution.
##
## @code{@var{x} = nctrnd (@var{p}, @var{df}, @var{mu})} returns an array of
## random numbers chosen from the noncentral T distribution with @var{df}
## degrees of freedom and noncentrality parameter @var{mu}.  The size of @var{r}
## is the common size of @var{df} and @var{mu}.  A scalar input functions as a
## constant matrix of the same size as the other input.
##
## @code{nctrnd} generates values using the definition of a noncentral T random
## variable, as the ratio of a normal with non-zero mean and the sqrt of a
## chi-square.
##
## When called with a single size argument, @code{nctrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the noncentral T distribution can be found at
## @url{https://en.wikipedia.org/wiki/Noncentral_t-distribution}
##
## @seealso{nctcdf, nctinv, nctpdf, nctstat, trnd, normrnd, chi2rnd}
## @end deftypefn

function r = nctrnd (df, mu, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("nctrnd: function called with too few input arguments.");
  endif

  ## Check for common size of DF and MU
  if (! isscalar (df) || ! isscalar (mu))
    [retval, df, mu] = common_size (df, mu);
    if (retval > 0)
      error ("nctrnd: DF and MU must be of common size or scalars.");
    endif
  endif

  ## Check for DF and MU being reals
  if (iscomplex (df) || iscomplex (mu))
    error ("nctrnd: DF and MU must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (df);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["nctrnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("nctrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (df) && ! isequal (size (df), sz))
    error ("nctrnd: DF and MU must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (df, "single") || isa (mu, "single"));
    cls = "single";
  else
    cls = "double";
  endif

  ## Return NaNs for out of range values of DF
  df(df <= 0) = NaN;

  ## Prevent Inf/Inf==NaN for the standardized chi-square in the denom.
  df(isinf (df)) = realmax;

  ## Generate random sample from noncentral F distribution
  r = (randn (sz) + mu) ./ sqrt (2 .* randg (df ./ 2, sz) ./ df);

  ## Cast to appropriate class
  r = cast (r, cls);

endfunction

## Test output
%!assert (size (nctrnd (1, 1)), [1 1])
%!assert (size (nctrnd (1, ones (2,1))), [2, 1])
%!assert (size (nctrnd (1, ones (2,2))), [2, 2])
%!assert (size (nctrnd (ones (2,1), 1)), [2, 1])
%!assert (size (nctrnd (ones (2,2), 1)), [2, 2])
%!assert (size (nctrnd (1, 1, 3)), [3, 3])
%!assert (size (nctrnd (1, 1, [4, 1])), [4, 1])
%!assert (size (nctrnd (1, 1, 4, 1)), [4, 1])
%!assert (size (nctrnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (nctrnd (1, 1, 0, 1)), [0, 1])
%!assert (size (nctrnd (1, 1, 1, 0)), [1, 0])
%!assert (size (nctrnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (nctrnd (1, 1)), "double")
%!assert (class (nctrnd (1, single (1))), "single")
%!assert (class (nctrnd (1, single ([1, 1]))), "single")
%!assert (class (nctrnd (single (1), 1)), "single")
%!assert (class (nctrnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<nctrnd: function called with too few input arguments.> nctrnd ()
%!error<nctrnd: function called with too few input arguments.> nctrnd (1)
%!error<nctrnd: DF and MU must be of common size or scalars.> ...
%! nctrnd (ones (3), ones (2))
%!error<nctrnd: DF and MU must be of common size or scalars.> ...
%! nctrnd (ones (2), ones (3))
%!error<nctrnd: DF and MU must not be complex.> nctrnd (i, 2)
%!error<nctrnd: DF and MU must not be complex.> nctrnd (1, i)
%!error<nctrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nctrnd (1, 2, -1)
%!error<nctrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nctrnd (1, 2, 1.2)
%!error<nctrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nctrnd (1, 2, ones (2))
%!error<nctrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nctrnd (1, 2, [2 -1 2])
%!error<nctrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nctrnd (1, 2, [2 0 2.5])
%!error<nctrnd: dimensions must be non-negative integers.> ...
%! nctrnd (1, 2, 2, -1, 5)
%!error<nctrnd: dimensions must be non-negative integers.> ...
%! nctrnd (1, 2, 2, 1.5, 5)
%!error<nctrnd: DF and MU must be scalars or of size SZ.> ...
%! nctrnd (2, ones (2), 3)
%!error<nctrnd: DF and MU must be scalars or of size SZ.> ...
%! nctrnd (2, ones (2), [3, 2])
%!error<nctrnd: DF and MU must be scalars or of size SZ.> ...
%! nctrnd (2, ones (2), 3, 2)
