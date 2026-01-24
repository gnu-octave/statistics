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
## @deftypefn  {statistics} {@var{r} =} ncx2rnd (@var{df}, @var{lambda})
## @deftypefnx {statistics} {@var{r} =} ncx2rnd (@var{df}, @var{lambda}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} ncx2rnd (@var{df}, @var{lambda}, [@var{sz}])
##
## Random arrays from the noncentral chi-squared distribution.
##
## @code{@var{r} = ncx2rnd (@var{df}, @var{lambda})} returns an array of random
## numbers chosen from the noncentral chi-squared distribution with @var{df}
## degrees of freedom and noncentrality parameter @var{lambda}.  The size of
## @var{r} is the common size of @var{df} and @var{lambda}.  A scalar input
## functions as a constant matrix of the same size as the other input.
##
## When called with a single size argument, @code{ncx2rnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the noncentral chi-squared distribution can be
## found at @url{https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution}
##
## @seealso{ncx2cdf, ncx2inv, ncx2pdf, ncx2stat}
## @end deftypefn

function r = ncx2rnd (df, lambda, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("ncx2rnd: function called with too few input arguments.");
  endif

  ## Check for common size of DF and LAMBDA
  if (! isscalar (df) || ! isscalar (lambda))
    [retval, df, lambda] = common_size (df, lambda);
    if (retval > 0)
      error ("ncx2rnd: DF and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Check for DF and LAMBDA being reals
  if (iscomplex (df) || iscomplex (lambda))
    error ("ncx2rnd: DF and LAMBDA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (df);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif ((isrow (varargin{1}) || isempty (varargin{1})) &&
            all (varargin{1} >= 0) && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("ncx2rnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("ncx2rnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (df) && ! isequal (size (df), sz))
    error ("ncx2rnd: DF and LAMBDA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (df, "single") || isa (lambda, "single"));
    cls = "single";
  else
    cls = "double";
  endif

  ## Return NaNs for out of range values of DF and LAMBDA
  df(df <= 0) = NaN;
  lambda(lambda <= 0) = NaN;

  ## Force DF and LAMBDA into the same size as SZ (if necessary)
  if (isscalar (df))
    df = repmat (df, sz);
  endif
  if (isscalar (lambda))
    lambda = repmat (lambda, sz);
  endif

  ## Generate random sample from noncentral chi-squared distribution
  r = randp (lambda ./ 2);
  r(r > 0) = 2 * randg (r(r > 0));
  r(df > 0) += 2 * randg (df(df > 0) / 2);

  ## Cast to appropriate class
  r = cast (r, cls);

endfunction

## Test output
%!assert (size (ncx2rnd (1, 1)), [1, 1])
%!assert (size (ncx2rnd (1, ones (2, 1))), [2, 1])
%!assert (size (ncx2rnd (1, ones (2, 2))), [2, 2])
%!assert (size (ncx2rnd (ones (2, 1), 1)), [2, 1])
%!assert (size (ncx2rnd (ones (2, 2), 1)), [2, 2])
%!assert (size (ncx2rnd (1, 1, 3)), [3, 3])
%!assert (size (ncx2rnd (1, 1, [4, 1])), [4, 1])
%!assert (size (ncx2rnd (1, 1, 4, 1)), [4, 1])
%!assert (size (ncx2rnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (ncx2rnd (1, 1, 0, 1)), [0, 1])
%!assert (size (ncx2rnd (1, 1, 1, 0)), [1, 0])
%!assert (size (ncx2rnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])
%!assert (size (ncx2rnd (1, 1, [])), [0, 0])
%!assert (size (ncx2rnd (1, 1, [2, 0, 2, 1])), [2, 0, 2])

## Test class of input preserved
%!assert (class (ncx2rnd (1, 1)), "double")
%!assert (class (ncx2rnd (1, single (1))), "single")
%!assert (class (ncx2rnd (1, single ([1, 1]))), "single")
%!assert (class (ncx2rnd (single (1), 1)), "single")
%!assert (class (ncx2rnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<ncx2rnd: function called with too few input arguments.> ncx2rnd ()
%!error<ncx2rnd: function called with too few input arguments.> ncx2rnd (1)
%!error<ncx2rnd: DF and LAMBDA must be of common size or scalars.> ...
%! ncx2rnd (ones (3), ones (2))
%!error<ncx2rnd: DF and LAMBDA must be of common size or scalars.> ...
%! ncx2rnd (ones (2), ones (3))
%!error<ncx2rnd: DF and LAMBDA must not be complex.> ncx2rnd (i, 2)
%!error<ncx2rnd: DF and LAMBDA must not be complex.> ncx2rnd (1, i)
%!error<ncx2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncx2rnd (1, 2, -1)
%!error<ncx2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncx2rnd (1, 2, 1.2)
%!error<ncx2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncx2rnd (1, 2, ones (2))
%!error<ncx2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncx2rnd (1, 2, [2 -1 2])
%!error<ncx2rnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncx2rnd (1, 2, [2 0 2.5])
%!error<ncx2rnd: dimensions must be non-negative integers.> ...
%! ncx2rnd (1, 2, 2, -1, 5)
%!error<ncx2rnd: dimensions must be non-negative integers.> ...
%! ncx2rnd (1, 2, 2, 1.5, 5)
%!error<ncx2rnd: DF and LAMBDA must be scalars or of size SZ.> ...
%! ncx2rnd (2, ones (2), 3)
%!error<ncx2rnd: DF and LAMBDA must be scalars or of size SZ.> ...
%! ncx2rnd (2, ones (2), [3, 2])
%!error<ncx2rnd: DF and LAMBDA must be scalars or of size SZ.> ...
%! ncx2rnd (2, ones (2), 3, 2)
