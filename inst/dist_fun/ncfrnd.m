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
## @deftypefn  {statistics} {@var{r} =} ncfrnd (@var{df1}, @var{df2}, @var{lambda})
## @deftypefnx {statistics} {@var{r} =} ncfrnd (@var{df1}, @var{df2}, @var{lambda}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} ncfrnd (@var{df1}, @var{df2}, @var{lambda}, [@var{sz}])
##
## Random arrays from the noncentral @math{F}-distribution.
##
## @code{@var{x} = ncfrnd (@var{p}, @var{df1}, @var{df2}, @var{lambda})} returns
## an array of random numbers chosen from the noncentral @math{F}-distribution with
## @var{df1} and @var{df2} degrees of freedom and noncentrality parameter
## @var{lambda}.  The size of @var{r} is the common size of @var{df1},
## @var{df2}, and @var{lambda}.  A scalar input functions as a constant matrix
## of the same size as the other input.
##
## @code{ncfrnd} generates values using the definition of a noncentral @math{F}
## random variable, as the ratio of a noncentral chi-squared distribution and a
## (central) chi-squared distribution.
##
## When called with a single size argument, @code{ncfrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the noncentral @math{F}-distribution can be found
## at @url{https://en.wikipedia.org/wiki/Noncentral_F-distribution}
##
## @seealso{ncfcdf, ncfinv, ncfpdf, ncfstat, frnd, ncx2rnd, chi2rnd}
## @end deftypefn

function r = ncfrnd (df1, df2, lambda, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("ncfrnd: function called with too few input arguments.");
  endif

  ## Check for common size of DF1, DF2, and LAMBDA
  if (! isscalar (df1) || ! isscalar (df2) || ! isscalar (lambda))
    [retval, df1, df2, lambda] = common_size (df1, df2, lambda);
    if (retval > 0)
      error ("ncfrnd: DF1, DF2, and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Check for DF1, DF2, and LAMBDA being reals
  if (iscomplex (df1) || iscomplex (df2) || iscomplex (lambda))
    error ("ncfrnd: DF1, DF2, and LAMBDA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 3)
    sz = size (df1);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif ((isrow (varargin{1}) || isempty (varargin{1})) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("ncfrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 4)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("ncfrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (df1) && ! isequal (size (df1), sz))
    error ("ncfrnd: DF1, DF2, and LAMBDA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (df1, "single") || isa (df2, "single") || isa (lambda, "single"));
    cls = "single";
  else
    cls = "double";
  endif

  ## Return NaNs for out of range values of DF1, DF2, and LAMBDA
  df1(df1 <= 0) = NaN;
  df2(df2 <= 0) = NaN;
  lambda(lambda <= 0) = NaN;

  ## Generate random sample from noncentral F distribution
  r = (ncx2rnd (df1, lambda, sz) ./ df1) ./ ...
      (2 .* randg (df2 ./ 2, sz) ./ df2);

  ## Cast to appropriate class
  r = cast (r, cls);

endfunction

## Test output
%!assert (size (ncfrnd (1, 1, 1)), [1 1])
%!assert (size (ncfrnd (1, ones (2,1), 1)), [2, 1])
%!assert (size (ncfrnd (1, ones (2,2), 1)), [2, 2])
%!assert (size (ncfrnd (ones (2,1), 1, 1)), [2, 1])
%!assert (size (ncfrnd (ones (2,2), 1, 1)), [2, 2])
%!assert (size (ncfrnd (1, 1, 1, 3)), [3, 3])
%!assert (size (ncfrnd (1, 1, 1, [4, 1])), [4, 1])
%!assert (size (ncfrnd (1, 1, 1, 4, 1)), [4, 1])
%!assert (size (ncfrnd (1, 1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (ncfrnd (1, 1, 1, 0, 1)), [0, 1])
%!assert (size (ncfrnd (1, 1, 1, 1, 0)), [1, 0])
%!assert (size (ncfrnd (1, 1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (ncfrnd (1, 1, 1)), "double")
%!assert (class (ncfrnd (1, single (1), 1)), "single")
%!assert (class (ncfrnd (1, 1, single (1))), "single")
%!assert (class (ncfrnd (1, single ([1, 1]), 1)), "single")
%!assert (class (ncfrnd (1, 1, single ([1, 1]))), "single")
%!assert (class (ncfrnd (single (1), 1, 1)), "single")
%!assert (class (ncfrnd (single ([1, 1]), 1, 1)), "single")

## Test input validation
%!error<ncfrnd: function called with too few input arguments.> ncfrnd ()
%!error<ncfrnd: function called with too few input arguments.> ncfrnd (1)
%!error<ncfrnd: function called with too few input arguments.> ncfrnd (1, 2)
%!error<ncfrnd: DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfrnd (ones (3), ones (2), ones (2))
%!error<ncfrnd: DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfrnd (ones (2), ones (3), ones (2))
%!error<ncfrnd: DF1, DF2, and LAMBDA must be of common size or scalars.> ...
%! ncfrnd (ones (2), ones (2), ones (3))
%!error<ncfrnd: DF1, DF2, and LAMBDA must not be complex.> ncfrnd (i, 2, 3)
%!error<ncfrnd: DF1, DF2, and LAMBDA must not be complex.> ncfrnd (1, i, 3)
%!error<ncfrnd: DF1, DF2, and LAMBDA must not be complex.> ncfrnd (1, 2, i)
%!error<ncfrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncfrnd (1, 2, 3, -1)
%!error<ncfrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncfrnd (1, 2, 3, 1.2)
%!error<ncfrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncfrnd (1, 2, 3, ones (2))
%!error<ncfrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncfrnd (1, 2, 3, [2 -1 2])
%!error<ncfrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! ncfrnd (1, 2, 3, [2 0 2.5])
%!error<ncfrnd: dimensions must be non-negative integers.> ...
%! ncfrnd (1, 2, 3, 2, -1, 5)
%!error<ncfrnd: dimensions must be non-negative integers.> ...
%! ncfrnd (1, 2, 3, 2, 1.5, 5)
%!error<ncfrnd: DF1, DF2, and LAMBDA must be scalars or of size SZ.> ...
%! ncfrnd (2, ones (2), 2, 3)
%!error<ncfrnd: DF1, DF2, and LAMBDA must be scalars or of size SZ.> ...
%! ncfrnd (2, ones (2), 2, [3, 2])
%!error<ncfrnd: DF1, DF2, and LAMBDA must be scalars or of size SZ.> ...
%! ncfrnd (2, ones (2), 2, 3, 2)
