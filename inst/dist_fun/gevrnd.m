## Copyright (C) 2012 Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{r} =} gevrnd (@var{k}, @var{sigma}, @var{mu})
## @deftypefnx {statistics} {@var{r} =} gevrnd (@var{k}, @var{sigma}, @var{mu}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} gevrnd (@var{k}, @var{sigma}, @var{mu}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} gevrnd (@var{k}, @var{sigma}, @var{mu}, [@var{sz}])
##
## Random arrays from the generalized extreme value (GEV) distribution.
##
## @code{@var{r} = gevrnd (@var{k}, @var{sigma}, @var{mu}} returns an array of
## random numbers chosen from the GEV distribution with shape parameter @var{k},
## scale parameter @var{sigma}, and location parameter @var{mu}.  The size of
## @var{r} is the common size of @var{k}, @var{sigma}, and @var{mu}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## When called with a single size argument, @code{gevrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## When @qcode{@var{k} < 0}, the GEV is the type III extreme value distribution.
## When @qcode{@var{k} > 0}, the GEV distribution is the type II, or Frechet,
## extreme value distribution.  If @var{W} has a Weibull distribution as
## computed by the @code{wblcdf} function, then @qcode{-@var{W}} has a type III
## extreme value distribution and @qcode{1/@var{W}} has a type II extreme value
## distribution.  In the limit as @var{k} approaches @qcode{0}, the GEV is the
## mirror image of the type I extreme value distribution as computed by the
## @code{evcdf} function.
##
## The mean of the GEV distribution is not finite when @qcode{@var{k} >= 1}, and
## the variance is not finite when @qcode{@var{k} >= 1/2}.  The GEV distribution
## has positive density only for values of @var{x} such that
## @qcode{@var{k} * (@var{x} - @var{mu}) / @var{sigma} > -1}.
##
## Further information about the generalized extreme value distribution can be
## found at
## @url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
##
## @seealso{gevcdf, gevinv, gevpdf, gevfit, gevlike, gevstat}
## @end deftypefn

function r = gevrnd (k, sigma, mu, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("gevrnd: function called with too few input arguments.");
  endif

  ## Check for common size of K, SIGMA, and MU
  if (! isscalar (k) || ! isscalar (sigma) || ! isscalar (mu))
    [retval, k, sigma, mu] = common_size (k, sigma, mu);
    if (retval > 0)
      error ("gevrnd: K, SIGMA, and MU must be of common size or scalars.");
    endif
  endif

  ## Check for K, SIGMA, and MU being reals
  if (iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gevrnd: K, SIGMA, and MU must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 3)
    sz = size (k);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0)
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif (isempty (varargin{1}))
      r = [];
      return;
    else
      error (strcat ("gevrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 4)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("gevrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  ## Use 'size (ones (sz))' to ignore any trailing singleton dimensions in SZ
  if (!isscalar (k) && ! isequal (size (k), size (ones (sz))))
    error ("gevrnd: K, SIGMA, and MU must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (k, "single") || isa (sigma, "single") || isa (mu, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Burr type XII distribution
  r = gevinv (rand(sz), k, sigma, mu);

  r = cast (r, cls);

endfunction

## Test output
%!assert (size (gevrnd (1, 2, 1)), [1, 1]);
%!assert (size (gevrnd (ones (2, 1), 2, 1)), [2, 1]);
%!assert (size (gevrnd (ones (2, 2), 2, 1)), [2, 2]);
%!assert (size (gevrnd (1, 2 * ones (2, 1), 1)), [2, 1]);
%!assert (size (gevrnd (1, 2 * ones (2, 2), 1)), [2, 2]);
%!assert (size (gevrnd (1, 2, 1, 3)), [3, 3]);
%!assert (size (gevrnd (1, 2, 1, [4, 1])), [4, 1]);
%!assert (size (gevrnd (1, 2, 1, 4, 1)), [4, 1]);
%!assert (size (gevrnd (1, 2, 1, [])), [0, 0])
%!assert (size (gevrnd (1, 2, 1, [2, 0, 2, 1])), [2, 0, 2])

## Test class of input preserved
%!assert (class (gevrnd (1,1,1)), "double")
%!assert (class (gevrnd (single (1),1,1)), "single")
%!assert (class (gevrnd (single ([1 1]),1,1)), "single")
%!assert (class (gevrnd (1,single (1),1)), "single")
%!assert (class (gevrnd (1,single ([1 1]),1)), "single")
%!assert (class (gevrnd (1,1,single (1))), "single")
%!assert (class (gevrnd (1,1,single ([1 1]))), "single")

## Test input validation
%!error<gevrnd: function called with too few input arguments.> gevrnd ()
%!error<gevrnd: function called with too few input arguments.> gevrnd (1)
%!error<gevrnd: function called with too few input arguments.> gevrnd (1, 2)
%!error<gevrnd: K, SIGMA, and MU must be of common size or scalars.> ...
%! gevrnd (ones (3), ones (2), ones (2))
%!error<gevrnd: K, SIGMA, and MU must be of common size or scalars.> ...
%! gevrnd (ones (2), ones (3), ones (2))
%!error<gevrnd: K, SIGMA, and MU must be of common size or scalars.> ...
%! gevrnd (ones (2), ones (2), ones (3))
%!error<gevrnd: K, SIGMA, and MU must not be complex.> gevrnd (i, 2, 3)
%!error<gevrnd: K, SIGMA, and MU must not be complex.> gevrnd (1, i, 3)
%!error<gevrnd: K, SIGMA, and MU must not be complex.> gevrnd (1, 2, i)
%!error<gevrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gevrnd (1, 2, 3, -1)
%!error<gevrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gevrnd (1, 2, 3, 1.2)
%!error<gevrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gevrnd (1, 2, 3, ones (2))
%!error<gevrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gevrnd (1, 2, 3, [2 -1 2])
%!error<gevrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gevrnd (1, 2, 3, [2 0 2.5])
%!error<gevrnd: dimensions must be non-negative integers.> ...
%! gevrnd (1, 2, 3, 2, -1, 5)
%!error<gevrnd: dimensions must be non-negative integers.> ...
%! gevrnd (1, 2, 3, 2, 1.5, 5)
%!error<gevrnd: K, SIGMA, and MU must be scalars or of size SZ.> ...
%! gevrnd (2, ones (2), 2, 3)
%!error<gevrnd: K, SIGMA, and MU must be scalars or of size SZ.> ...
%! gevrnd (2, ones (2), 2, [3, 2])
%!error<gevrnd: K, SIGMA, and MU must be scalars or of size SZ.> ...
%! gevrnd (2, ones (2), 2, 3, 2)
