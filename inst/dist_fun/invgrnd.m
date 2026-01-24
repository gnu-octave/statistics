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
## @deftypefn  {statistics} {@var{r} =} invgrnd (@var{mu}, @var{lambda})
## @deftypefnx {statistics} {@var{r} =} invgrnd (@var{mu}, @var{lambda}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} invgrnd (@var{mu}, @var{lambda}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} invgrnd (@var{mu}, @var{lambda}, [@var{sz}])
##
## Random arrays from the inverse Gaussian distribution.
##
## @code{@var{r} = invgrnd (@var{mu}, @var{lambda})} returns an array of random
## numbers chosen from the inverse Gaussian distribution with location parameter
## @var{mu} and scale parameter @var{lambda}.  The size of @var{r} is the common
## size of @var{mu} and @var{lambda}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## When called with a single size argument, @code{invgrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## The inverse Gaussian CDF is only defined for @qcode{@var{mu} > 0} and
## @qcode{@var{lambda} > 0}.
##
## Further information about the inverse Gaussian distribution can be found at
## @url{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}
##
## @seealso{invgcdf, invginv, invgpdf, invgfit, invglike, invgstat}
## @end deftypefn

function r = invgrnd (mu, lambda, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("invgrnd: function called with too few input arguments.");
  endif

  ## Check for common size of MU, and LAMBDA
  if (! isscalar (mu) || ! isscalar (lambda))
    [retval, mu, lambda] = common_size (mu, lambda);
    vec = true;
    if (retval > 0)
      error ("invgrnd: MU and LAMBDA must be of common size or scalars.");
    endif
  else
    vec = false;
  endif

  ## Check for X, MU, and LAMBDA being reals
  if (iscomplex (mu) || iscomplex (lambda))
    error ("invgrnd: MU and LAMBDA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (mu);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif ((isrow (varargin{1}) || isempty (varargin{1})) &&
            all (varargin{1} >= 0) && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("invgrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("invgrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("invgrnd: MU and LAMBDA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single") || isa (lambda, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Expand parameters (if needed)
  if (! vec)
    mu = repmat (mu, sz);
    lambda = repmat (lambda, sz);
  endif

  ## Generate random sample from inverse Gaussian distribution
  v = randn (sz, cls);
  y = v .^ 2;
  r = mu + (mu .^ 2 .* y) ./ (2 .* lambda) - (mu ./ (2 .* lambda)) .* ...
           sqrt (4 * mu .* lambda .* y + mu .* mu .* y .* y);

  inver = (rand (sz) .* (mu + r) > mu);
  r(inver) = mu(inver) .^2 ./ r(inver);

  ## Force output to NaN for invalid parameters MU and LAMBDA
  k = (mu <= 0 | lambda <= 0);
  r(k) = NaN;

endfunction

## Test results
%!assert (size (invgrnd (1, 1, 1)), [1, 1])
%!assert (size (invgrnd (1, 1, 2)), [2, 2])
%!assert (size (invgrnd (1, 1, [2, 1])), [2, 1])
%!assert (size (invgrnd (1, zeros (2, 2))), [2, 2])
%!assert (size (invgrnd (1, ones (2, 1))), [2, 1])
%!assert (size (invgrnd (1, ones (2, 2))), [2, 2])
%!assert (size (invgrnd (ones (2, 1), 1)), [2, 1])
%!assert (size (invgrnd (ones (2, 2), 1)), [2, 2])
%!assert (size (invgrnd (1, 1, 3)), [3, 3])
%!assert (size (invgrnd (1, 1, [4 1])), [4, 1])
%!assert (size (invgrnd (1, 1, 4, 1)), [4, 1])
%!assert (size (invgrnd (1, 1, [])), [0, 0])
%!assert (size (invgrnd (1, 1, [2, 0, 2, 1])), [2, 0, 2])
%!test
%! r =  invgrnd (1, [1, 0, -1]);
%! assert (r([2:3]), [NaN, NaN])

## Test class of input preserved
%!assert (class (invgrnd (1, 0)), "double")
%!assert (class (invgrnd (1, single (0))), "single")
%!assert (class (invgrnd (1, single ([0, 0]))), "single")
%!assert (class (invgrnd (1, single (1))), "single")
%!assert (class (invgrnd (1, single ([1, 1]))), "single")
%!assert (class (invgrnd (single (1), 1)), "single")
%!assert (class (invgrnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<invgrnd: function called with too few input arguments.> invgrnd ()
%!error<invgrnd: function called with too few input arguments.> invgrnd (1)
%!error<invgrnd: MU and LAMBDA must be of common size or scalars.> ...
%! invgrnd (ones (3), ones (2))
%!error<invgrnd: MU and LAMBDA must be of common size or scalars.> ...
%! invgrnd (ones (2), ones (3))
%!error<invgrnd: MU and LAMBDA must not be complex.> invgrnd (i, 2, 3)
%!error<invgrnd: MU and LAMBDA must not be complex.> invgrnd (1, i, 3)
%!error<invgrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! invgrnd (1, 2, -1)
%!error<invgrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! invgrnd (1, 2, 1.2)
%!error<invgrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! invgrnd (1, 2, ones (2))
%!error<invgrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! invgrnd (1, 2, [2 -1 2])
%!error<invgrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! invgrnd (1, 2, [2 0 2.5])
%!error<invgrnd: dimensions must be non-negative integers.> ...
%! invgrnd (1, 2, 2, -1, 5)
%!error<invgrnd: dimensions must be non-negative integers.> ...
%! invgrnd (1, 2, 2, 1.5, 5)
%!error<invgrnd: MU and LAMBDA must be scalars or of size SZ.> ...
%! invgrnd (2, ones (2), 3)
%!error<invgrnd: MU and LAMBDA must be scalars or of size SZ.> ...
%! invgrnd (2, ones (2), [3, 2])
%!error<invgrnd: MU and LAMBDA must be scalars or of size SZ.> ...
%! invgrnd (2, ones (2), 3, 2)
