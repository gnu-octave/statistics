## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
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
## @deftypefn  {statistics} {@var{r} =} gprnd (@var{k}, @var{sigma}, @var{theta})
## @deftypefnx {statistics} {@var{r} =} gprnd (@var{k}, @var{sigma}, @var{theta}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} gprnd (@var{k}, @var{sigma}, @var{theta}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} gprnd (@var{k}, @var{sigma}, @var{theta}, [@var{sz}])
##
## Random arrays from the generalized Pareto distribution.
##
## @code{@var{r} = gprnd (@var{k}, @var{sigma}, @var{theta})} returns an array
## of random numbers chosen from the generalized Pareto distribution with shape
## parameter @var{k}, scale parameter @var{sigma}, and location parameter
## @var{theta}.  The size of @var{r} is the common size of @var{k}, @var{sigma},
## and @var{theta}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## When called with a single size argument, @code{gprnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## When @qcode{@var{k} = 0} and @qcode{@var{theta} = 0}, the Generalized Pareto
## is equivalent to the exponential distribution.  When @qcode{@var{k} > 0} and
## @code{@var{theta} = @var{k} / @var{k}} the Generalized Pareto is equivalent
## to the Pareto distribution.  The mean of the Generalized Pareto is not finite
## when @qcode{@var{k} >= 1} and the variance is not finite when
## @qcode{@var{k} >= 1/2}.  When @qcode{@var{k} >= 0}, the Generalized Pareto
## has positive density for @qcode{@var{x} > @var{theta}}, or, when
## @qcode{@var{theta} < 0}, for
## @qcode{0 <= (@var{x} - @var{theta}) / @var{sigma} <= -1 / @var{k}}.
##
## Further information about the generalized Pareto distribution can be found at
## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
##
## @seealso{gpcdf, gpinv, gppdf, gpfit, gplike, gpstat}
## @end deftypefn

function r = gprnd (k, sigma, theta, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("gprnd: function called with too few input arguments.");
  endif

  ## Check for common size of K, SIGMA, and THETA
  if (! isscalar (k) || ! isscalar (sigma) || ! isscalar (theta))
    [retval, k, sigma, theta] = common_size (k, sigma, theta);
    if (retval > 0)
      error ("gprnd: K, SIGMA, and THETA must be of common size or scalars.");
    endif
  endif

  ## Check for K, SIGMA, and THETA being reals
  if (iscomplex (k) || iscomplex (sigma) || iscomplex (theta))
    error ("gprnd: K, SIGMA, and THETA must not be complex.");
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
      error (strcat ("gprnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 4)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("gprnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  ## Use 'size (ones (sz))' to ignore any trailing singleton dimensions in SZ
  if (!isscalar (k) && ! isequal (size (k), size (ones (sz))))
    error ("gprnd: K, SIGMA, and THETA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (k, "single") || isa (sigma, "single") || isa (theta, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from generalized Pareto distribution
  r = rand (sz, cls);

  ## Find valid parameters
  vr = (isfinite (r)) & (theta > -Inf) & (theta < Inf) ...
                      & (sigma > 0) & (sigma < Inf) ...
                      & (-Inf < k) & (k < Inf);

  ## Force invalid parameters to NaN
  r(! vr) = NaN;

  if (isscalar (k))
    if (k == 0)
      r(vr) = theta - (sigma .* log (1 - r(vr)));
    else
      r(vr) = theta + ((sigma .* ((r(vr) .^ -k) - 1)) ./ k);
    endif
  else
    if (any (k == 0))
      r(vr) = theta(vr) - (sigma(vr) .* log (1 - r(vr)));
    endif
    if (any (k < 0 | k > 0))
      r(vr) = theta(vr) + ((sigma(vr) .* ((r(vr) .^ -k(vr)) - 1)) ./ k(vr));
    endif
  endif

endfunction

## Test output
%!assert (size (gprnd (0, 1, 0)), [1, 1])
%!assert (size (gprnd (0, 1, zeros (2,1))), [2, 1])
%!assert (size (gprnd (0, 1, zeros (2,2))), [2, 2])
%!assert (size (gprnd (0, ones (2,1), 0)), [2, 1])
%!assert (size (gprnd (0, ones (2,2), 0)), [2, 2])
%!assert (size (gprnd (zeros (2,1), 1, 0)), [2, 1])
%!assert (size (gprnd (zeros (2,2), 1, 0)), [2, 2])
%!assert (size (gprnd (0, 1, 0, 3)), [3, 3])
%!assert (size (gprnd (0, 1, 0, [4 1])), [4, 1])
%!assert (size (gprnd (0, 1, 0, 4, 1)), [4, 1])
%!assert (size (gprnd (1,1,0)), [1, 1])
%!assert (size (gprnd (1, 1, zeros (2,1))), [2, 1])
%!assert (size (gprnd (1, 1, zeros (2,2))), [2, 2])
%!assert (size (gprnd (1, ones (2,1), 0)), [2, 1])
%!assert (size (gprnd (1, ones (2,2), 0)), [2, 2])
%!assert (size (gprnd (ones (2,1), 1, 0)), [2, 1])
%!assert (size (gprnd (ones (2,2), 1, 0)), [2, 2])
%!assert (size (gprnd (1, 1, 0, 3)), [3, 3])
%!assert (size (gprnd (1, 1, 0, [4 1])), [4, 1])
%!assert (size (gprnd (1, 1, 0, 4, 1)), [4, 1])
%!assert (size (gprnd (-1, 1, 0)), [1, 1])
%!assert (size (gprnd (-1, 1, zeros (2,1))), [2, 1])
%!assert (size (gprnd (1, -1, zeros (2,2))), [2, 2])
%!assert (size (gprnd (-1, ones (2,1), 0)), [2, 1])
%!assert (size (gprnd (-1, ones (2,2), 0)), [2, 2])
%!assert (size (gprnd (-ones (2,1), 1, 0)), [2, 1])
%!assert (size (gprnd (-ones (2,2), 1, 0)), [2, 2])
%!assert (size (gprnd (-1, 1, 0, 3)), [3, 3])
%!assert (size (gprnd (-1, 1, 0, [4, 1])), [4, 1])
%!assert (size (gprnd (-1, 1, 0, 4, 1)), [4, 1])
%!assert (size (gprnd (-1, 1, 0, [])), [0, 0])
%!assert (size (gprnd (-1, 1, 0, [2, 0, 2, 1])), [2, 0, 2])

## Test class of input preserved
%!assert (class (gprnd (0, 1, 0)), "double")
%!assert (class (gprnd (0, 1, single (0))), "single")
%!assert (class (gprnd (0, 1, single ([0, 0]))), "single")
%!assert (class (gprnd (0, single (1),0)), "single")
%!assert (class (gprnd (0, single ([1, 1]),0)), "single")
%!assert (class (gprnd (single (0), 1, 0)), "single")
%!assert (class (gprnd (single ([0, 0]), 1, 0)), "single")

## Test input validation
%!error<gprnd: function called with too few input arguments.> gprnd ()
%!error<gprnd: function called with too few input arguments.> gprnd (1)
%!error<gprnd: function called with too few input arguments.> gprnd (1, 2)
%!error<gprnd: K, SIGMA, and THETA must be of common size or scalars.> ...
%! gprnd (ones (3), ones (2), ones (2))
%!error<gprnd: K, SIGMA, and THETA must be of common size or scalars.> ...
%! gprnd (ones (2), ones (3), ones (2))
%!error<gprnd: K, SIGMA, and THETA must be of common size or scalars.> ...
%! gprnd (ones (2), ones (2), ones (3))
%!error<gprnd: K, SIGMA, and THETA must not be complex.> gprnd (i, 2, 3)
%!error<gprnd: K, SIGMA, and THETA must not be complex.> gprnd (1, i, 3)
%!error<gprnd: K, SIGMA, and THETA must not be complex.> gprnd (1, 2, i)
%!error<gprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gprnd (1, 2, 3, -1)
%!error<gprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gprnd (1, 2, 3, 1.2)
%!error<gprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gprnd (1, 2, 3, ones (2))
%!error<gprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gprnd (1, 2, 3, [2 -1 2])
%!error<gprnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! gprnd (1, 2, 3, [2 0 2.5])
%!error<gprnd: dimensions must be non-negative integers.> ...
%! gprnd (1, 2, 3, 2, -1, 5)
%!error<gprnd: dimensions must be non-negative integers.> ...
%! gprnd (1, 2, 3, 2, 1.5, 5)
%!error<gprnd: K, SIGMA, and THETA must be scalars or of size SZ.> ...
%! gprnd (2, ones (2), 2, 3)
%!error<gprnd: K, SIGMA, and THETA must be scalars or of size SZ.> ...
%! gprnd (2, ones (2), 2, [3, 2])
%!error<gprnd: K, SIGMA, and THETA must be scalars or of size SZ.> ...
%! gprnd (2, ones (2), 2, 3, 2)
