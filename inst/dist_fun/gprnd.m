## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} gprnd (@var{k}, @var{sigma}, @var{mu})
## @deftypefnx {statistics} {@var{r} =} gprnd (@var{k}, @var{sigma}, @var{mu}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} gprnd (@var{k}, @var{sigma}, @var{mu}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} gprnd (@var{k}, @var{sigma}, @var{mu}, [@var{sz}])
##
## Random arrays from the generalized Pareto distribution.
##
## @code{@var{r} = gprnd (@var{k}, @var{sigma}, @var{mu})} returns an array of
## random numbers chosen from the generalized Pareto distribution with shape
## parameter @var{k}, scale parameter @var{sigma}, and location parameter
## @var{mu}.  The size of @var{r} is the common size of the input arguments.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## When @qcode{@var{k} = 0} and @qcode{@var{mu} = 0}, the Generalized Pareto CDF
## is equivalent to the exponential distribution.  When @qcode{@var{k} > 0} and
## @code{@var{mu} = @var{k} / @var{k}} the Generalized Pareto is equivalent to
## the Pareto distribution.  The mean of the Generalized Pareto is not finite
## when @qcode{@var{k} >= 1} and the variance is not finite when
## @qcode{@var{k} >= 1/2}.  When @qcode{@var{k} >= 0}, the Generalized Pareto
## has positive density for @qcode{@var{x} > @var{mu}}, or, when
## @qcode{@var{mu} < 0},for
## @qcode{0 <= (@var{x} - @var{mu}) / @var{sigma} <= -1 / @var{k}}.
##
## Further information about the generalized Pareto distribution can be found at
## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
##
## @seealso{gpcdf, gpinv, gppdf, gpfit, gplike, gpstat}
## @end deftypefn

function r = gprnd (k, sigma, mu, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("gprnd: function called with too few input arguments.");
  endif

  ## Check for common size of K, SIGMA, and MU
  if (! isscalar (k) || ! isscalar (sigma) || ! isscalar (mu))
    [retval, k, sigma, mu] = common_size (k, sigma, mu);
    if (retval > 0)
      error ("gprnd: K, SIGMA, and MU must be of common size or scalars.");
    endif
  endif

  ## Check for K, SIGMA, and MU being reals
  if (iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gprnd: K, SIGMA, and MU must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 3)
    sz = size (k);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["gprnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 4)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("gprnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (!isscalar (k) && ! isequal (size (k), sz))
    error ("gprnd: K, SIGMA, and MU must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (k, "single") || isa (sigma, "single") || isa (mu, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from generalized Pareto distribution
  r = NaN (sz, cls);

  kr = (-Inf < mu) & (mu < Inf) & (sigma > 0) & (sigma < Inf) ...
                                & (-Inf < k) & (k < Inf);
  r(kr(:)) = rand (1, sum(kr(:)), cls);
  if (any (k == 0))
      r(kr) = mu(kr) - (sigma(kr) .* log(1 - r(kr)));
  elseif (any (k < 0 | k > 0))
    r(kr) = mu(kr) + ((sigma(kr) .* ((r(kr) .^ -k(kr)) - 1)) ./ k(kr));
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
%!error<gprnd: K, SIGMA, and MU must be of common size or scalars.> ...
%! gprnd (ones (3), ones (2), ones (2))
%!error<gprnd: K, SIGMA, and MU must be of common size or scalars.> ...
%! gprnd (ones (2), ones (3), ones (2))
%!error<gprnd: K, SIGMA, and MU must be of common size or scalars.> ...
%! gprnd (ones (2), ones (2), ones (3))
%!error<gprnd: K, SIGMA, and MU must not be complex.> gprnd (i, 2, 3)
%!error<gprnd: K, SIGMA, and MU must not be complex.> gprnd (1, i, 3)
%!error<gprnd: K, SIGMA, and MU must not be complex.> gprnd (1, 2, i)
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
%!error<gprnd: K, SIGMA, and MU must be scalar or of size SZ.> ...
%! gprnd (2, ones (2), 2, 3)
%!error<gprnd: K, SIGMA, and MU must be scalar or of size SZ.> ...
%! gprnd (2, ones (2), 2, [3, 2])
%!error<gprnd: K, SIGMA, and MU must be scalar or of size SZ.> ...
%! gprnd (2, ones (2), 2, 3, 2)
