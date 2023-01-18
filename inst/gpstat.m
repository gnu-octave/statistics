## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
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
## @deftypefn  {statistics} [@var{m}, @var{v}] = gpstat (@var{shape}, @var{scale}, @var{location})
##
## Mean and variance of the generalized Pareto distribution.
##
## @code{[@var{m}, @var{v}] = gpstat (@var{shape}, @var{scale}, @var{location})}
## returns the mean and variance of the generalized Pareto distribution with
## @var{shape}, @var{scale}, and @var{location} parameeters.
##
## The default value for @var{location} is 0.
##
## When @var{shape} = 0 and @var{location} = 0, the generalized Pareto
## distribution is equivalent to the exponential distribution.  When
## @code{@var{shape} > 0} and @code{@var{location} = @var{scale} / @var{shape}},
## the generalized Pareto distribution is equivalent to the Pareto distribution.
## The mean of the generalized Pareto distribution is not finite when
## @code{@var{shape} >= 1}, and the variance is not finite when
## @code{@var{shape} >= 1/2}.  When @code{@var{shape} >= 0}, the generalized
## Pareto distribution has positive density for @code{@var{x} > @var{location}},
## or, when @code{@var{shape} < 0}, for
## @code{0 <= (@var{x} -  @var{location}) / @var{scale} <= -1 / @var{shape}}.
##
## @seealso{gpcdf, gpinv, gppdf, gprnd, gpfit, gplike}
## @end deftypefn

function [m, v] = gpstat (shape, scale, location)

  ## Check input arguments
  if (nargin < 2)
      error ("gpstat: too few input arguments.");
  endif
  if (nargin < 3 || isempty(location))
    location = 0;
  endif

  [err, shape, scalelocationmu] = common_size (shape, scale, location);
  if (err > 0)
    error ("gpstat: inputs must be of common size or scalars.");
  endif

  ## Return NaNs for out of range SCALE parameters.
  scale(scale <= 0) = NaN;

  ## Check for appropriate class
  if (isa (shape, "single") || isa (scale, "single") || ...
                               isa (location, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (shape) || iscomplex (scale) || iscomplex (location))
    error ("gpstat: SHAPE, SCALE, and LOCATION must not be complex.");
  endif

  ## Prepare output
  m = NaN (size (shape), is_class);
  v = NaN (size (shape), is_class);

  ## Compute cases for SHAPE == 0
  k = (abs (shape) < eps (is_class));
  m(k) = 1;
  v(k) = 1;

  ## Compute cases for SHAPE != 0
  k = ! k;
  ## SHAPE < 1
  kk = k & (shape < 1);
  m(kk) = 1 ./ (1 - shape(kk));
  ## SHAPE > 1
  m(shape >= 1) = Inf;
  ## SHAPE < 1/2
  % Find the shape~=0 cases and fill in the variance.
  kk = k & (shape < 1/2);
  v(kk) = 1 ./ ((1-shape(kk)).^2 .* (1-2.*shape(kk)));
  ## SHAPE > 1/2
  v(shape >= 1/2) = Inf;

  ## Compute mean and variance
  m = location + scale .* m;
  v = scale.^2 .* v;

endfunction

## Test results
%!shared x, y
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! y = [0, 0.5, 1, 2, Inf, Inf];
%!assert (gpstat (x, ones (1,6), zeros (1,6)), y, eps)

## Test class of input preserved
%!assert (gpstat (single (x), 1, 0), single (y), eps("single"))
%!assert (gpstat (x, single (1), 0), single (y), eps("single"))
%!assert (gpstat (x, single (1)), single (y), eps("single"))
%!assert (gpstat (x, 1, single (0)), single (y), eps("single"))
%!assert (gpstat (single ([x, NaN]), 1, 0), single ([y, NaN]), eps("single"))
%!assert (gpstat ([x, NaN], single (1), 0), single ([y, NaN]), eps("single"))
%!assert (gpstat ([x, NaN], single (1)), single ([y, NaN]), eps("single"))
%!assert (gpstat ([x, NaN], 1, single (0)), single ([y, NaN]), eps("single"))

## Test input validation
%!error<gpstat: too few input arguments.> gpstat ()
%!error<gpstat: too few input arguments.> gpstat (1)
%!error<gpstat: inputs must be of common size or scalars.> ...
%! gpstat (1, ones (2), ones (3))
%!error<gpstat: inputs must be of common size or scalars.> ...
%! gpstat (rand (3), ones (2), ones (3))
%!error<gpstat: SHAPE, SCALE, and LOCATION must not be complex.> gpstat (i, 2, 2)
%!error<gpstat: SHAPE, SCALE, and LOCATION must not be complex.> gpstat (2, i, 2)
%!error<gpstat: SHAPE, SCALE, and LOCATION must not be complex.> gpstat (2, 2, i)
