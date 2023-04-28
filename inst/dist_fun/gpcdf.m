## Copyright (C) 1997-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2018 John Donoghue
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} gpcdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} gpcdf (@var{x}, @var{shape})
## @deftypefnx {statistics} {@var{p} =} gpcdf (@var{x}, @var{shape}, @var{scale})
## @deftypefnx {statistics} {@var{p} =} gpcdf (@var{x}, @var{shape}, @var{scale}, @var{location})
## @deftypefnx {statistics} {@var{p} =} gpcdf (@dots{}, "upper")
##
## Generalized Pareto cumulative distribution function (cdf).
##
## Compute the cumulative distribution function (CDF) at @var{x} of the
## generalized Pareto distribution with parameters @var{location}, @var{scale},
## and @var{shape}.  The size of P is the common size of the input arguments.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Default values for @var{shape}, @var{scale}, and @var{location} are 0, 1,
## and 0, respectively.
##
## When @code{@var{shape} = 0} and @code{@var{location} = 0}, the Generalized
## Pareto CDF is equivalent to the exponential distribution.
## When @code{@var{shape} > 0} and @code{@var{location} = @var{shape} /
## @var{shape}} the Generalized Pareto is equivalent to the Pareto distribution.
## The mean of the Generalized Pareto is not finite when @code{@var{shape} >= 1}
## and the variance is not finite when @code{@var{shape} >= 1/2}.  When
## @code{@var{shape} >= 0}, the Generalized Pareto has positive density for
## @code{@var{x} > @var{location}}, or, when @code{@var{location} < 0},for
## @code{0 <= (@var{x} - @var{location}) / @var{scale} <= -1 / @var{shape}}.
##
## @code{[@dots{}] = gpcdf(@dots{}, "upper")} returns the upper tail probability
## of the generalized Pareto distribution.
##
## @seealso{gpinv, gppdf, gprnd, gpfit, gplike, gpstat}
## @end deftypefn

function p = gpcdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 5)
    error ("gpcdf: invalid number of input arguments.");
  endif

  ## Check for 'upper' flag
  if (nargin > 1 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 1  && ischar (varargin{end}) && ...
          ! strcmpi (varargin{end}, "upper"))
    error ("gpcdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Get extra arguments (if they exist) or add defaults
  if (numel (varargin) > 0)
    shape = varargin{1};
  else
    shape = 0;
  endif
  if (numel (varargin) > 1)
    scale = varargin{2};
  else
    scale = 1;
  endif
  if (numel (varargin) > 2)
    location = varargin{3};
  else
    location = 0;
  endif

  ## Check for common size of X, SHAPE, SCALE, and LOCATION
  if (! isscalar (x) || ! isscalar (shape) || ! isscalar (scale) || ...
      ! isscalar (location))
    [err, x, shape, scale, location] = common_size (x, shape, scale, location);
    if (err > 0)
      error (strcat (["gpcdf: X, SHAPE, SCALE, and LOCATION"], ...
                     [" must be of common size or scalars."]));
    endif
  endif

  ## Check for X, SHAPE, SCALE, and LOCATION being reals
  if (iscomplex (x) || iscomplex (shape) || iscomplex (scale) || ...
      iscomplex (location))
    error ("gpcdf: X, SHAPE, SCALE, and LOCATION must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (shape, "single") || ...
      isa (scale, "single") || isa (location, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  p = zeros (size (x), is_class);

  ## Return NaNs for out of range values of scale parameter
  scale(scale <= 0) = NaN;

  ## Calculate (x-location)/scale => 0 and force zero below that
  z = (x - location) ./ scale;
  z(z < 0) = 0;

  ## Compute cases for SHAPE == 0
  kz = (abs (shape) < eps (is_class));
  if (uflag)
    p(kz) = exp (-z(kz));
  else
    p(kz) = -expm1 (-z(kz));
  endif

  ## For SHAPE < 0, calculate 0 <= x/sigma <= -1/k and force zero below that
  t = z .* shape;
  kt = (t <= -1 & shape < -eps (is_class));
  t(kt) = 0;

  ## Compute cases for SHAPE != 0
  kz = ! kz;
  if (uflag)
    p(kz) = exp ((-1 ./ shape(kz)) .* log1p (t(kz)));
  else
    p(kz) = -expm1 ((-1 ./ shape(kz)) .* log1p (t(kz)));
  endif
  if (uflag)
    p(kt) = 0;
  else
    p(kt) = 1;
  endif

  ## For SHAPE == NaN force p = NaN
  p(isnan (shape)) = NaN;

endfunction

## Test results
%!shared x, y1, y1u, y2, y2u, y3, y3u
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! y1 = [0, 0, 0, 0.3934693402873666, 0.6321205588285577, 1];
%! y1u = [1, 1, 1, 0.6065306597126334, 0.3678794411714423, 0];
%! y2 = [0, 0, 0, 1/3, 1/2, 1];
%! y2u = [1, 1, 1, 2/3, 1/2, 0];
%! y3 = [0, 0, 0, 1/2, 1, 1];
%! y3u = [1, 1, 1, 1/2, 0, 0];
%!assert (gpcdf (x, zeros (1,6), ones (1,6), zeros (1,6)), y1, eps)
%!assert (gpcdf (x, 0, 1, zeros (1,6)), y1, eps)
%!assert (gpcdf (x, 0, ones (1,6), 0), y1, eps)
%!assert (gpcdf (x, zeros (1,6), 1, 0), y1, eps)
%!assert (gpcdf (x, 0, 1, 0), y1, eps)
%!assert (gpcdf (x, 0, 1, [0, 0, 0, NaN, 0, 0]), [y1(1:3), NaN, y1(5:6)], eps)
%!assert (gpcdf (x, 0, [1, 1, 1, NaN, 1, 1], 0), [y1(1:3), NaN, y1(5:6)], eps)
%!assert (gpcdf (x, [0, 0, 0, NaN, 0, 0], 1, 0), [y1(1:3), NaN, y1(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], 0, 1, 0), [y1(1:3), NaN, y1(5:6)], eps)
%!assert (gpcdf (x, zeros (1,6), ones (1,6), zeros (1,6), "upper"), y1u, eps)
%!assert (gpcdf (x, 0, 1, zeros (1,6), "upper"), y1u, eps)
%!assert (gpcdf (x, 0, ones (1,6), 0, "upper"), y1u, eps)
%!assert (gpcdf (x, zeros (1,6), 1, 0, "upper"), y1u, eps)
%!assert (gpcdf (x, 0, 1, 0, "upper"), y1u, eps)

%!assert (gpcdf (x, ones (1,6), ones (1,6), zeros (1,6)), y2, eps)
%!assert (gpcdf (x, 1, 1, zeros (1,6)), y2, eps)
%!assert (gpcdf (x, 1, ones (1,6), 0), y2, eps)
%!assert (gpcdf (x, ones (1,6), 1, 0), y2, eps)
%!assert (gpcdf (x, 1, 1, 0), y2, eps)
%!assert (gpcdf (x, 1, 1, [0, 0, 0, NaN, 0, 0]), [y2(1:3), NaN, y2(5:6)], eps)
%!assert (gpcdf (x, 1, [1, 1, 1, NaN, 1, 1], 0), [y2(1:3), NaN, y2(5:6)], eps)
%!assert (gpcdf (x, [1, 1, 1, NaN, 1, 1], 1, 0), [y2(1:3), NaN, y2(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], 1, 1, 0), [y2(1:3), NaN, y2(5:6)], eps)
%!assert (gpcdf (x, ones (1,6), ones (1,6), zeros (1,6), "upper"), y2u, eps)
%!assert (gpcdf (x, 1, 1, zeros (1,6), "upper"), y2u, eps)
%!assert (gpcdf (x, 1, ones (1,6), 0, "upper"), y2u, eps)
%!assert (gpcdf (x, ones (1,6), 1, 0, "upper"), y2u, eps)
%!assert (gpcdf (x, 1, 1, 0, "upper"), y2u, eps)
%!assert (gpcdf (x, 1, 1, [0, 0, 0, NaN, 0, 0], "upper"), ...
%!                        [y2u(1:3), NaN, y2u(5:6)], eps)
%!assert (gpcdf (x, 1, [1, 1, 1, NaN, 1, 1], 0, "upper"), ...
%!                        [y2u(1:3), NaN, y2u(5:6)], eps)
%!assert (gpcdf (x, [1, 1, 1, NaN, 1, 1], 1, 0, "upper"), ...
%!                        [y2u(1:3), NaN, y2u(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], 1, 1, 0, "upper"), ...
%!                        [y2u(1:3), NaN, y2u(5:6)], eps)
%!assert (gpcdf (x, -ones (1,6), ones (1,6), zeros (1,6)), y3, eps)
%!assert (gpcdf (x, -1, 1, zeros (1,6)), y3, eps)
%!assert (gpcdf (x, -1, ones (1,6), 0), y3, eps)
%!assert (gpcdf (x, -ones (1,6), 1, 0), y3, eps)
%!assert (gpcdf (x, -1, 1, 0), y3, eps)
%!assert (gpcdf (x, -1, 1, [0, 0, 0, NaN, 0, 0]), [y3(1:3), NaN, y3(5:6)], eps)
%!assert (gpcdf (x, -1, [1, 1, 1, NaN, 1, 1], 0), [y3(1:3), NaN, y3(5:6)], eps)
%!assert (gpcdf (x, [-1, -1, -1, NaN, -1, -1], 1, 0), [y3(1:3), NaN, y3(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], -1, 1, 0), [y3(1:3), NaN, y3(5:6)], eps)
%!assert (gpcdf (x, -ones (1,6), ones (1,6), zeros (1,6), "upper"), y3u, eps)
%!assert (gpcdf (x, -1, 1, zeros (1,6), "upper"), y3u, eps)
%!assert (gpcdf (x, -1, ones (1,6), 0, "upper"), y3u, eps)
%!assert (gpcdf (x, -ones (1,6), 1, 0, "upper"), y3u, eps)
%!assert (gpcdf (x, -1, 1, 0, "upper"), y3u, eps)
%!assert (gpcdf (x, -1, 1, [0, 0, 0, NaN, 0, 0], "upper"), ...
%!                         [y3u(1:3), NaN, y3u(5:6)], eps)
%!assert (gpcdf (x, -1, [1, 1, 1, NaN, 1, 1], 0, "upper"), ...
%!                          [y3u(1:3), NaN, y3u(5:6)], eps)
%!assert (gpcdf (x, [-1, -1, -1, NaN, -1, -1], 1, 0, "upper"), ...
%!                          [y3u(1:3), NaN, y3u(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], -1, 1, 0, "upper"), ...
%!                          [y3u(1:3), NaN, y3u(5:6)], eps)

## Test class of input preserved
%!assert (gpcdf (single ([x, NaN]), 0, 1, 0), single ([y1, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], 0, 1, single (0)), single ([y1, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], 0, single (1), 0), single ([y1, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], single (0), 1, 0), single ([y1, NaN]), eps("single"))

%!assert (gpcdf (single ([x, NaN]), 1, 1, 0), single ([y2, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], 1, 1, single (0)), single ([y2, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], 1, single (1), 0), single ([y2, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], single (1), 1, 0), single ([y2, NaN]), eps("single"))

%!assert (gpcdf (single ([x, NaN]), -1, 1, 0), single ([y3, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], -1, 1, single (0)), single ([y3, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], -1, single (1), 0), single ([y3, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], single (-1), 1, 0), single ([y3, NaN]), eps("single"))

## Test input validation
%!error gpcdf ()
%!error gpcdf (1, 2, 3, 4, 5, 6)
%!error gpcdf (ones (3), ones (2), ones (2), ones (2))
%!error gpcdf (ones (2), ones (2), ones (2), ones (3))
%!error gpcdf (ones (2), ones (2), ones (3), ones (2))
%!error gpcdf (ones (2), ones (3), ones (2), ones (2))
%!error gpcdf (i, 2, 2, 2)
%!error gpcdf (2, i, 2, 2)
%!error gpcdf (2, 2, i, 2)
%!error gpcdf (2, 2, 2, i)

