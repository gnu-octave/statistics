## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{r} =} stblrnd (@var{alpha}, @var{beta}, @var{gam}, @var{delta})
## @deftypefnx {statistics} {@var{r} =} stblrnd (@var{alpha}, @var{beta}, @var{gam}, @var{delta}, @var{m})
## @deftypefnx {statistics} {@var{r} =} stblrnd (@var{alpha}, @var{beta}, @var{gam}, @var{delta}, @var{m}, @var{n}, @dots{})
## @deftypefnx {statistics} {@var{r} =} stblrnd (@var{alpha}, @var{beta}, @var{gam}, @var{delta}, [@var{m}, @var{n}, @dots{}])
##
## Random arrays from the stable distribution.
##
## @code{@var{r} = stblrnd (@var{alpha}, @var{beta}, @var{gam}, @var{delta})}
## returns a random value drawn from the stable distribution with tail index
## (first shape parameter) @var{alpha}, skewness (second shape parameter)
## @var{beta}, scale parameter @var{gam}, and location parameter @var{delta}, in
## the Nolan @qcode{S0} parameterization.
##
## @var{alpha} must be in the range @math{(0, 2]}, @var{beta} in @math{[-1, 1]},
## @var{gam} positive, and @var{delta} real.  The parameters must be scalars.
##
## @code{stblrnd (@var{alpha}, @var{beta}, @var{gam}, @var{delta}, @var{m},
## @var{n}, @dots{})} or @code{stblrnd (@dots{}, [@var{m}, @var{n}, @dots{}])}
## returns an @var{m}-by-@var{n}-by-@dots{} array, following the size
## conventions of @code{rand}.
##
## The values are generated with the Chambers-Mallows-Stuck method.
##
## @seealso{stblpdf, stblcdf, stblinv, makedist}
## @end deftypefn

function r = stblrnd (alpha, beta, gam, delta, varargin)

  if (nargin < 4)
    print_usage ();
  endif

  msg = __stable_checkparams__ (alpha, beta, gam, delta);
  if (! isempty (msg))
    error ("stblrnd: %s", msg);
  endif

  ## Output size, following rand's conventions
  if (numel (varargin) == 0)
    sz = [1, 1];
  elseif (numel (varargin) == 1)
    a = varargin{1};
    if (isscalar (a))
      sz = [a, a];
    else
      sz = a(:)';
    endif
  else
    sz = [varargin{:}];
  endif
  if (! all (sz >= 0 & sz == fix (sz)))
    error ("stblrnd: dimensions must be non-negative integers.");
  endif

  ## Chambers-Mallows-Stuck: draw a uniform on (-pi/2, pi/2) and a unit
  ## exponential, form a standard S1 variate, then shift to S0 and scale.
  V = (rand (sz) - 0.5) .* pi;
  W = -log (rand (sz));

  if (alpha == 1)
    X = (2 ./ pi) .* ((pi ./ 2 + beta .* V) .* tan (V) ...
                      - beta .* log ((pi ./ 2 .* W .* cos (V)) ...
                                     ./ (pi ./ 2 + beta .* V)));
    r = gam .* X + delta;
  else
    ct = tan (pi .* alpha ./ 2);
    B = atan (beta .* ct) ./ alpha;
    S = (1 + beta .^ 2 .* ct .^ 2) .^ (1 ./ (2 .* alpha));
    X = S .* sin (alpha .* (V + B)) ./ (cos (V)) .^ (1 ./ alpha) ...
          .* (cos (V - alpha .* (V + B)) ./ W) .^ ((1 - alpha) ./ alpha);
    r = gam .* (X - beta .* ct) + delta;
  endif

endfunction

%!demo
%! ## Draw a large stable sample and overlay the theoretical density
%! r = stblrnd (1.5, 0.5, 1, 0, 1, 1e5);
%! r = r(abs (r) < 15);
%! hist (r, 100, 1);
%! hold on;
%! x = linspace (-15, 15, 400);
%! plot (x, stblpdf (x, 1.5, 0.5, 1, 0), "r-", "linewidth", 2);
%! hold off;

## The empirical cdf of a large sample matches stblcdf (seeded, repeatable)
%!test
%! rand ("state", 42);
%! r = stblrnd (1.5, 0.5, 1, 0, 1, 200000);
%! xs = [-2, -0.5, 0.5, 2];
%! ec = arrayfun (@(x) mean (r <= x), xs);
%! assert (ec, stblcdf (xs, 1.5, 0.5, 1, 0), 0.01);
%!test  # heavy-tailed alpha < 1, scaled and shifted
%! rand ("state", 7);
%! r = stblrnd (0.8, -0.3, 2, 1, 1, 200000);
%! xs = [-3, 0, 1, 4];
%! ec = arrayfun (@(x) mean (r <= x), xs);
%! assert (ec, stblcdf (xs, 0.8, -0.3, 2, 1), 0.01);
%!test  # alpha = 1 with beta ~= 0
%! rand ("state", 99);
%! r = stblrnd (1, 0.5, 1.5, -2, 1, 200000);
%! xs = [-5, -2, 0, 3];
%! ec = arrayfun (@(x) mean (r <= x), xs);
%! assert (ec, stblcdf (xs, 1, 0.5, 1.5, -2), 0.01);
%!test  # alpha = 2 is normal with variance 2*gam^2
%! rand ("state", 1);
%! r = stblrnd (2, 0, 1, 0, 1, 200000);
%! assert (mean (r), 0, 0.02);
%! assert (var (r), 2, 0.05);

## Size handling follows rand
%!test
%! assert (size (stblrnd (1.5, 0.5, 1, 0, 3, 4)), [3, 4]);
%! assert (size (stblrnd (1.5, 0.5, 1, 0, [2, 5])), [2, 5]);
%! assert (isscalar (stblrnd (1.5, 0.5, 1, 0)));

## Test input validation
%!error <Invalid call to stblrnd> stblrnd (1.5, 0.5, 1)
%!error <stblrnd: ALPHA must be a scalar in the range \(0, 2\].> ...
%! stblrnd (2.5, 0, 1, 0)
%!error <stblrnd: BETA must be a scalar in the range \[-1, 1\].> ...
%! stblrnd (1.5, 2, 1, 0)
%!error <stblrnd: GAM must be a positive scalar.> stblrnd (1.5, 0, 0, 0)
%!error <stblrnd: DELTA must be a real scalar.> stblrnd (1.5, 0, 1, 1i)
%!error <stblrnd: dimensions must be non-negative integers.> ...
%! stblrnd (1.5, 0, 1, 0, 2.5)
