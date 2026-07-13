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
## @deftypefn  {statistics} {@var{f} =} ksdensity (@var{x})
## @deftypefnx {statistics} {@var{f} =} ksdensity (@var{x}, @var{pts})
## @deftypefnx {statistics} {[@var{f}, @var{xi}] =} ksdensity (@dots{})
## @deftypefnx {statistics} {[@var{f}, @var{xi}, @var{bw}] =} ksdensity (@dots{})
## @deftypefnx {statistics} {[@dots{}] =} ksdensity (@dots{}, @var{Name}, @var{Value})
##
## Kernel smoothing density estimate.
##
## @code{@var{f} = ksdensity (@var{x})} computes a probability density estimate
## of the sample in the vector @var{x}, evaluated at 100 equally spaced points
## @var{xi} that span the range of the data.  @code{[@var{f}, @var{xi}] =
## ksdensity (@var{x})} also returns those points.  When called without output
## arguments, the estimate is plotted instead.
##
## @code{@var{f} = ksdensity (@var{x}, @var{pts})} evaluates the estimate at the
## values in @var{pts} instead; @var{f} is then the same size as @var{pts}.  For
## @qcode{'Function'} equal to @qcode{'icdf'} the entries of @var{pts} are
## probabilities in @math{[0, 1]}.
##
## @code{[@var{f}, @var{xi}, @var{bw}] = ksdensity (@dots{})} additionally
## returns the bandwidth @var{bw} of the smoothing kernel.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.18 0.82
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Kernel'} @tab The smoothing kernel: @qcode{'normal'} (default),
## @qcode{'box'}, @qcode{'triangle'}, @qcode{'epanechnikov'}, or a function
## handle @code{@@(z)} evaluating a kernel density at the standardized distance
## @var{z}.
##
## @item @qcode{'Bandwidth'} @tab The kernel bandwidth, a positive scalar.  The
## default is the value that is optimal for estimating a normal density,
## @math{@var{bw} = @var{sigma} * (4 / (3 * @var{n})) ^ (1 / 5)}, with
## @var{sigma} a robust estimate of the standard deviation of @var{x}.
##
## @item @qcode{'Function'} @tab The function to estimate: @qcode{'pdf'}
## (default), @qcode{'cdf'}, @qcode{'icdf'}, @qcode{'survivor'}, or
## @qcode{'cumhazard'}.
##
## @item @qcode{'Weights'} @tab A vector of non-negative weights, one for each
## element of @var{x}.  The default weights are all equal.
##
## @item @qcode{'NumPoints'} @tab The number of equally spaced points @var{xi}
## at which to evaluate the estimate when @var{pts} is not given.  The default
## is @math{100}.
## @end multitable
##
## @seealso{hist, histc, ecdf}
## @end deftypefn

function [f, xi, bw] = ksdensity (x, varargin)

  if (nargin < 1)
    print_usage ();
  endif
  if (! (isnumeric (x) && isreal (x) && isvector (x)))
    error ("ksdensity: X must be a vector of real values.");
  endif
  x = x(! isnan (x));
  x = x(:)';
  n = numel (x);
  if (n < 2)
    error ("ksdensity: X must contain at least two non-missing values.");
  endif

  ## Optional positional PTS argument (a numeric vector before any Name-Value).
  pts = [];
  ptssz = [];
  if (numel (varargin) >= 1 && ! ischar (varargin{1}))
    pts = varargin{1};
    varargin(1) = [];
    if (! (isnumeric (pts) && isreal (pts) && isvector (pts)))
      error ("ksdensity: PTS must be a vector of real values.");
    endif
    ptssz = size (pts);
    pts = pts(:);
  endif

  ## Defaults and Name-Value parsing.
  kernel = 'normal';
  bw = [];
  func = 'pdf';
  weights = [];
  npoints = 100;
  support = 'unbounded';
  bc = 'log';
  if (mod (numel (varargin), 2) != 0)
    error ("ksdensity: optional arguments must be Name-Value pairs.");
  endif
  for k = 1:2:numel (varargin)
    if (! ischar (varargin{k}))
      error ("ksdensity: parameter names must be character vectors.");
    endif
    switch (lower (varargin{k}))
      case 'kernel'
        kernel = varargin{k+1};
      case 'bandwidth'
        bw = varargin{k+1};
      case 'function'
        func = lower (varargin{k+1});
      case 'weights'
        weights = varargin{k+1};
      case 'numpoints'
        npoints = varargin{k+1};
      case 'support'
        support = varargin{k+1};
      case 'boundarycorrection'
        bc = lower (varargin{k+1});
      case 'censoring'
        error ("ksdensity: 'Censoring' is not yet supported.");
      otherwise
        error ("ksdensity: unknown parameter name '%s'.", varargin{k});
    endswitch
  endfor

  ## Validate the kernel.
  if (ischar (kernel))
    kernel = lower (kernel);
    if (! any (strcmp (kernel, {'normal', 'box', 'triangle', 'epanechnikov'})))
      error ("ksdensity: unrecognised 'Kernel' value.");
    endif
  elseif (! is_function_handle (kernel))
    error ("ksdensity: 'Kernel' must be a name or a function handle.");
  endif
  if (! any (strcmp (func, {'pdf', 'cdf', 'icdf', 'survivor', 'cumhazard'})))
    error ("ksdensity: unrecognised 'Function' value.");
  endif

  ## Support and boundary-correction validation -> lower/upper bounds L, U.
  if (ischar (support))
    switch (lower (support))
      case 'unbounded'
        L = -Inf;  U = Inf;
      case 'positive'
        L = 0;  U = Inf;
      otherwise
        error ("ksdensity: unrecognised 'Support' value.");
    endswitch
  elseif (isnumeric (support) && isreal (support) && numel (support) == 2
          && support(1) < support(2))
    L = support(1);  U = support(2);
  else
    error ("ksdensity: 'Support' must be 'unbounded', 'positive', or [L U].");
  endif
  unbounded = (L == -Inf && U == Inf);
  if (! any (strcmp (bc, {'log', 'reflection'})))
    error ("ksdensity: unrecognised 'BoundaryCorrection' value.");
  endif
  if (! unbounded && (any (x <= L) || any (x >= U)))
    error ("ksdensity: X must lie strictly inside the specified 'Support'.");
  endif

  ## Weights (normalised to sum to one).
  if (isempty (weights))
    w = ones (1, n) / n;
  else
    if (! (isnumeric (weights) && isreal (weights) && isvector (weights)
           && numel (weights) == n && all (weights >= 0) && any (weights > 0)))
      error ("ksdensity: 'Weights' must be a non-negative vector with one element per X.");
    endif
    w = weights(:)' / sum (weights);
  endif

  ## Working space.  With the default 'log' boundary correction a bounded
  ## support is mapped to the whole real line by a log/logit transform; the
  ## bandwidth and default grid are then formed in that space.  The unbounded
  ## case and the 'reflection' correction work in the native space.
  dolog = (! unbounded && strcmp (bc, 'log'));
  if (dolog)
    xw = ksdensity_fwd_ (x, L, U);
  else
    xw = x;
  endif

  ## Bandwidth: robust normal-reference rule unless supplied (working space).
  if (isempty (bw))
    sigma = median (abs (xw - median (xw))) / 0.6745;
    if (sigma <= 0)
      sigma = std (xw);
    endif
    if (sigma <= 0)
      sigma = 1;
    endif
    bw = sigma * (4 / (3 * n)) ^ (1 / 5);
  elseif (! (isnumeric (bw) && isscalar (bw) && isreal (bw) && bw > 0))
    error ("ksdensity: 'Bandwidth' must be a positive scalar.");
  endif

  ## Evaluation points.  For 'icdf' the points are probabilities.
  if (strcmp (func, 'icdf'))
    if (isempty (pts))
      xi = (((1:npoints) - 0.5) / npoints)';
    else
      xi = pts;
    endif
    f = ksdensity_icdf_ (xi(:), x, w, bw, kernel, L, U, bc);
  else
    if (isempty (pts))
      if (dolog)
        gw = linspace (min (xw) - 3 * bw, max (xw) + 3 * bw, npoints);
        xi = ksdensity_inv_ (gw, L, U)';
      else
        lo = min (x) - 3 * bw;
        hi = max (x) + 3 * bw;
        if (isfinite (L))
          lo = max (lo, L);
        endif
        if (isfinite (U))
          hi = min (hi, U);
        endif
        xi = linspace (lo, hi, npoints)';
      endif
    else
      xi = pts;
    endif
    switch (func)
      case 'pdf'
        f = ksdensity_eval_ (xi(:), x, w, bw, kernel, L, U, bc, 'pdf');
      case 'cdf'
        f = ksdensity_eval_ (xi(:), x, w, bw, kernel, L, U, bc, 'cdf');
      case 'survivor'
        f = 1 - ksdensity_eval_ (xi(:), x, w, bw, kernel, L, U, bc, 'cdf');
      case 'cumhazard'
        Fc = ksdensity_eval_ (xi(:), x, w, bw, kernel, L, U, bc, 'cdf');
        f = -log (1 - Fc);
    endswitch
  endif

  ## Match the orientation of a supplied PTS argument.
  if (! isempty (pts))
    f = reshape (f, ptssz);
    xi = reshape (xi, ptssz);
  endif

  ## With no output requested, plot the estimate (as MATLAB does) and return
  ## nothing.
  if (nargout == 0)
    plot (xi, f);
    clear f xi bw
  endif

endfunction

## Kernel pdf or cdf estimate at native query points Q, honouring the support
## bounds L, U and the boundary-correction method BC.  WANT is 'pdf' or 'cdf'.
function v = ksdensity_eval_ (q, x, w, h, kernel, L, U, bc, want)
  q = q(:);
  ispdf = strcmp (want, 'pdf');
  if (L == -Inf && U == Inf)
    ## Unbounded: plain kernel sum.
    z = (q - x) / h;
    if (ispdf)
      v = (kernelpdf (z, kernel) * w(:)) / h;
    else
      v = kernelcdf (z, kernel) * w(:);
    endif
  elseif (strcmp (bc, 'reflection'))
    ## Reflect the sample across each finite bound (single reflection).
    a = x;  aw = w;
    if (isfinite (L))
      a = [a, 2 * L - x];  aw = [aw, w];
    endif
    if (isfinite (U))
      a = [a, 2 * U - x];  aw = [aw, w];
    endif
    if (ispdf)
      v = (kernelpdf ((q - a) / h, kernel) * aw(:)) / h;
      v(q < L | q > U) = 0;
    else
      ## Integrate the reflected density from the lower edge of the support.
      v = (kernelcdf ((q - a) / h, kernel) ...
           - kernelcdf ((L - a) / h, kernel)) * aw(:);
    endif
  else
    ## Log/logit transform: estimate in the transformed space and map back.
    z = (ksdensity_fwd_ (q, L, U) - ksdensity_fwd_ (x, L, U)) / h;
    if (ispdf)
      g = (kernelpdf (z, kernel) * w(:)) / h;
      v = g .* ksdensity_logjac_ (q, L, U);
    elseif (isfinite (U) && ! isfinite (L))    ## upper bound only: t decreases
      v = 1 - kernelcdf (z, kernel) * w(:);
    else
      v = kernelcdf (z, kernel) * w(:);
    endif
  endif
  if (! ispdf)
    v = min (max (v, 0), 1);
  endif
endfunction

## Forward transform mapping the support (L, U) onto the whole real line.
function t = ksdensity_fwd_ (x, L, U)
  if (U == Inf)
    t = log (x - L);
  elseif (L == -Inf)
    t = log (U - x);
  else
    t = log ((x - L) ./ (U - x));
  endif
endfunction

## Inverse of ksdensity_fwd_.
function x = ksdensity_inv_ (t, L, U)
  if (U == Inf)
    x = L + exp (t);
  elseif (L == -Inf)
    x = U - exp (t);
  else
    x = L + (U - L) ./ (1 + exp (-t));
  endif
endfunction

## |dt/dx| of the forward transform (the density change-of-variables factor).
function j = ksdensity_logjac_ (x, L, U)
  if (U == Inf)
    j = 1 ./ (x - L);
  elseif (L == -Inf)
    j = 1 ./ (U - x);
  else
    j = (U - L) ./ ((x - L) .* (U - x));
  endif
endfunction

## Inverse cdf at probabilities P by monotone interpolation of the cdf over a
## fine grid spanning the support.
function q = ksdensity_icdf_ (p, x, w, h, kernel, L, U, bc)
  if (! (L == -Inf && U == Inf) && strcmp (bc, 'log'))
    tx = ksdensity_fwd_ (x, L, U);
    grid = ksdensity_inv_ (linspace (min (tx) - 10 * h, ...
                                     max (tx) + 10 * h, 4000), L, U)';
  else
    lo = min (x) - 10 * h;
    hi = max (x) + 10 * h;
    if (isfinite (L))
      lo = max (lo, L);
    endif
    if (isfinite (U))
      hi = min (hi, U);
    endif
    grid = linspace (lo, hi, 4000)';
  endif
  F = ksdensity_eval_ (grid, x, w, h, kernel, L, U, bc, 'cdf');
  ## Keep a strictly increasing (F, grid) relation for interp1.
  dF = diff (F);
  keep = [true; dF > 0];
  q = interp1 (F(keep), grid(keep), p(:), 'linear', NA);
  q(p(:) <= min (F)) = grid(1);
  q(p(:) >= max (F)) = grid(end);
endfunction

%!demo
%! ## Kernel density estimate of a small sample, with a histogram for reference
%! x = [1 1.5 2 2 2.5 3 3.5 3.5 4 6];
%! [f, xi] = ksdensity (x);
%! hist (x, 6, 6 / numel (x));
%! hold on;  plot (xi, f, 'r-', 'LineWidth', 2);  hold off;

%!test  # density integrates to ~1 over a wide grid (normal kernel)
%! x = [2.1 0.3 1.2 -0.7 0.9 1.5 2.8 0.1 0.4 1.1 3.2 0.6 2.0 0.9 1.7];
%! [f, xi] = ksdensity (x, "NumPoints", 4000);
%! assert_equal (trapz (xi, f), 1, 5e-3);

%!test  # every named compact kernel integrates to ~1
%! x = randn (1, 200);
%! for k = {"box", "triangle", "epanechnikov"}
%!   xi = linspace (-8, 8, 6000)';
%!   f = ksdensity (x, xi, "Kernel", k{1});
%!   assert_equal (trapz (xi, f), 1, 1e-2);
%! endfor

%!test  # cdf is monotone from 0 to 1 and matches the analytic normal-kernel sum
%! x = [2.1 0.3 1.2 -0.7 0.9 1.5 2.8 0.1 0.4 1.1 3.2 0.6 2.0 0.9 1.7];
%! [F, xi] = ksdensity (x, "Function", "cdf", "NumPoints", 500);
%! assert_equal (all (diff (F) >= -1e-12), true);
%! assert_equal (F(1), 0, 5e-3);
%! assert_equal (F(end), 1, 5e-3);
%! [~, ~, bw] = ksdensity (x);
%! Fdirect = mean (normcdf ((xi - x) / bw), 2);
%! assert_equal (F, Fdirect, 1e-12);

%!test  # survivor and cumhazard are consistent with the cdf
%! x = [2.1 0.3 1.2 -0.7 0.9 1.5 2.8 0.1 0.4 1.1 3.2 0.6 2.0 0.9 1.7];
%! pts = [-1 0 1 2 3]';
%! F = ksdensity (x, pts, "Function", "cdf");
%! S = ksdensity (x, pts, "Function", "survivor");
%! H = ksdensity (x, pts, "Function", "cumhazard");
%! assert_equal (S, 1 - F, 1e-12);
%! assert_equal (H, -log (1 - F), 1e-12);

%!test  # icdf inverts the cdf
%! x = [2.1 0.3 1.2 -0.7 0.9 1.5 2.8 0.1 0.4 1.1 3.2 0.6 2.0 0.9 1.7];
%! p = [0.1 0.25 0.5 0.75 0.9]';
%! q = ksdensity (x, p, "Function", "icdf");
%! Fq = ksdensity (x, q, "Function", "cdf");
%! assert_equal (Fq, p, 5e-3);

%!test  # evaluation at supplied points preserves shape
%! x = randn (1, 50);
%! pts = [-1 0 1];
%! f = ksdensity (x, pts);
%! assert_equal (size (f), size (pts));

%!test  # weights: a duplicated point equals a doubled weight
%! x = [0 1 2 3];
%! pts = linspace (-2, 5, 40)';
%! f1 = ksdensity ([x, 3], pts, "Bandwidth", 0.5);
%! f2 = ksdensity (x, pts, "Bandwidth", 0.5, "Weights", [1 1 1 2]);
%! assert_equal (f1, f2, 1e-12);

%!test  # MATLAB parity: default bandwidth, grid size and range
%! x = [2.1 0.3 1.2 -0.7 0.9 1.5 2.8 0.1 0.4 1.1 3.2 0.6 2.0 0.9 1.7];
%! [~, xi, bw] = ksdensity (x);
%! assert_equal (bw, 0.6396, 5e-4);
%! assert_equal (numel (xi), 100);
%! assert_equal ([xi(1), xi(end)], [-2.6187, 5.1187], 1e-3);

%!test  # MATLAB parity: pdf (normal and box kernels) and cdf at fixed points
%! x = [2.1 0.3 1.2 -0.7 0.9 1.5 2.8 0.1 0.4 1.1 3.2 0.6 2.0 0.9 1.7];
%! pts = [-0.5 0 0.5 1 1.5 2 2.5 3];
%! assert_equal (ksdensity (x, pts), ...
%!         [0.1214 0.2141 0.3051 0.3394 0.3061 0.2375 0.1697 0.1166], 2e-3);
%! assert_equal (ksdensity (x, pts, "Kernel", "box"), ...
%!         [0.1505 0.2407 0.2708 0.3611 0.3009 0.2708 0.1805 0.1204], 2e-3);
%! assert_equal (ksdensity (x, pts, "Function", "cdf"), ...
%!         [0.0710 0.1538 0.2850 0.4492 0.6128 0.7494 0.8506 0.9217], 2e-3);

%!test  # MATLAB parity: positive support (log) bandwidth, grid and pdf
%! y = [0.2 0.5 0.7 1.1 1.4 2 2.6 3.3 4.1 5.5];
%! ypts = [0.1 0.3 0.6 1 2 3 4 5];
%! [~, xy, by] = ksdensity (y, "Support", "positive");
%! assert_equal (by, 0.7682, 5e-4);
%! assert_equal ([xy(1), xy(end)], [0.0200, 55.1111], 1e-3);
%! assert_equal (ksdensity (y, ypts, "Support", "positive"), ...
%!         [0.4300 0.4620 0.3620 0.2740 0.1570 0.0990 0.0660 0.0450], 2e-3);

%!test  # MATLAB parity: positive support with reflection boundary correction
%! y = [0.2 0.5 0.7 1.1 1.4 2 2.6 3.3 4.1 5.5];
%! ypts = [0.1 0.3 0.6 1 2 3 4 5];
%! f = ksdensity (y, ypts, "Support", "positive", ...
%!                "BoundaryCorrection", "reflection");
%! assert_equal (f, [0.2920 0.2900 0.2810 0.2620 0.2010 0.1470 0.1070 0.0740], 2e-3);

%!test  # bounded support integrates to ~1 and vanishes outside the bounds
%! y = [0.2 0.5 0.7 1.1 1.4 2 2.6 3.3];
%! [f, xi] = ksdensity (y, "Support", [0 4], "NumPoints", 4000);
%! assert_equal (trapz (xi, f), 1, 5e-3);
%! assert_equal (all (xi >= 0 & xi <= 4), true);

%!test  # reflection density is confined to the support
%! y = [0.2 0.5 0.7 1.1 1.4 2 2.6 3.3];
%! f = ksdensity (y, [-1 -0.1 5], "Support", "positive", ...
%!                "BoundaryCorrection", "reflection");
%! assert_equal (f(1:2), [0 0]);

## Test input validation
%!error <Invalid call to ksdensity> ksdensity ()
%!error <ksdensity: X must be a vector of real values.> ksdensity (ones (3, 3))
%!error <ksdensity: X must contain at least two non-missing values.> ...
%! ksdensity (5)
%!error <ksdensity: unrecognised 'Kernel' value.> ...
%! ksdensity (1:10, "Kernel", "cosine")
%!error <ksdensity: unrecognised 'Function' value.> ...
%! ksdensity (1:10, "Function", "hazard")
%!error <ksdensity: 'Bandwidth' must be a positive scalar.> ...
%! ksdensity (1:10, "Bandwidth", -1)
%!error <ksdensity: 'Weights' must be a non-negative vector with one element per X.> ...
%! ksdensity (1:10, "Weights", [1 2 3])
%!error <ksdensity: unrecognised 'Support' value.> ...
%! ksdensity (1:10, "Support", "half")
%!error <ksdensity: 'Support' must be 'unbounded', 'positive', or .L U..> ...
%! ksdensity (1:10, "Support", [2 1])
%!error <ksdensity: X must lie strictly inside the specified 'Support'.> ...
%! ksdensity ([-1 1 2 3], "Support", "positive")
%!error <ksdensity: unrecognised 'BoundaryCorrection' value.> ...
%! ksdensity (1:10, "Support", "positive", "BoundaryCorrection", "linear")
%!error <ksdensity: 'Censoring' is not yet supported.> ...
%! ksdensity (1:10, "Censoring", ones (1, 10))
