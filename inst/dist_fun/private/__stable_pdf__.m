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
## @deftypefn {Function File} {@var{y} =} __stable_pdf__ (@var{x}, @var{alpha}, @var{beta}, @var{gam}, @var{delta})
##
## Fast stable probability density for the fitting functions.
##
## Computes the stable density in the Nolan @qcode{S0} parameterization at all
## the points in @var{x} by a single shared-grid inversion of the characteristic
## function.  Unlike @code{stblpdf}, which invokes an adaptive @code{quadgk} per
## point, the characteristic function is evaluated once on a common quadrature
## grid and combined against every data point at once, which is roughly an order
## of magnitude faster.  This makes it suitable as the inner density of the
## maximum-likelihood loop in @code{stblfit}/@code{stbllike}, where accuracy of a
## few parts in @math{10^{-8}} across the bulk is ample for the estimator.
##
## The closed forms are used for the normal (@var{alpha} equal to @code{2}) and
## the Cauchy (@var{alpha} equal to @code{1} with @var{beta} equal to @code{0})
## special cases.  The parameters must be valid scalars; no checking is done
## here, as this is a private helper.
##
## @seealso{stblpdf, stblfit, stbllike}
## @end deftypefn

function y = __stable_pdf__ (x, alpha, beta, gam, delta)

  z = (x - delta) ./ gam;
  y = nan (size (z));
  ok = ! isnan (z);
  zk = z(ok)(:);

  if (isempty (zk))
    return;
  endif

  if (alpha == 2)
    ## Normal with variance 2
    v = exp (-zk .^ 2 ./ 4) ./ (2 .* sqrt (pi));
  elseif (alpha == 1 && beta == 0)
    ## Cauchy
    v = 1 ./ (pi .* (1 + zk .^ 2));
  else
    v = cf_invert (zk, alpha, beta);
  endif

  y(ok) = max (v, 0) ./ gam;

endfunction

## Density of the standard stable S(alpha, beta, 1, 0) at the points Z by
## Gil-Pelaez inversion of the characteristic function on a shared grid.
function v = cf_invert (z, alpha, beta)

  ## Upper limit: |phi(t)| = exp (-t^alpha), so t^alpha = L makes it negligible.
  L = 30;                          # exp (-30) ~ 9e-14
  T = L ^ (1 / alpha);

  ## Uniform grid.  The inversion integrand is the slowly decaying envelope
  ## |phi(t)| = exp (-t^alpha) times the oscillation exp (-i t z), whose
  ## frequency |z| is constant in t; a uniform mesh therefore resolves it most
  ## efficiently.  The node count scales with the peak total phase over [0, T],
  ## i.e. with the largest |z| in the data, with a floor that also covers the
  ## mild t^alpha cusp at the origin when alpha < 1.  The driving |z| is capped:
  ## far-tail points carry negligible density, so leaving their oscillation
  ## under-resolved keeps N bounded for heavy-tailed samples without affecting
  ## the likelihood.
  zmax = min (max (abs (z)), 30);
  N = ceil (6 * (zmax + 1) * T) + 512;
  N = min (N, 200000);

  t = linspace (0, T, N)';
  phi = __stable_cf__ (t, alpha, beta);

  ## Uniform trapezoidal weights
  dt = T / (N - 1);
  w = dt .* ones (N, 1);
  w(1) /= 2;
  w(end) /= 2;

  ## Vectorized inversion: real (exp (-i t z') .* phi) integrated over t, per z
  M = real (exp (-1i .* (t * z(:)')) .* phi);
  v = (w' * M)' ./ pi;

endfunction
