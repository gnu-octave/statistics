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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Private Function} {[@var{F}, @var{U}] =} __stdr__ (@var{x}, @var{k}, @var{df}, @qcode{'cdf'})
## @deftypefnx {Private Function} {@var{y} =} __stdr__ (@var{x}, @var{k}, @var{df}, @qcode{'pdf'})
##
## Studentized range distribution at positive finite values.
##
## @var{x} is a vector of positive finite values, @var{k} a scalar integer of
## at least 2 and @var{df} a positive scalar, possibly @code{Inf}.  With
## @qcode{'cdf'}, @var{F} and @var{U} are the lower and upper tail
## probabilities, each computed directly so that neither loses its relative
## accuracy to the other; with @qcode{'pdf'}, @var{y} is the density.
##
## The range @math{W} of @var{k} standard normals has
## @math{P(W <= w) = k int phi(z) (Phi(z) - Phi(z - w))^(k-1) dz}, evaluated on
## a fixed composite Gauss-Legendre grid in @math{z}.  The studentized range is
## @math{W / s} with @math{s^2} a chi-square over @var{df}, so a finite
## @var{df} adds an adaptive integral over @math{log s}, weighted by the
## density of @math{s}.
## @end deftypefn

function [F, U] = __stdr__ (x, k, df, what)

  x = x(:)';
  if (strcmp (what, 'pdf'))
    if (isinf (df))
      F = rangepdf (x, k);
    else
      F = zeros (size (x));
      for i = 1:numel (x)
        F(i) = outer (@(w) rangepdf (w, k), x(i), k, df, 1);
      endfor
    endif
    return;
  endif

  if (isinf (df))
    [F, U] = rangecdf (x, k);
    return;
  endif

  ## Each tail is integrated on its own; the larger one is the complement of
  ## the smaller, which is where the accuracy matters.
  F = zeros (size (x));
  U = zeros (size (x));
  for i = 1:numel (x)
    U(i) = outer (@(w) nthargout (2, @rangecdf, w, k), x(i), k, df, 0);
    if (U(i) < 0.5)
      F(i) = 1 - U(i);
    else
      F(i) = outer (@(w) rangecdf (w, k), x(i), k, df, 0);
      U(i) = 1 - F(i);
    endif
  endfor

endfunction

## The integral over t = log (s) of H (x * s) weighted by the density of s,
## times s ^ P.  Its integrand peaks near t = 0, where the density of s does,
## near 0.5 * log (df / (df + x^2 / 2)) for the upper tail and near
## 0.5 * log ((df + k - 1) / df) for the lower tail, and no feature of it is
## narrower than 1 / sqrt (2 * df); the window reaches far enough past both
## peaks for the density of s to fall below exp (-40).  A fixed grid of
## panels half that width is exact to rounding, where an adaptive rule loses
## the peak once df is large.
function I = outer (H, x, k, df, P)

  ## The density of t is exp (C + df * (t - expm1 (2 * t) / 2)), written so
  ## that no two large terms cancel inside the integrand, and C from
  ## Stirling's series once df is large enough for the direct form to lose
  ## digits.
  a = df / 2;
  if (a > 1e3)
    C = log (2) + 0.5 * log (a / (2 * pi)) - 1 / (12 * a) ...
        + 1 / (360 * a ^ 3) - 1 / (1260 * a ^ 5);
  else
    C = log (2) + a * log (a) - gammaln (a) - a;
  endif
  sd = 1 / sqrt (2 * df);
  tu = 0.5 * log (df / (df + x ^ 2 / 2));
  tl = 0.5 * log ((df + k - 1) / df);
  ## Left of both peaks the integrand is exp (df * t) times a constant, so
  ## panels there need only be 2 / df wide
  mid = min (tu, 0) - 12 * sd;
  lo = mid - 40 / df;
  hi = max (tl + 12 * sd, ...
            0.5 * log ((df + k + 40 * sqrt (2 * (df + k)) + 80) / df));
  [t1, w1] = legendre ((mid - lo) / max (sd / 2, 2 / df), lo, mid);
  [t2, w2] = legendre ((hi - mid) / (sd / 2), mid, hi);
  t = [t1, t2];
  wt = [w1, w2];
  I = sum (wt .* exp (C + df * (t - expm1 (2 * t) / 2) + P * t) ...
           .* H (x * exp (t)));

endfunction

## Nodes and weights of a composite 16-point Gauss-Legendre rule over
## [LO, HI] in at least N panels, as row vectors.
function [t, wt] = legendre (n, lo, hi)

  persistent xg wg;
  if (isempty (xg))
    m = 16;
    b = 0.5 ./ sqrt (1 - (2 * (1:m-1)) .^ (-2));
    [V, D] = eig (diag (b, 1) + diag (b, -1));
    xg = diag (D);
    wg = 2 * V(1,:)' .^ 2;
  endif
  np = max (1, ceil (n));
  h = (hi - lo) / np;
  c = lo + h * ((0:np-1) + 0.5);
  t = reshape (c + (h / 2) * xg, 1, []);
  wt = repmat ((h / 2) * wg', 1, np);

endfunction

## Lower and upper tail probabilities of the range of K standard normals.
function [F, U] = rangecdf (w, k)

  [z, wz, phi, a, Pz, Qz, Pzw] = grid (w, k);
  F = k * sum (wz .* phi .* a .^ (k - 1), 1);
  ## Phi(z)^(k-1) - a^(k-1), written so that no two close numbers are
  ## subtracted: a = Phi(z) - Phi(z-w).
  ## The ratio is at most 1 in exact arithmetic, but where w is tiny against z
  ## a libm whose erfc is not monotone in the last bit can push it above,
  ## and log1p would then turn complex.
  logPz = log (Pz);
  logPz(z > 0) = log1p (- Qz(z > 0));
  r = min (Pzw ./ Pz, 1);
  d = exp ((k - 1) * logPz) .* (- expm1 ((k - 1) * log1p (- r)));
  U = k * sum (wz .* phi .* d, 1);
  ## Beyond w = 40 the upper tail is below exp (-400).
  far = (w > 40);
  F(far) = 1;
  U(far) = 0;
  F = reshape (F, size (w));
  U = reshape (U, size (w));

endfunction

## Density of the range of K standard normals.
function f = rangepdf (w, k)

  [z, wz, phi, a] = grid (w, k);
  phizw = exp (- (z - w(:)') .^ 2 / 2) / sqrt (2 * pi);
  f = k * (k - 1) * sum (wz .* phi .* phizw .* a .^ (k - 2), 1);
  f(w > 40) = 0;
  f = reshape (f, size (w));

endfunction

## The z grid, one column of nodes shared by every W, and the quantities
## on it.  Panels are at most 4 / sqrt (k) wide, the width of the peak
## the integrands have for large K.
function [z, wz, phi, a, Pz, Qz, Pzw] = grid (w, k)

  persistent xg wg;
  if (isempty (xg))
    n = 16;
    b = 0.5 ./ sqrt (1 - (2 * (1:n-1)) .^ (-2));
    [V, D] = eig (diag (b, 1) + diag (b, -1));
    xg = diag (D);
    wg = 2 * V(1,:)' .^ 2;
  endif

  w = w(:)';
  zlo = -9;
  zhi = min (29, max (9, max (w(w <= 40)) / 2 + 9));
  if (isempty (zhi))
    zhi = 9;
  endif
  np = ceil ((zhi - zlo) / min (1, 4 / sqrt (k)));
  h = (zhi - zlo) / np;
  c = zlo + h * ((0:np-1) + 0.5);
  z = reshape (c + (h / 2) * xg, [], 1);
  wz = repmat ((h / 2) * wg, np, 1);
  phi = exp (- z .^ 2 / 2) / sqrt (2 * pi);
  Pz = 0.5 * erfc (- z / sqrt (2));
  Qz = 0.5 * erfc (z / sqrt (2));
  zw = z - w;
  Pzw = 0.5 * erfc (- zw / sqrt (2));
  ## a = Phi(z) - Phi(z-w) from whichever tail keeps it accurate
  a = Pz - Pzw;
  r = (zw > 0);
  Qzw = 0.5 * erfc (zw / sqrt (2));
  Qzm = repmat (Qz, 1, numel (w));
  a(r) = Qzw(r) - Qzm(r);
  s = (z > 0) & (zw <= 0);
  E = erf (z / sqrt (2)) - erf (zw / sqrt (2));
  a(s) = 0.5 * E(s);
  a = max (a, 0);

endfunction
