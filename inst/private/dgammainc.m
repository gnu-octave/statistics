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
## @deftypefn {Private Function} [@var{y}, @var{dy}, @var{d2y] = exact2xkCT (@var{x}, @var{k})
## @deftypefn {Private Function} [@var{y}, @var{dy}, @var{d2y] = exact2xkCT (@var{x}, @var{k}, @qcode{"upper"})
##
## Compute the incomplete Gamma function along with its first and second
## derivatives.
##
## @end deftypefn

function [y, dy, d2y] = dgammainc (x, k, uflag)

  ## Force X and K to common size
  [err, x, k] = common_size (x, k);
  if (err > 0)
    error ("dgammainc: X and K must be of common size or scalars.");
  endif

  ## Check for negative K
  if (any (k) < 0)
    error ("dgammainc: K cannot be negative.");
  endif

  ## Initialize return variables
  y = nan (size (x));
  dy = y;
  d2y = y;

  ## Check for valid "upper" flag
  if (nargin > 2)
    if (! strcmpi (uflag, "upper"))
      error ("dgammainc: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Set upper limit for series and continued fraction
  ulim = 2^20;

  ## Use approximation for K > 2^20
  is_lim = find (k > ulim);
  if (! isempty (is_lim))
    x(is_lim) = max (ulim - 1/3 + sqrt (ulim ./ k(is_lim)) .* ...
                     (x(is_lim) - (k(is_lim) - 1/3)), 0);
    k(is_lim) = ulim;
  endif

  ## For x < k+1
  is_lo = find(x < k + 1 & x != 0);
  if (! isempty (is_lo))
    x_lo = x(is_lo);
    k_lo = k(is_lo);
    k_1 = k_lo;
    step = 1;
    d1st = 0;
    d2st = 0;
    stsum = step;
    d1sum = d1st;
    d2sum = d2st;
    while norm (step, "inf") >= 100 * eps (norm (stsum, "inf"))
      k_1 += 1;
      step = step .* x_lo ./ k_1;
      d1st = (d1st .* x_lo - step) ./ k_1;
      d2st = (d2st .* x_lo - 2 .* d1st) ./ k_1;
      stsum = stsum + step;
      d1sum = d1sum + d1st;
      d2sum = d2sum + d2st;
    endwhile
    fklo = exp (-x_lo + k_lo .* log (x_lo) - gammaln (k_lo + 1));
    y_lo = fklo .* stsum;
    ## Fix very small k
    y_lo(x_lo > 0 & y_lo > 1) = 1;
    ## Compute 1st derivative
    dlogfklo = (log(x_lo) - psi(k_lo+1));
    d1fklo = fklo .* dlogfklo;
    d1y_lo = d1fklo.*stsum + fklo.*d1sum;
    ## Compute 2nd derivative
    d2fklo = d1fklo.*dlogfklo - fklo.*psi(1,k_lo+1);
    d2y_lo = d2fklo.*stsum + 2.*d1fklo.*d1sum + fklo.*d2sum;
    if (uflag)
      y(is_lo) = 1 - y_lo;
      dy(is_lo) = -d1y_lo;
      d2y(is_lo) = -d2y_lo;
    else
      y(is_lo) = y_lo;
      dy(is_lo) = d1y_lo;
      d2y(is_lo) = d2y_lo;
    endif
  endif

  ## For x >= k+1
  is_hi = find(x >= k+1); % & x ~= 0
  if ~isempty(is_hi)
    x_hi = x(is_hi);
    k_hi = k(is_hi);
    zc = 0;
    k0 = 0;
    k1 = k_hi;
    x0 = 1;
    x1 = x_hi;
    d1k0 = 0;
    d1k1 = 1;
    d1x0 = 0;
    d1x1 = 0;
    d2k0 = 0;
    d2k1 = 0;
    d2x0 = 0;
    d2x2 = 0;
    kx = k_hi ./ x_hi;
    d1kx = 1 ./ x_hi;
    d2kx = 0;
    start = 1;
    while norm (d2kx - start, "Inf") > 100 * eps (norm (d2kx, "Inf"))
      rescale = 1 ./ x1;
      zc += 1;
      nminusa = zc - k_hi;
      d2k0 = (d2k1 + d2k0 .* nminusa - 2 .* d1k0) .* rescale;
      d2x0 = (d2x2 + d2x0 .* nminusa - 2 .* d1x0) .* rescale;
      d1k0 = (d1k1 + d1k0 .* nminusa - k0) .* rescale;
      d1x0 = (d1x1 + d1x0 .* nminusa - x0) .* rescale;
      k0 = (k1 + k0 .* nminusa) .* rescale;
      x0 = 1 + (x0 .* nminusa) .* rescale;
      nrescale = zc .* rescale;
      d2k1 = d2k0 .* x_hi + d2k1 .* nrescale;
      d2x2 = d2x0 .* x_hi + d2x2 .* nrescale;
      d1k1 = d1k0 .* x_hi + d1k1 .* nrescale;
      d1x1 = d1x0 .* x_hi + d1x1 .* nrescale;
      k1 = k0 .* x_hi + k1 .* nrescale;
      x1 = x0 .* x_hi + zc;
      start = d2kx;
      kx = k1 ./ x1;
      d1kx = (d1k1 - kx.*d1x1) ./ x1;
      d2kx = (d2k1 - d1kx.*d1x1 - kx.*d2x2 - d1kx.*d1x1) ./ x1;
    endwhile
    fkhi = exp(-x_hi + k_hi.*log(x_hi) - gammaln(k_hi+1));
    y_hi = fkhi.*kx;
    ## Compute 1st derivative
    dlogfkhi = (log(x_hi) - psi(k_hi+1));
    d1fkhi = fkhi .* dlogfkhi;
    d1y_hi = d1fkhi.*kx + fkhi.*d1kx;
    ## Compute 2nd derivative
    d2fkhi = d1fkhi.*dlogfkhi - fkhi.*psi(1,k_hi+1);
    d2y_hi = d2fkhi.*kx + 2.*d1fkhi.*d1kx + fkhi.*d2kx;
    if (uflag)
      y(is_hi) = y_hi;
      dy(is_hi) = d1y_hi;
      d2y(is_hi) = d2y_hi;
    else
      y(is_hi) = 1 - y_hi;
      dy(is_hi) = -d1y_hi;
      d2y(is_hi) = -d2y_hi;
    endif
  endif

  ## Handle x == 0
  is_x0 = find(x == 0);
  if (! isempty (is_x0))
    if uflag
      y(is_x0) = 1;
    else
      y(is_x0) = 0;
    endif
    dy(is_x0) = 0;
    d2y(is_x0) = 0;
  endif

  ## Handle k == 0
  is_k0 = find(k == 0);
  if (! isempty (is_k0))
    is_k0x0 = find(k == 0 & x == 0);
    if (uflag)
      y(is_k0) = 0;
      dy(is_k0x0) = Inf;
      d2y(is_k0x0) = -Inf;
    else
      y(is_k0) = 1;
      dy(is_k0x0) = -Inf;
      d2y(is_k0x0) = Inf;
    endif
  endif
endfunction
