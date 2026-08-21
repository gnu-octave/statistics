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

## Build, once per predictor, a basis of the spline space that every round of a
## GAM fit works in.  A round fits the same X at the same knots and order and
## varies only the response, so the design is fixed: TMPL is a piecewise
## polynomial of that space with zero coefficients, and CACHE holds the basis
## values at X, the basis coefficients and a factorisation of the design, which
## reduce a round to two small products.  CACHE is empty where X does not
## determine the space, and the round has to be fitted outright.
function [tmpl, cache] = splineCache (x, brk, ord)

  n = numel (x);
  tmpl = splinefit (x, zeros (n, 1), brk, 'order', ord);
  tmpl.coefs(:) = 0;
  cache = [];

  ## Dimension of the spline space for these breaks and this order
  m = tmpl.pieces + tmpl.order - 1;
  if (n < m)
    return;
  endif

  ## Fit m independent responses in one call, which builds the basis once, and
  ## take the fitted splines as a basis of the space.  Any response is fitted
  ## by least squares onto the same space, so its fit is the projection onto
  ## the span of these, and its piecewise polynomial the same combination of
  ## theirs.
  Yb = cos ((0:m-1)' * pi * ((1:n) - 0.5) / n);
  ppb = splinefit (x, Yb, brk, 'order', ord);
  F = ppval (ppb, x)';
  if (rank (F) < m)
    return;
  endif

  [Q, R] = qr (F, 0);
  cache = struct ('F', F, 'Q', Q, 'R', R, ...
                  'C', reshape (ppb.coefs, [m, ppb.pieces, ppb.order]));

endfunction
