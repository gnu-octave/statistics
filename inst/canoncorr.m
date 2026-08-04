## Copyright (C) 2016-2019 by Nir Krakauer <mail@nirkrakauer.net>
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
## @deftypefn  {statistics} {[@var{A}, @var{B}, @var{r}, @var{U}, @var{V}] =} canoncorr (@var{X}, @var{Y})
##
## Canonical correlation analysis.
##
## Given @var{X} (size @var{k}*@var{m}) and @var{Y} (@var{k}*@var{n}), returns
## projection matrices of canonical coefficients @var{A} (size @var{m}*@var{d},
## where @var{d} is the smallest of @var{m}, @var{n}, @var{d}) and @var{B}
## (size @var{m}*@var{d}); the canonical correlations @var{r} (1*@var{d},
## arranged in decreasing order); the canonical variables @var{U}, @var{V}
## (both @var{k}*@var{d}, with orthonormal columns); and @var{stats},
## a structure containing results from Bartlett's chi-square and Rao's F tests
## of significance.
##
## @seealso{princomp}
## @end deftypefn

function [A,B,r,U,V,stats] = canoncorr (X,Y)

  k = size (X, 1); # should also be size (Y, 1)
  m = size (X, 2);
  n = size (Y, 2);

  X = center (X);
  Y = center (Y);

  ## Factor with column pivoting and work in the numerical rank.  Without this
  ## a rank-deficient input is solved against a singular triangular factor:
  ## the coefficients run away to 1e15 and the canonical correlations come back
  ## wrong rather than merely imprecise.
  [Qx, Rx, px] = qr (X, 0);
  rankX = rank_of (Rx, k, m);
  if (rankX == 0)
    error ("canoncorr: X must contain at least one non-constant column.");
  elseif (rankX < m)
    warning ("canoncorr:NotFullRank", "canoncorr: X is not full rank.");
    Qx = Qx(:, 1:rankX);
    Rx = Rx(1:rankX, 1:rankX);
  endif

  [Qy, Ry, py] = qr (Y, 0);
  rankY = rank_of (Ry, k, n);
  if (rankY == 0)
    error ("canoncorr: Y must contain at least one non-constant column.");
  elseif (rankY < n)
    warning ("canoncorr:NotFullRank", "canoncorr: Y is not full rank.");
    Qy = Qy(:, 1:rankY);
    Ry = Ry(1:rankY, 1:rankY);
  endif

  d = min (rankX, rankY);

  [U S V] = svd (Qx' * Qy, 'econ');

  A = Rx \ U(:, 1:d);
  B = Ry \ V(:, 1:d);

  ## A, B are scaled to make the covariance matrices of the outputs U, V
  ## identity matrices
  f = sqrt (k-1);
  A .*= f;
  B .*= f;

  ## Put the coefficients back to their full height in the original column
  ## order, the dropped columns contributing nothing.
  A(px, :) = [A; zeros(m - rankX, d)];
  B(py, :) = [B; zeros(n - rankY, d)];

  if (nargout > 2)
    r = max (0, min (diag (S)(1:d), 1))';
  endif
  if (nargout > 3)
    U = X * A;
  endif
  if (nargout > 4)
    V = Y * B;
  endif

  if (nargout > 5)
    ## The degrees of freedom count the dimensions actually fitted, which is
    ## the rank rather than the number of columns supplied.
    Wilks = fliplr (cumprod (fliplr ((1 - r .^ 2))));
    chisq = - (k - 1 - (rankX + rankY + 1)/2) * log (Wilks);
    df1 = (rankX - (1:d) + 1) .* (rankY - (1:d) + 1);
    pChisq = 1 - chi2cdf (chisq, df1);
    s = sqrt ((df1.^2 - 4) ./ ((rankX - (1:d) + 1).^2 + ...
                               (rankY - (1:d) + 1).^2 - 5));
    df2 = (k - 1 - (rankX + rankY + 1)/2) * s - df1/2 + 1;
    ls = Wilks .^ (1 ./ s);
    F = (1 ./ ls  -  1) .* (df2 ./ df1);
    pF = 1 - fcdf (F, df1, df2);
    stats.Wilks = Wilks;
    stats.df1 = df1;
    stats.df2 = df2;
    stats.F = F;
    stats.pF = pF;
    stats.chisq = chisq;
    stats.pChisq = pChisq;
  endif
endfunction

## Numerical rank of a triangular QR factor, on the scale of its leading entry.
function rk = rank_of (R, nrows, ncols)
  if (isempty (R))
    rk = 0;
  else
    rk = sum (abs (diag (R)) > eps (abs (R(1))) * max (nrows, ncols));
  endif
endfunction

%!shared X, Y, A, B, r, U, V, k, Cuv
%! k = 10;
%! X = [1:k; sin(1:k); cos(1:k)]';
%! Y = [tan(1:k); tanh((1:k)/k)]';
%! [A, B, r, U, V, stats] = canoncorr (X, Y);
%! Cuv = (U' * V) / (k - 1);
%!assert_equal (diag (Cuv)', r, 10 * eps);
%!assert_equal (diag (diag (Cuv)), Cuv, 2 * eps);
%!assert_equal (r, [0.99590, 0.26754], 1E-5);
%!assert_equal (U, center (X) * A, 10 * eps);
%!assert_equal (V, center (Y) * B, 10 * eps);
%!assert_equal (cov (U), eye (size (U, 2)), 10 * eps);
%!assert_equal (cov (V), eye (size (V, 2)), 10 * eps);
%! rand ('state', 1); [A, B, r] = canoncorr (rand (5, 10), rand (5, 20));
%! ## Four, not five: centring a five-row matrix leaves rank at most four, so
%! ## there is no fifth canonical correlation to report.  The count used to
%! ## come from the number of rows rather than the rank.
%!assert_equal (r, ones (1, 4), 10*eps);

## Rank-deficient input is reduced to its rank rather than solved against a
## singular factor.  Verified against MATLAB R2024a.
%!test
%! Xr = [X(:,1), X(:,1), X(:,2)];
%! warning ("off", "canoncorr:NotFullRank", "local");
%! [Ar, Br, rr] = canoncorr (Xr, Y);
%! ## the duplicated column contributes nothing
%! assert_equal (Ar(2,:), zeros (1, columns (Ar)));
%! ## and the coefficients stay finite, where they used to reach 1e15
%! assert_equal (all (isfinite (Ar(:))), true);
%! assert_equal (max (abs (Ar(:))) < 1e3, true);
%! ## dropping the duplicate leaves the same fit as never having had it
%! [A2, B2, r2] = canoncorr (X(:,1:2), Y);
%! assert_equal (rr, r2, 1e-12);
%! assert_equal (Ar([1, 3], :), A2, 1e-10);

%!test  # a constant column carries no information and is dropped
%! Xc = [X(:,1:2), ones(rows (X), 1)];
%! warning ("off", "canoncorr:NotFullRank", "local");
%! [Ac, Bc, rc] = canoncorr (Xc, Y);
%! assert_equal (Ac(3,:), zeros (1, columns (Ac)));
%! [A2, B2, r2] = canoncorr (X(:,1:2), Y);
%! assert_equal (rc, r2, 1e-12);

%!test  # the deficiency is reported
%! Xr = [X(:,1), X(:,1), X(:,2)];
%! fail ("canoncorr (Xr, Y)", "warning", "X is not full rank");

%!test  # the returned coefficients satisfy the definition whatever the rank
%! Xr = [X(:,1), X(:,1), X(:,2)];
%! warning ("off", "canoncorr:NotFullRank", "local");
%! [Ar, Br, rr, Ur, Vr] = canoncorr (Xr, Y);
%! kk = rows (Xr);
%! assert_equal ((Ur' * Ur) / (kk - 1), eye (numel (rr)), 1e-10);
%! assert_equal ((Vr' * Vr) / (kk - 1), eye (numel (rr)), 1e-10);
%! assert_equal ((Ur' * Vr) / (kk - 1), diag (rr), 1e-10);

%!error <X must contain at least one non-constant column.> ...
%! canoncorr (ones (10, 2), [tan(1:10); tanh((1:10)/10)]')
