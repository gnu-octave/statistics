## Copyright (C) 2025 Jayant Chauhan <0001jayant@gmail.com>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.



function [b, se, pval, finalmodel, stats, nextstep, history] = ...
         stepwisefit1 (X, y)

  if (nargin != 2)
    error ("stepwisefit1: exactly two input arguments required");
  endif

  if (! ismatrix (X) || ! isvector (y))
    error ("stepwisefit1: invalid input dimensions");
  endif

  y = y(:);

  ## Handle missing values
  wasnan = any (isnan ([X y]), 2);
  Xc = X(!wasnan, :);
  yc = y(!wasnan);

  n = rows (Xc);
  p = columns (Xc);

  ## Stepwise selection (legacy implementation)
  X_use = stepwisefit (yc, Xc);

  ## Final regression on selected predictors
  Xfinal = [ones(n,1), Xc(:, X_use)];
  [B, BINT, R, RINT, regstats] = regress (yc, Xfinal);

  ## Allocate outputs
  b    = zeros (p,1);
  se   = zeros (p,1);
  pval = zeros (p,1);

  df = n - columns (Xfinal);

  ## Included predictors
  b(X_use) = B(2:end);
  se(X_use) = (BINT(2:end,2) - B(2:end)) ./ tinv (0.975, df);
  pval(X_use) = 2 * (1 - tcdf (abs (B(2:end) ./ se(X_use)), df));

  ## Excluded predictors: conditional refit
  excluded = setdiff (1:p, X_use);

  for j = excluded
    Xj = [ones(n,1), Xc(:, [X_use j])];
    [Bj, BjINT] = regress (yc, Xj);

    bj  = Bj(end);
    sej = (BjINT(end,2) - bj) ./ tinv (0.975, n - columns (Xj));

    b(j) = bj;
    se(j) = sej;
    pval(j) = 2 * (1 - tcdf (abs (bj / sej), n - columns (Xj)));
  endfor

  ## Final model indicator
  finalmodel = false (1,p);
  finalmodel(X_use) = true;

  ## Stats structure
  stats = struct ();
  stats.source    = "stepwisefit";
  stats.df0       = numel (X_use);
  stats.dfe       = n - stats.df0 - 1;
  stats.SStotal   = sum ((yc - mean (yc)).^2);
  stats.SSresid   = sum (R.^2);
  stats.rmse      = sqrt (stats.SSresid / stats.dfe);
  stats.intercept = B(1);
  stats.wasnan    = wasnan;

  ## Placeholders for future phases
  nextstep = 0;
  history  = struct ();

endfunction


%!test
%! % S1.2: default full outputs (numeric contract)
%! X = [7 26 6 60;
%!      1 29 15 52;
%!      11 56 8 20;
%!      11 31 8 47;
%!      7 52 6 33;
%!      11 55 9 22;
%!      3 71 17 6;
%!      1 31 22 44;
%!      2 54 18 22;
%!      21 47 4 26;
%!      1 40 23 34;
%!      11 66 9 12;
%!      10 68 8 12];
%! y = [78.5; 74.3; 104.3; 87.6; 95.9; 109.2;
%!      102.7; 72.5; 93.1; 115.9; 83.8; 113.3; 109.4];
%! [b,se,pval,finalmodel,stats] = stepwisefit1 (X,y);
%! assert (finalmodel, [true false false true]);
%! assert (b, [1.4400; 0.4161; -0.4100; -0.6140], 1e-4);
%! assert (se, [0.1384; 0.1856; 0.1992; 0.0486], 1e-4);
%! assert (pval, [0; 0.0517; 0.0697; 0], 1e-4);
%! assert (stats.rmse, 2.7343, 1e-4);
%! assert (stats.SStotal, 2715.7631, 1e-3);
%! assert (stats.SSresid, 74.7621, 1e-4);
%! assert (stats.df0, 2);
%! assert (stats.dfe, 10);
%! assert (stats.intercept, 103.0974, 1e-4);

