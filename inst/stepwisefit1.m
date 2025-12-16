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

# -*- texinfo -*-
## @deftypefn {statistics} {@var{b}, @var{se}, @var{pval}, @var{finalmodel}, @var{stats}, @var{nextstep}, @var{history}} = stepwisefit1 (@var{X}, @var{y})
##
## Stepwise linear regression with expanded diagnostic outputs.
##
## This function performs stepwise variable selection using linear regression
## and returns coefficient estimates, standard errors, p-values, and model
## diagnostics in a consolidated form.
##
## Predictor selection is performed using an internal stepwise procedure.
## After selection, the final regression model is fit explicitly to compute
## inferential statistics for both included and excluded predictors.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{X} is an @var{n}-by-@var{p} matrix of predictor variables.
## @item
## @var{y} is an @var{n}-by-1 response vector.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{b} is a @var{p}-by-1 vector of regression coefficients.
## @item
## @var{se} is a @var{p}-by-1 vector of standard errors.
## @item
## @var{pval} is a @var{p}-by-1 vector of two-sided p-values.
## @item
## @var{finalmodel} is a logical row vector indicating which predictors are
## included in the final model.
## @item
## @var{stats} is a structure containing regression diagnostics, including
## degrees of freedom, sums of squares, root mean squared error, and the
## intercept term.
## @item
## @var{nextstep} is a scalar indicating whether an additional step is
## recommended.
## @item
## @var{history} is a structure recording intermediate model states during
## the stepwise procedure.
## @end itemize
##
## @seealso{stepwisefit, regress}
## @end deftypefn

function [b, se, pval, finalmodel, stats, nextstep, history] = ...
         stepwisefit1 (X, y, varargin)

  ## Input validation (positional)

  if (nargin < 2)
    error ("stepwisefit1: at least two input arguments required");
  endif

  if (! ismatrix (X) || ! isvector (y))
    error ("stepwisefit1: X must be a matrix and y a vector");
  endif

  y = y(:);

  ## Parse Name–Value pairs
  InModel  = [];
  Display  = "on";

  if (mod (numel (varargin), 2) != 0)
    error ("stepwisefit1: Name–Value arguments must come in pairs");
  endif

  for k = 1:2:numel (varargin)
    name  = varargin{k};
    value = varargin{k+1};

    if (! (ischar (name) || isstring (name)))
      error ("stepwisefit1: Name–Value keys must be strings");
    endif
    name = char (name);

    switch lower (name)
      case "inmodel"
        if (! islogical (value))
          error ("stepwisefit1: InModel must be a logical vector");
        endif
        InModel = value(:).';
      case "display"
        if (! any (strcmpi (value, {"on", "off"})))
          error ("stepwisefit1: Display must be 'on' or 'off'");
        endif
        Display = lower (value);
      otherwise
        error ("stepwisefit1: Name–Value option '%s' not supported", name);
    endswitch
  endfor

  ## Handle missing values
  wasnan = any (isnan ([X y]), 2);
  Xc = X(!wasnan, :);
  yc = y(!wasnan);

  n = rows (Xc);
  p = columns (Xc);

  %% ----------------------------
  %% Validate InModel
  %% ----------------------------
  if (! isempty (InModel))
    if (numel (InModel) != p)
      error ("stepwisefit1: InModel length must match number of predictors");
    endif
  endif

  ## Stepwise variable selection
  if (isempty (InModel))
  X_use = stepwisefit (yc, Xc);
  else
    % Force initial model, then allow expansion
    X_init = find (InModel);
    X_use  = X_init;

    % Run legacy stepwisefit to allow expansion
    X_more = stepwisefit (yc, Xc(:, !InModel));
    X_use  = unique ([X_init, find (!InModel)(X_more)]);
  endif


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

  stats.yr = R;
  stats.B     = b;
  stats.SE    = se;
  stats.TSTAT = b ./ se;
  stats.PVAL  = pval;
  stats.TSTAT (!isfinite (stats.TSTAT)) = NaN;
  excluded = setdiff (1:p, X_use);
  xr = zeros (n, numel (excluded));

  for k = 1:numel (excluded)
    j = excluded(k);
    Xproj = [ones(n,1), Xc(:, X_use)];
    bj = regress (Xc(:,j), Xproj);
    xr(:,k) = Xc(:,j) - Xproj * bj;
  endfor

  stats.xr = xr;
  covB = (stats.rmse^2) * inv (Xfinal' * Xfinal);
  covb = NaN (p+1, p+1);
  covb(1:rows (covB), 1:columns (covB)) = covB;
  stats.covb = covb;

  stats.fstat = ((stats.SStotal - stats.SSresid) / stats.df0) ...
                / (stats.SSresid / stats.dfe);

  stats.pval = 1 - fcdf (stats.fstat, stats.df0, stats.dfe);

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

%!test
%! % S3.1 (numeric kernel only): forced baseline selection
%! X = [
%!   12.0 4 120 95 2600;
%!   11.5 6 200 110 3000;
%!   10.5 8 300 150 3600;
%!   13.0 4 140 100 2800;
%!   12.5 6 180 120 3200;
%!   11.0 8 250 140 3500;
%!   14.0 4 130 98 2700;
%!   13.5 6 210 115 3100;
%!   12.2 8 320 160 3800;
%!   11.8 4 150 105 2900
%! ];
%! y = [28; 22; 18; 27; 23; 19; 29; 21; 17; 26];
%!
%! [b,se,pval,finalmodel,stats] = stepwisefit1 (X,y);
%!
%! assert (islogical (finalmodel));
%! assert (numel (finalmodel) == 5);
%! assert (sum (finalmodel) >= 1);
%! assert (isnumeric (b));
%! assert (isnumeric (se));
%! assert (isnumeric (pval));
%! assert (stats.rmse > 0);
%! assert (isfinite (stats.intercept));

%!test
%! X = randn (30, 4);
%! y = randn (30, 1);
%! [~,~,~,~,stats] = stepwisefit1 (X, y);
%!
%! required_fields = {
%!   "source", "df0", "dfe", "SStotal", "SSresid", "fstat", "pval", ...
%!   "rmse", "xr", "yr", "B", "SE", "TSTAT", "PVAL", "covb", ...
%!   "intercept", "wasnan"
%! };
%!
%! for k = 1:numel (required_fields)
%!   assert (isfield (stats, required_fields{k}));
%! endfor

%!test
%! X = randn (40, 5);
%! y = randn (40, 1);
%! [b,se,pval,finalmodel,stats] = stepwisefit1 (X, y);
%!
%! p = columns (X);
%! n = rows (X(~stats.wasnan, :));
%!
%! assert (size (stats.yr), [n, 1]);
%! assert (rows (stats.B) == p);
%! assert (rows (stats.SE) == p);
%! assert (rows (stats.TSTAT) == p);
%! assert (rows (stats.PVAL) == p);
%! assert (size (stats.covb), [p+1, p+1]);

%!test
%! X = randn (25, 3);
%! y = randn (25, 1);
%! [~,~,~,~,stats] = stepwisefit1 (X, y);
%!
%! SSresid_calc = sum (stats.yr .^ 2);
%! assert (SSresid_calc, stats.SSresid, 1e-10);
%!
%! rmse_calc = sqrt (stats.SSresid / stats.dfe);
%! assert (rmse_calc, stats.rmse, 1e-10);

%!test
%! X = randn (50, 6);
%! y = randn (50, 1);
%! [~,~,~,~,stats] = stepwisefit1 (X, y);
%!
%! if (stats.df0 > 0)
%!   F_calc = ((stats.SStotal - stats.SSresid) / stats.df0) ...
%!            / (stats.SSresid / stats.dfe);
%!
%!   assert (F_calc, stats.fstat, 1e-10);
%!   assert (stats.pval >= 0 && stats.pval <= 1);
%! else
%!   assert (isnan (stats.fstat));
%!   assert (isnan (stats.pval));
%! endif


%!test
%! X = randn (35, 4);
%! y = randn (35, 1);
%! [~,~,~,finalmodel,stats] = stepwisefit1 (X, y);
%!
%! Xc = X(~stats.wasnan, :);
%! Xfinal = [ones(rows (Xc),1), Xc(:, finalmodel)];
%!
%! for j = 1:columns (stats.xr)
%!   corrval = corr (stats.xr(:,j), Xfinal(:,2:end));
%!   assert (max (abs (corrval(:))) < 1e-10);
%! endfor

%!test
%! X = randn (35, 4);
%! y = randn (35, 1);
%! [~,~,~,finalmodel,stats] = stepwisefit1 (X, y);
%!
%! Xc = X(~stats.wasnan, :);
%! Xfinal = [ones(rows (Xc),1), Xc(:, finalmodel)];
%!
%! for j = 1:columns (stats.xr)
%!   ortho = Xfinal' * stats.xr(:,j);
%!   assert (max (abs (ortho(:))) < 1e-8);
%! endfor
