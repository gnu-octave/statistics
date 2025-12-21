## Copyright (C) 2013-2021 Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Copyright (C) 2014 Mikael Kurula <mkurula@abo.fi>
## Copyright (C) 2025 Jayant Chauhan <0001jayant@gmail.com>
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
## @deftypefn {statistics} {@var{b}, @var{se}, @var{pval}, @var{finalmodel}, @var{stats}, @var{nextstep}, @var{history}} = stepwisefit (@var{X}, @var{y}, @var{varargin})
##
## Perform stepwise linear regression with extended diagnostic outputs.
##
## @code{stepwisefit} provides additional regression diagnostics, structured
## outputs, and partial MATLAB-compatible Name–Value argument handling.
##
## Predictor selection is performed using the existing @code{stepwisefit}
## algorithm. After variable selection, the final regression model is
## explicitly refit using @code{regress} to compute coefficient estimates
## and inferential statistics for both included and excluded predictors.
##
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{X} is an @var{n}-by-@var{p} numeric matrix of predictor variables.
##
## @item
## @var{y} is an @var{n}-by-1 numeric response vector.
##
## @item
## Optional Name–Value pairs may be provided to control model selection.
## @end itemize
##
## @subheading Name–Value Arguments
##
## @table @asis
## @item @qcode{"InModel"}
## Logical row vector of length @var{p} specifying an initial model.
## Predictors marked @code{true} are treated as initially included.
##
## @item @qcode{"Keep"}
## Logical row vector of length @var{p} specifying predictors that must
## always be included in the final model.
##
## @item @qcode{"PEnter"}
## Scalar significance level in the open interval (0,1) specifying the
## entry threshold for stepwise selection.  Default is @code{0.05}.
##
## @item @qcode{"PRemove"}
## Scalar significance level in the open interval (0,1) specifying the
## removal threshold for stepwise selection.  If not specified, a default
## value greater than or equal to @qcode{"PEnter"} is used.
##
## @item @qcode{"MaxIter"}
## Positive integer specifying the maximum number of stepwise iterations.
## Default is @code{Inf}.
##
## @item @qcode{"Scale"}
## Either @qcode{"on"} or @qcode{"off"}.  When enabled, predictors are
## standardized prior to stepwise selection only.  Final coefficients
## are always reported on the original data scale.
##
## @item @qcode{"Display"}
## Either @qcode{"on"} or @qcode{"off"}.  Currently accepted for interface
## compatibility but does not affect output.
## @end table
##
## @subheading Return Values
##
## @itemize @bullet
## @item
## @var{b} is a @var{p}-by-1 vector of regression coefficients.  Coefficients
## corresponding to excluded predictors are estimated conditionally.
##
## @item
## @var{se} is a @var{p}-by-1 vector of standard errors.
##
## @item
## @var{pval} is a @var{p}-by-1 vector of two-sided p-values.
##
## @item
## @var{finalmodel} is a logical row vector indicating predictors selected
## in the final model.
##
## @item
## @var{stats} is a structure containing regression diagnostics, including
## sums of squares, degrees of freedom, residuals, covariance estimates,
## F-statistic, and related quantities.
##
## @item
## @var{nextstep} is a scalar placeholder indicating whether an additional
## stepwise iteration is recommended.  Currently always zero.
##
## @item
## @var{history} is a structure summarizing the final model state, including
## selected predictors and coefficient history.
## @end itemize
##
## @seealso{regress}
## @end deftypefn

function [b, se, pval, finalmodel, stats, nextstep, history] = ...
         stepwisefit (X, y, varargin)

  ## Input validation (positional)

  if (nargin < 2)
    error ("stepwisefit: at least two input arguments required");
  endif

  if (! ismatrix (X) || ! isvector (y))
    error ("stepwisefit: X must be a matrix and y a vector");
  endif

  y = y(:);

  ## Parse Name–Value pairs
  InModel  = [];
  Display  = "on";
  % MATLAB-compatible defaults
  PEnter  = 0.05;
  PRemove = [];
  Scale   = "off";
  MaxIter = Inf;
  Keep    = [];


  if (mod (numel (varargin), 2) != 0)
    error ("stepwisefit: Name–Value arguments must come in pairs");
  endif

  for k = 1:2:numel (varargin)
    name  = varargin{k};
    value = varargin{k+1};

    if (! (ischar (name) || isstring (name)))
      error ("stepwisefit: Name–Value keys must be strings");
    endif
    name = char (name);

    switch lower (name)
      case "inmodel"
        if (! islogical (value))
          error ("stepwisefit: InModel must be a logical vector");
        endif
        InModel = value(:).';
      case "display"
        if (! any (strcmpi (value, {"on", "off"})))
          error ("stepwisefit: Display must be 'on' or 'off'");
        endif
        Display = lower (value);

      case "penter"
        if (! isscalar (value) || ! isnumeric (value) || value <= 0 || value >= 1)
          error ("stepwisefit: PEnter must be a scalar strictly between 0 and 1");
        endif
        PEnter = value;

      case "premove"
        if (! isscalar (value) || ! isnumeric (value) || value <= 0 || value >= 1)
          error ("stepwisefit: PRemove must be a scalar strictly between 0 and 1");
        endif
        PRemove = value;

      case "scale"
        if (! any (strcmpi (value, {"on", "off"})))
          error ("stepwisefit: Scale must be 'on' or 'off'");
        endif
        Scale = lower (value);

      case "maxiter"
        if (! isscalar (value) || value <= 0 || fix (value) != value)
          error ("stepwisefit: MaxIter must be a positive integer");
        endif
        MaxIter = value;

      case "keep"
        if (! islogical (value))
          error ("stepwisefit: Keep must be a logical vector");
        endif
        Keep = value(:).';

      otherwise
        error ("stepwisefit: Name–Value option '%s' not supported", name);

    endswitch
  endfor

  ## Handle missing values
  wasnan = any (isnan ([X y]), 2);
  Xc = X(!wasnan, :);
  yc = y(!wasnan);

  n = rows (Xc);
  p = columns (Xc);

  if (! isempty (Keep))
    X_forced = find (Keep);
  else
    X_forced = [];
  endif

  free_idx = setdiff (1:p, X_forced);

  if (strcmp (Scale, "on"))
    muX = mean (Xc, 1);
    sigX = std (Xc, 0, 1);
    sigX(sigX == 0) = 1;   % prevent division by zero
    Xs = (Xc - muX) ./ sigX;
  else
    Xs = Xc;
  endif

  if (isempty (free_idx))
    X_use = X_forced;
  else
    X_step_prev = [];
    X_step = [];

    iter = 0;
    while (iter < MaxIter)
      iter++;
      X_step = stepwisefit (yc, Xs(:, free_idx), PEnter, PRemove);

      if (isequal (X_step, X_step_prev))
        break;
      endif

      X_step_prev = X_step;
    endwhile

    X_use = sort ([X_forced, free_idx(X_step)]);
  endif

  ## Validate InModel
  
  if (! isempty (InModel))
    if (numel (InModel) != p)
      error ("stepwisefit: InModel length must match number of predictors");
    endif
  endif

    if (isempty (PRemove))
    PRemove = max (PEnter, 0.1);
  endif

  if (PRemove < PEnter)
    error ("stepwisefit: PRemove must be greater than or equal to PEnter");
  endif


  ## Final regression on selected predictors
  Xfinal = [ones(n,1), Xc(:, X_use)];
  [B, BINT, R, RINT, regstats] = regress (yc, Xfinal);

  Rresid = R(:);   % freeze residual vector

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
  stats.SSresid   = sum (Rresid.^2);
  stats.rmse      = sqrt (stats.SSresid / stats.dfe);
  stats.intercept = B(1);
  stats.wasnan    = wasnan;

  stats.yr = Rresid;
  stats.B     = b;
  stats.SE    = se;
  stats.TSTAT = b ./ se;
  stats.PVAL  = pval;
  stats.TSTAT (!isfinite (stats.TSTAT)) = NaN;
  
  excluded = setdiff (1:p, X_use);
xr = zeros (n, numel (excluded));

if (! isempty (X_use))
  Z = [ones(n,1), Xc(:, X_use)];
  P = Z / (Z' * Z) * Z';   % projection matrix
  for k = 1:numel (excluded)
    j = excluded(k);
    xr(:,k) = Xc(:,j) - P * Xc(:,j);
  endfor
else
  % intercept-only case
  for k = 1:numel (excluded)
    j = excluded(k);
    xr(:,k) = Xc(:,j) - mean (Xc(:,j));
  endfor
endif

stats.xr = xr;
  
  covb = NaN (p+1, p+1);
  covB = (stats.rmse^2) * inv (Xfinal' * Xfinal);

  idx = [1, X_use + 1];
  covb(idx, idx) = covB;

  stats.covb = covb;

  stats.fstat = ((stats.SStotal - stats.SSresid) / stats.df0) ...
                / (stats.SSresid / stats.dfe);

  stats.pval = 1 - fcdf (stats.fstat, stats.df0, stats.dfe);

  history = struct ();
  history.in = finalmodel;
  history.df0 = stats.df0;
  history.rmse = stats.rmse;

  % Coefficient history (excluding intercept)
  % MATLAB stores this as p-by-k; here k = 1
  Bhist = zeros (p, 1);
  Bhist(finalmodel) = b(finalmodel);
  history.B = Bhist;
 
  ## Placeholders for future phases
  nextstep = 0;

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
%! [b,se,pval,finalmodel,stats] = stepwisefit (X,y);
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
%! [b,se,pval,finalmodel,stats] = stepwisefit (X,y);
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
%! [~,~,~,~,stats] = stepwisefit (X, y);
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
%! [b,se,pval,finalmodel,stats] = stepwisefit (X, y);
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
%! [~,~,~,~,stats] = stepwisefit (X, y);
%!
%! SSresid_calc = sum (stats.yr .^ 2);
%! assert (SSresid_calc, stats.SSresid, 1e-10);
%!
%! rmse_calc = sqrt (stats.SSresid / stats.dfe);
%! assert (rmse_calc, stats.rmse, 1e-10);

%!test
%! X = randn (50, 6);
%! y = randn (50, 1);
%! [~,~,~,~,stats] = stepwisefit (X, y);
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
%! [~,~,~,finalmodel,stats] = stepwisefit (X, y);
%! p = columns (X);
%! k = sum (finalmodel);
%! assert (size (stats.xr, 2) == p - k);
%! assert (all (isfinite (stats.xr(:))));

%!test
%! X = randn (35, 4);
%! y = randn (35, 1);
%! [~,~,~,finalmodel,stats] = stepwisefit (X, y);
%!
%! Xc = X(~stats.wasnan, :);
%! Xfinal = [ones(rows (Xc),1), Xc(:, finalmodel)];
%!
%! for j = 1:columns (stats.xr)
%!   ortho = Xfinal' * stats.xr(:,j);
%!   assert (max (abs (ortho(:))) < 1e-6);
%! endfor

%!test
%! X = randn (40, 5);
%! y = randn (40, 1);
%! [~,~,~,finalmodel,stats,nextstep,history] = stepwisefit (X, y);
%!
%! assert (nextstep == 0);
%! assert (isstruct (history));
%! assert (isfield (history, "in"));
%! assert (isfield (history, "df0"));
%! assert (isfield (history, "rmse"));
%! assert (isfield (history, "B"));
%!
%! assert (isequal (history.in, finalmodel));
%! assert (history.df0 == stats.df0);
%! assert (history.rmse == stats.rmse);
%! assert (rows (history.B) == columns (X));

%!test
%! X = randn (20,4);
%! y = randn (20,1);
%! stepwisefit (X,y,'Keep',[true false true false]);

%!test
%! X = randn (20,4);
%! y = randn (20,1);
%! try
%!   stepwisefit (X,y,'Keep',[true false]);
%!   error ("Expected error not thrown");
%! catch
%!   assert (true);
%! end_try_catch

%!test
%! X = randn (30, 4);
%! y = randn (30, 1);
%! keep = [true false false false];
%! [~,~,~,finalmodel] = stepwisefit (X, y, "Keep", keep);
%! assert (finalmodel(1) == true);

%!test
%! X = randn (40, 6);
%! y = randn (40, 1);
%! [~,~,~,finalmodel] = stepwisefit (X, y, "MaxIter", 1);
%! assert (islogical (finalmodel));

%!test
%! X = randn (50, 5);
%! y = randn (50, 1);
%! [b1] = stepwisefit (X, y);
%! [b2] = stepwisefit (X, y, "Scale", "on");
%! assert (rows (b1) == rows (b2));
