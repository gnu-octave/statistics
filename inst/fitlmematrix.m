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
## @deftypefn  {statistics} {@var{lme} =} fitlmematrix (@var{X}, @var{y}, @var{Z}, @var{G})
## @deftypefnx {statistics} {@var{lme} =} fitlmematrix (@dots{}, @var{name}, @var{value})
##
## Fit a linear mixed-effects model from design matrices.
##
## @code{fitlmematrix (@var{X}, @var{y}, @var{Z}, @var{G})} fits the linear
## mixed-effects model
## @tex
## $y = X\beta + Zb + \varepsilon$
## @end tex
## @ifnottex
## @code{y = X*beta + Z*b + e}
## @end ifnottex
## with fixed-effects design @var{X}, response @var{y}, random-effects design
## @var{Z}, and grouping variable @var{G}.  The random effects @var{b} are
## normally distributed with mean zero and an unstructured covariance @var{Psi}
## (shared across the levels of the grouping variable), and the observation
## errors are independent @code{N(0, sigma2)}.
##
## @var{X} is an @var{n}-by-@var{p} numeric matrix and @var{y} an @var{n}-by-1
## response vector.  @var{Z} is an @var{n}-by-@var{q} random-effects design and
## @var{G} an @var{n}-by-1 grouping variable (numeric, logical, char, cell array
## of strings, or categorical).  To specify several grouping terms, pass @var{Z}
## and @var{G} as cell arrays of the same length, one design and one grouping
## variable per term.
##
## The following @var{name}/@var{value} pairs are accepted:
##
## @table @asis
## @item @qcode{"FitMethod"}
## The estimation criterion, either @qcode{"ML"} (maximum likelihood, the
## default) or @qcode{"REML"} (restricted maximum likelihood).
##
## @item @qcode{"FixedEffectPredictors"}
## A cell array of @var{p} names for the columns of @var{X} (default
## @code{@{"x1", @dots{}, "xp"@}}).
##
## @item @qcode{"RandomEffectPredictors"}
## A cell array (one entry per grouping term) of cell arrays naming the columns
## of each @var{Z} (default @code{z1, z2, @dots{}}).
##
## @item @qcode{"RandomEffectGroups"}
## A cell array of names for the grouping terms (default @code{g1, g2}, etc.).
## @end table
##
## The returned @var{lme} is a @code{LinearMixedModel} object describing the
## fitted model: the estimated fixed effects and their statistics
## (@code{lme.Coefficients}), the covariance parameters
## (@code{covarianceParameters}), the random-effects BLUPs (@code{randomEffects}),
## the log-likelihood, and methods for prediction, residuals, and hypothesis
## tests.
##
## Only the full (unstructured) random-effects covariance is currently
## supported.
##
## @seealso{LinearMixedModel, fitlm, parseWilkinsonFormula}
## @end deftypefn

function lme = fitlmematrix (X, y, Z, G, varargin)

  if (nargin < 4)
    print_usage ();
  endif

  ## --- fixed-effects design and response ---
  if (! (isnumeric (X) && ismatrix (X) && isreal (X)))
    error ("fitlmematrix: X must be a real numeric matrix.");
  endif
  if (! (isnumeric (y) && isreal (y) && isvector (y)))
    error ("fitlmematrix: y must be a real numeric vector.");
  endif
  y = y(:);
  n = rows (X);
  if (numel (y) != n)
    error ("fitlmematrix: y must have as many elements as X has rows.");
  endif

  ## --- random-effects designs and grouping (normalise to cells) ---
  if (! iscell (Z))
    Z = {Z};
  endif
  if (! iscell (G))
    G = {G};
  endif
  if (numel (Z) != numel (G))
    error ("fitlmematrix: Z and G must have the same number of terms.");
  endif
  nt = numel (Z);
  for k = 1:nt
    if (! (isnumeric (Z{k}) && ismatrix (Z{k}) && isreal (Z{k})))
      error ("fitlmematrix: each Z must be a real numeric matrix.");
    endif
    if (rows (Z{k}) != n)
      error ("fitlmematrix: each Z must have as many rows as X.");
    endif
    G{k} = G{k}(:);
    if (numel (G{k}) != n)
      error ("fitlmematrix: each G must have as many elements as X has rows.");
    endif
  endfor

  ## --- name/value options ---
  method = "ML";
  fenames = {};
  renames = {};
  grnames = {};
  formula_str = "";
  respname = "";
  if (mod (numel (varargin), 2) != 0)
    error ("fitlmematrix: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    name = varargin{i};
    value = varargin{i+1};
    if (! (ischar (name) || (isstring (name) && isscalar (name))))
      error ("fitlmematrix: option names must be strings.");
    endif
    switch (lower (char (name)))
      case "fitmethod"
        method = upper (char (value));
        if (! any (strcmp (method, {"ML", "REML"})))
          error ("fitlmematrix: FitMethod must be 'ML' or 'REML'.");
        endif
      case "fixedeffectpredictors"
        fenames = cellstr (value);
      case "randomeffectpredictors"
        renames = value;
      case "randomeffectgroups"
        grnames = cellstr (value);
      case "formula"
        formula_str = char (value);
      case "responsename"
        respname = char (value);
      case "covariancepattern"
        pat = lower (char (value));
        if (! any (strcmp (pat, {"fullcholesky", "full"})))
          error (strcat ("fitlmematrix: only the 'FullCholesky' covariance", ...
                         " pattern is currently supported."));
        endif
      otherwise
        error ("fitlmematrix: unknown option '%s'.", char (name));
    endswitch
  endfor

  p = columns (X);
  if (isempty (fenames))
    fenames = arrayfun (@(j) sprintf ("x%d", j), 1:p, "UniformOutput", false);
  elseif (numel (fenames) != p)
    error ("fitlmematrix: FixedEffectPredictors must name every column of X.");
  endif
  if (isempty (grnames))
    grnames = arrayfun (@(k) sprintf ("g%d", k), 1:nt, "UniformOutput", false);
  endif
  if (isempty (renames))
    renames = cell (1, nt);
    for k = 1:nt
      renames{k} = arrayfun (@(j) sprintf ("z%d", j), 1:columns (Z{k}), ...
                             "UniformOutput", false);
    endfor
  endif

  ## --- fit ---
  fit = __lmefit__ (X, y, Z, G, method);

  ## --- build the info struct and wrap it in a LinearMixedModel object ---
  info = fit;
  info.X = X;             info.y = y;
  info.Zcell = Z;         info.Gcell = G;
  info.CoefficientNames = fenames;
  info.GroupNames = grnames;
  info.REPred = renames;
  info.method = method;
  info.Formula = formula_str;
  info.ResponseName = respname;

  lme = LinearMixedModel (info);

endfunction

%!demo
%! ## A random-intercept model on a small balanced data set: five subjects
%! ## measured at four values of a predictor x.
%! x = repmat ([1 2 3 4]', 5, 1);
%! subject = reshape (repmat (1:5, 4, 1), [], 1);
%! y = 2 + 0.8 * x + reshape (repmat ([1 -1 0.5 -0.5 0]', 1, 4)', [], 1) ...
%!         + 0.1 * sin (1:20)';
%! X = [ones(20,1), x];
%! lme = fitlmematrix (X, y, ones (20, 1), subject, "FitMethod", "REML");
%! beta = fixedEffects (lme);
%! [psi, mse] = covarianceParameters (lme);
%! printf ("intercept = %.4f, slope = %.4f\n", beta);
%! printf ("between-subject var = %.4f, residual var = %.4f\n", psi{1}, mse);

## --- Balanced one-way random-intercept model: closed-form REML/ANOVA check ---
## For a balanced one-way design (a groups of r observations) the REML estimates
## coincide with the ANOVA estimates: sigma2 = MSE, tau2 = (MSB - MSE)/r, and
## the fixed intercept is the grand mean.
%!test
%! a = 5; r = 6; n = a*r;
%! grp = reshape (repmat (1:a, r, 1), [], 1);
%! yv = [ 4.1 4.5 3.8 4.3 4.0 4.2, 6.0 5.7 6.3 5.9 6.1 6.2, ...
%!        2.2 2.5 2.0 2.4 2.1 2.3, 5.1 4.9 5.3 5.0 5.2 4.8, ...
%!        3.3 3.6 3.1 3.4 3.2 3.5 ]';
%! X = ones (n, 1);
%! lme = fitlmematrix (X, yv, ones (n, 1), grp, "FitMethod", "REML");
%! gm = mean (yv);
%! gmeans = accumarray (grp, yv, [], @mean);
%! SSB = r * sum ((gmeans - gm) .^ 2);   MSB = SSB / (a - 1);
%! SSE = sum ((yv - gmeans(grp)) .^ 2);  MSE = SSE / (n - a);
%! tau2 = (MSB - MSE) / r;
%! [psi, mse] = covarianceParameters (lme);
%! assert (fixedEffects (lme), gm, 1e-8);
%! assert (mse, MSE, 1e-6);
%! assert (psi{1}, tau2, 1e-6);

## --- Internal consistency: sigma2 and covbeta satisfy their defining GLS
## relations at the returned fit. ---
%!test
%! a = 4; r = 5; n = a*r;
%! grp = reshape (repmat (1:a, r, 1), [], 1);
%! x = (1:n)' / n;
%! yv = 1 + 2*x + reshape (repmat ([0.5 -0.5 0.2 -0.2], r, 1), [], 1) ...
%!        + 0.05 * cos (1:n)';
%! X = [ones(n,1), x];
%! lme = fitlmematrix (X, yv, ones (n, 1), grp, "FitMethod", "REML");
%! ## rebuild V = Psi over groups + sigma2 I and check the GLS identities
%! [psi, mse] = covarianceParameters (lme);
%! Zx = zeros (n, a);
%! for l = 1:a, Zx(grp==l, l) = 1; end
%! V = psi{1} * (Zx * Zx') + mse * eye (n);
%! Vi = inv (V);
%! beta = (X' * Vi * X) \ (X' * Vi * yv);
%! assert (fixedEffects (lme), beta, 1e-6);
%! assert (lme.CoefficientCovariance, inv (X' * Vi * X), 1e-5);

## --- MATLAB-verified parity: unbalanced random-intercept fit (ML). ---
## Reference values are MATLAB fitlmematrix (R2026a) on this exact data.
%!shared X, yL, grp, xL, x2
%! xL = [0.032760004 0.70410822 -0.8646718 -0.28869454 0.51276678 -1.4975462 ...
%!  -1.4527871 -0.80013541 -1.644209 1.5137701 0.72905543 0.20880758 1.0856145 ...
%!  0.62862577 -0.87409978 1.9178276 0.09748204 0.50697633 1.0247569 ...
%!  -0.92789896 -0.88921018 -0.98322849 -0.031378913 0.86875961 -0.91481141 ...
%!  0.034324163 -0.25025257 -1.0575644 -0.86131607 -0.35355444 0.82950729 ...
%!  -0.36874363 0.061580868 0.55803564 -0.1763803 1.0482413 1.0137831 ...
%!  -0.94876976 -0.010703972 -0.35149845 -1.6828735 -1.0493301]';
%! x2 = [0.68979276 0.0074354814 -0.45697437 -0.5636481 1.4567202 -0.97829955 ...
%!  -1.12922 -0.030542479 1.5847779 -0.87837755 0.24121762 0.68747601 ...
%!  -0.56728765 0.98895053 -0.39350661 0.85326015 0.36524343 0.15824977 ...
%!  -1.7665212 0.59808246 -0.55763708 -1.1982294 -2.1473319 0.22521416 ...
%!  0.37034398 -1.880586 0.052941033 -0.70016994 0.2174853 -1.7797082 ...
%!  0.51971317 -0.35551286 1.9845963 -1.3498848 -0.63514097 -0.78794714 ...
%!  1.3681179 1.4423152 -0.51233905 0.30238864 2.0458136 0.17326323]';
%! yL = [3.5635971 -0.36763498 1.4192715 3.4471073 2.2760668 3.6420287 ...
%!  4.4481522 1.5292777 5.437158 0.44016873 0.7259756 2.9771749 1.2941347 ...
%!  0.26407508 1.9673466 0.77589984 1.2176078 0.8528384 -0.86522876 2.8651388 ...
%!  2.0127931 3.1242841 -0.49164361 1.3192295 5.3782287 -0.40680701 2.3989882 ...
%!  3.606108 1.9193197 1.6713765 2.4383124 1.9314844 3.5112706 0.91644553 ...
%!  0.034923706 -0.15510033 2.1765206 2.5327412 3.1421219 3.3472574 5.343425 ...
%!  3.8775899]';
%! grp = [1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 ...
%!  5 6 1 2 3 4 5 6]';
%! X = [ones(42,1), xL, x2];

%!test  # returns a LinearMixedModel object
%! lme = fitlmematrix (X, yL, ones (42, 1), grp, "FitMethod", "ML");
%! assert (isa (lme, "LinearMixedModel"));

%!test  # random intercept, ML -- matches MATLAB fitlmematrix
%! lme = fitlmematrix (X, yL, ones (42, 1), grp, "FitMethod", "ML");
%! [psi, mse] = covarianceParameters (lme);
%! assert (fixedEffects (lme), [1.9685543; -1.3771207; 0.81016147], 1e-4);
%! assert (mse, 0.36422151, 1e-4);
%! assert (psi{1}, 0.65265828, 1e-4);
%! assert (lme.LogLikelihood, -46.203282, 1e-4);

%!test  # correlated random intercept + slope, REML -- matches MATLAB
%! Z = [ones(42,1), xL];
%! lme = fitlmematrix (X, yL, Z, grp, "FitMethod", "REML");
%! [psi, mse] = covarianceParameters (lme);
%! assert (fixedEffects (lme), [2.0034354; -1.3374715; 0.83788802], 1e-3);
%! assert (mse, 0.36649172, 1e-3);
%! assert (lme.LogLikelihood, -48.280962, 1e-3);
%! assert (psi{1}, [0.7846112, -0.14287556; -0.14287556, 0.026017248], 1e-3);

%!test  # ML and REML give different (both sensible) variance components
%! lme_ml = fitlmematrix (X, yL, ones (42, 1), grp, "FitMethod", "ML");
%! lme_re = fitlmematrix (X, yL, ones (42, 1), grp, "FitMethod", "REML");
%! pml = covarianceParameters (lme_ml);
%! pre = covarianceParameters (lme_re);
%! assert (pre{1} > pml{1});
%! assert (fixedEffects (lme_ml), fixedEffects (lme_re), 5e-3);  # close, not equal

%!test  # default FitMethod is ML
%! l1 = fitlmematrix (X, yL, ones (42, 1), grp);
%! l2 = fitlmematrix (X, yL, ones (42, 1), grp, "FitMethod", "ML");
%! assert (l1.LogLikelihood, l2.LogLikelihood, 1e-10);

%!test  # names propagate to the fitted object
%! lme = fitlmematrix (X, yL, ones (42, 1), grp, ...
%!                     "FixedEffectPredictors", {"Int", "x", "x2"});
%! assert (lme.CoefficientNames, {"Int", "x", "x2"});

## Input validation
%!error <Invalid call> fitlmematrix (1, 2)
%!error <X must be a real numeric matrix> fitlmematrix ({1}, 1, 1, 1)
%!error <y must be a real numeric vector> fitlmematrix (ones (3), {1}, 1, 1)
%!error <y must have as many elements> fitlmematrix (ones (3, 2), [1 2], ones (3, 1), [1 2 3])
%!error <Z and G must have the same number> fitlmematrix (ones (3, 2), [1;2;3], {ones(3,1), ones(3,1)}, {[1;2;3]})
%!error <each Z must have as many rows> fitlmematrix (ones (3, 2), [1;2;3], ones (2, 1), [1;2;3])
%!error <FitMethod must be 'ML' or 'REML'> fitlmematrix (ones (3, 2), [1;2;3], ones (3, 1), [1;2;3], "FitMethod", "xxx")
%!error <name/value arguments must come in pairs> fitlmematrix (ones (3, 2), [1;2;3], ones (3, 1), [1;2;3], "FitMethod")
%!error <unknown option 'bogus'> fitlmematrix (ones (3, 2), [1;2;3], ones (3, 1), [1;2;3], "bogus", 1)
%!error <only the 'FullCholesky' covariance> fitlmematrix (ones (3, 2), [1;2;3], ones (3, 1), [1;2;3], "CovariancePattern", "Diagonal")
