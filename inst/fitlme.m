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
## @deftypefn  {statistics} {@var{lme} =} fitlme (@var{tbl}, @var{formula})
## @deftypefnx {statistics} {@var{lme} =} fitlme (@dots{}, @var{name}, @var{value})
##
## Fit a linear mixed-effects model specified by a formula.
##
## @code{fitlme (@var{tbl}, @var{formula})} fits the linear mixed-effects model
## described by @var{formula} to the variables in the table @var{tbl}, and
## returns a @code{LinearMixedModel} object.
##
## @var{formula} is a character vector in Wilkinson notation extended with
## random-effects terms, for example @qcode{"y ~ x1 + x2 + (1 | g)"}.  The part
## to the left of @code{~} names the response; the fixed-effects part uses the
## usual operators (@code{+}, @code{*}, @code{:}, @code{^}, and @code{-1} to
## drop the intercept); and each random-effects term
## @code{(@var{expr} | @var{group})}
## adds random intercepts and slopes @var{expr} grouped by the factor
## @var{group} (or an interaction of factors, e.g.@: @code{g1:g2}).  As with
## fixed effects, a random intercept is implicit unless suppressed with
## @code{0} or @code{-1}.
##
## Rows of @var{tbl} with missing values in any model variable are removed
## before fitting.
##
## The following @var{name}/@var{value} pairs are accepted:
##
## @table @asis
## @item @qcode{"FitMethod"}
## The estimation criterion, @qcode{"ML"} (maximum likelihood, the default) or
## @qcode{"REML"} (restricted maximum likelihood).
## @end table
##
## @seealso{LinearMixedModel, fitlmematrix, fitlm, parseWilkinsonFormula}
## @end deftypefn

function lme = fitlme (tbl, formula, varargin)

  if (nargin < 2)
    print_usage ();
  endif
  if (! isa (tbl, "table"))
    error ("fitlme: TBL must be a table.");
  endif
  if (! (ischar (formula) || (isstring (formula) && isscalar (formula))))
    error ("fitlme: FORMULA must be a character vector.");
  endif
  formula = char (formula);

  ## --- options ---
  method = "ML";
  if (mod (numel (varargin), 2) != 0)
    error ("fitlme: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    switch (lower (char (varargin{i})))
      case "fitmethod"
        method = upper (char (varargin{i+1}));
        if (! any (strcmp (method, {"ML", "REML"})))
          error ("fitlme: FitMethod must be 'ML' or 'REML'.");
        endif
      otherwise
        error ("fitlme: unknown option '%s'.", char (varargin{i}));
    endswitch
  endfor

  ## --- decompose the formula ---
  S = parseWilkinsonFormula (formula, "mixed");
  if (! S.HasRandom)
    error (strcat ("fitlme: FORMULA must contain a random-effects term", ...
                   " '(...|...)'; use fitlm for fixed-effects models."));
  endif
  if (isempty (S.Response))
    error ("fitlme: FORMULA must specify a response variable.");
  endif
  if (! isempty (strfind (S.Response, ",")))
    error ("fitlme: multiple responses are not supported.");
  endif

  ## --- drop rows with missing values in any model variable ---
  vars = collect_vars (S);
  tvars = tbl.Properties.VariableNames;
  mask = true (height (tbl), 1);
  for i = 1:numel (vars)
    if (! ismember (vars{i}, tvars))
      error ("fitlme: variable '%s' is not in the table.", vars{i});
    endif
    col = tbl.(vars{i});
    if (isnumeric (col))
      mask = mask & all (! isnan (col), 2);
    endif
  endfor
  tbl = tbl(mask, :);

  ## --- fixed-effects design (reordered from the alphabetical 'model_matrix'
  ## column order to the formula's term order, to match MATLAB) ---
  [X, y, fenames] = parseWilkinsonFormula (S.FixedFormula, "model_matrix", tbl);
  perm = formula_order (fenames, S.FixedTerms, S.FixedIntercept);
  X = X(:, perm);
  fenames = fenames(perm);

  ## --- random-effects designs and grouping ---
  nt = numel (S.Random);
  Z = cell (1, nt);
  G = cell (1, nt);
  renames = cell (1, nt);
  grnames = cell (1, nt);
  for k = 1:nt
    rhs = terms_to_rhs (S.Random(k).Terms, S.Random(k).Intercept);
    [Zk, ~, znames] = parseWilkinsonFormula (["~ ", rhs], "model_matrix", tbl);
    Z{k} = Zk;
    renames{k} = znames;
    G{k} = combine_groups (tbl, S.Random(k).GroupVars);
    grnames{k} = S.Random(k).Group;
  endfor

  ## --- fit ---
  lme = fitlmematrix (X, y, Z, G, "FitMethod", method, ...
    "FixedEffectPredictors", fenames, "RandomEffectPredictors", renames, ...
    "RandomEffectGroups", grnames, "Formula", formula, ...
    "ResponseName", S.Response);

endfunction

## All variable names referenced by the decomposed formula.
function vars = collect_vars (S)
  vars = {};
  if (! isempty (S.Response))
    vars{end+1} = strtrim (S.Response);
  endif
  vars = [vars, flatten_terms(S.FixedTerms)];
  for k = 1:numel (S.Random)
    vars = [vars, flatten_terms(S.Random(k).Terms), S.Random(k).GroupVars];
  endfor
  vars = unique (vars);
endfunction

function names = flatten_terms (terms)
  names = {};
  for i = 1:numel (terms)
    names = [names, terms{i}];
  endfor
endfunction

## Permutation putting the alphabetically-ordered model_matrix columns into the
## formula's term order: intercept first, then each fixed term (matched by its
## set of variables), with any unmatched columns (e.g. categorical dummies)
## appended in their original order.
function perm = formula_order (names, terms, intercept)
  assigned = false (size (names));
  perm = [];
  if (intercept)
    j = find (strcmp (names, "(Intercept)"), 1);
    if (! isempty (j))
      perm(end+1) = j;
      assigned(j) = true;
    endif
  endif
  for t = 1:numel (terms)
    tvars = sort (terms{t});
    for j = 1:numel (names)
      if (assigned(j))
        continue;
      endif
      cvars = sort (strsplit (names{j}, ":"));
      if (isequal (cvars, tvars))
        perm(end+1) = j;
        assigned(j) = true;
      endif
    endfor
  endfor
  perm = [perm, find(! assigned)'];
endfunction

## Formula RHS string for a random-effects design expression.
function rhs = terms_to_rhs (terms, intercept)
  parts = cell (1, numel (terms));
  for i = 1:numel (terms)
    parts{i} = strjoin (terms{i}, ":");
  endfor
  if (intercept)
    if (isempty (parts))
      rhs = "1";
    else
      rhs = strjoin (parts, " + ");
    endif
  else
    rhs = [strjoin(parts, " + "), " - 1"];
  endif
endfunction

## Grouping index for one random-effects term: a single grouping variable, or
## the observed combinations of an interaction of grouping variables.
function g = combine_groups (tbl, gvars)
  if (numel (gvars) == 1)
    g = tbl.(gvars{1});
  else
    keys = [];
    for i = 1:numel (gvars)
      [~, ~, idx] = unique (tbl.(gvars{i}));
      keys = [keys, idx(:)];
    endfor
    [~, ~, g] = unique (keys, "rows");
  endif
endfunction

%!demo
%! ## Random-intercept model: sleep-study-like data with per-subject intercepts.
%! subject = reshape (repmat (1:6, 5, 1), [], 1);
%! days = repmat ((0:4)', 6, 1);
%! b0 = reshape (repmat ([1 -1 0.5 -0.5 0.2 -0.2], 5, 1), [], 1);
%! y = 250 + 10 * days + 15 * b0 + 3 * sin (1:30)';
%! tbl = table (y, days, subject);
%! lme = fitlme (tbl, "y ~ days + (1 | subject)", "FitMethod", "REML");
%! disp (lme.Coefficients);

## --- MATLAB-verified parity: the formula path reproduces the design-matrix fit
## and MATLAB's fitlme (R2026a) on the reference 42-row data set. ---
%!shared tbl, xL
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
%! g = [1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 ...
%!  5 6 1 2 3 4 5 6]';
%! g2 = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 ...
%!  2 3 1 2 3 1 2 3]';
%! tbl = table (yL, xL, x2, g, g2);

%!test  # every random-effects syntax form reproduces MATLAB's log-likelihood
%! forms = { "yL ~ xL + x2 + (1|g)", -49.016009; ...
%!           "yL ~ xL + x2 + (1 + xL|g)", -48.280962; ...
%!           "yL ~ xL + x2 + (1|g) + (xL-1|g)", -48.947341; ...
%!           "yL ~ xL + x2 + (xL-1|g)", -61.635296; ...
%!           "yL ~ xL + x2 + (-1 + xL|g)", -61.635296; ...
%!           "yL ~ xL + x2 + (1|g) + (1|g2)", -45.924711; ...
%!           "yL ~ xL + x2 + (1|g) + (1|g:g2)", -49.016009; ...
%!           "yL ~ xL*x2 + (1|g)", -50.366415; ...
%!           "yL ~ xL + x2 + (1 + xL + x2|g)", -48.162461 };
%! for k = 1:rows (forms)
%!   lme = fitlme (tbl, forms{k,1}, "FitMethod", "REML");
%!   assert (lme.LogLikelihood, forms{k,2}, 1e-4);
%! endfor

%!test  # random intercept via formula -- matches MATLAB fitlme
%! lme = fitlme (tbl, "yL ~ xL + x2 + (1 | g)", "FitMethod", "REML");
%! assert (isa (lme, "LinearMixedModel"));
%! assert (lme.Coefficients.Estimate, [1.96839; -1.37926; 0.811747], 1e-4);
%! assert (lme.Coefficients.SE, [0.376509; 0.113264; 0.0966571], 1e-5);
%! assert (lme.LogLikelihood, -49.016009, 1e-4);
%! [psi, mse] = covarianceParameters (lme);
%! assert (psi{1}, 0.79429917, 1e-3);
%! assert (mse, 0.38539058, 1e-3);

%!test  # correlated random intercept + slope via formula -- matches MATLAB
%! lme = fitlme (tbl, "yL ~ xL + x2 + (1 + xL | g)", "FitMethod", "REML");
%! assert (lme.Coefficients.Estimate, [2.0034354; -1.3374715; 0.83788802], 1e-3);
%! assert (lme.LogLikelihood, -48.280962, 1e-3);
%! psi = covarianceParameters (lme);
%! assert (psi{1}, [0.7846112, -0.14287556; -0.14287556, 0.026017248], 1e-3);

%!test  # formula fit equals the equivalent fitlmematrix fit
%! lme_f = fitlme (tbl, "yL ~ xL + x2 + (1 | g)", "FitMethod", "REML");
%! X = [ones(42,1), xL, tbl.x2];
%! lme_m = fitlmematrix (X, tbl.yL, ones (42, 1), tbl.g, "FitMethod", "REML");
%! assert (lme_f.Coefficients.Estimate, lme_m.Coefficients.Estimate, 1e-8);
%! assert (lme_f.LogLikelihood, lme_m.LogLikelihood, 1e-8);

%!test  # slope-only random term (no random intercept)
%! lme = fitlme (tbl, "yL ~ xL + x2 + (xL - 1 | g)", "FitMethod", "REML");
%! [psi, ~] = covarianceParameters (lme);
%! assert (isscalar (psi{1}));        # 1x1 covariance (slope only)

%!test  # metadata: Formula and ResponseName are stored
%! lme = fitlme (tbl, "yL ~ xL + (1 | g)");
%! assert (lme.Formula, "yL ~ xL + (1 | g)");
%! assert (lme.ResponseName, "yL");
%! assert (lme.FitMethod, "ML");      # default

%!test  # rows with missing values are dropped before fitting
%! t2 = tbl;
%! t2.yL(3) = NaN;  t2.xL(10) = NaN;
%! lme = fitlme (t2, "yL ~ xL + x2 + (1 | g)");
%! assert (lme.NumObservations, 40);

## Input validation
%!error <Invalid call> fitlme (table ())
%!error <TBL must be a table> fitlme (magic (3), "y ~ x + (1|g)")
%!error <FORMULA must be a character vector> fitlme (table (), 5)
%!error <must contain a random-effects term> fitlme (table ((1:3)', "VariableNames", {"y"}), "y ~ 1")
%!error <unknown option 'bogus'> fitlme (table ((1:3)', "VariableNames", {"y"}), "y ~ (1|y)", "bogus", 1)
%!error <FitMethod must be 'ML' or 'REML'> fitlme (table ((1:3)', "VariableNames", {"y"}), "y ~ (1|y)", "FitMethod", "xxx")
