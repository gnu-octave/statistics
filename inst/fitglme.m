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
## @deftypefn  {statistics} {@var{glme} =} fitglme (@var{tbl}, @var{formula})
## @deftypefnx {statistics} {@var{glme} =} fitglme (@dots{}, @var{name}, @var{value})
##
## Fit a generalized linear mixed-effects model specified by a formula.
##
## @code{fitglme (@var{tbl}, @var{formula})} fits the generalized linear
## mixed-effects model described by @var{formula} to the table @var{tbl} and
## returns a @code{GeneralizedLinearMixedModel} object.
##
## @var{formula} uses the same syntax as @code{fitlme}: a response, a
## fixed-effects part, and one or more random-effects terms
## @code{(@var{expr} | @var{group})}, for example
## @qcode{"y ~ x + (1 | g)"}.  The model is fitted by penalized
## quasi-likelihood.
##
## The following @var{name}/@var{value} pairs are accepted:
##
## @table @asis
## @item @qcode{"Distribution"}
## The response distribution: @qcode{"normal"} (default), @qcode{"binomial"}, or
## @qcode{"poisson"}.
##
## @item @qcode{"Link"}
## The link function: @qcode{"identity"}, @qcode{"logit"}, or @qcode{"log"}.
## The default is the canonical link of the chosen distribution.
##
## @item @qcode{"FitMethod"}
## @qcode{"MPL"} (maximum pseudo-likelihood, the default), @qcode{"REMPL"}
## (restricted MPL), @qcode{"Laplace"}, or @qcode{"ApproximateLaplace"}.  The
## first two differ in the pseudo-likelihood used for the covariance parameters;
## the last two report the Laplace-approximated marginal log-likelihood.
## @end table
##
## Only the canonical links and the full (unstructured) random-effects
## covariance are currently supported.
##
## @seealso{GeneralizedLinearMixedModel, fitlme, fitglm}
## @end deftypefn

function glme = fitglme (tbl, formula, varargin)

  if (nargin < 2)
    print_usage ();
  endif
  if (! isa (tbl, "table"))
    error ("fitglme: TBL must be a table.");
  endif
  if (! (ischar (formula) || (isstring (formula) && isscalar (formula))))
    error ("fitglme: FORMULA must be a character vector.");
  endif
  formula = char (formula);

  ## --- options ---
  distr = "normal";
  link = "";
  method = "MPL";
  if (mod (numel (varargin), 2) != 0)
    error ("fitglme: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    switch (lower (char (varargin{i})))
      case "distribution"
        distr = lower (char (varargin{i+1}));
      case "link"
        link = lower (char (varargin{i+1}));
      case "fitmethod"
        method = char (varargin{i+1});
      otherwise
        error ("fitglme: unknown option '%s'.", char (varargin{i}));
    endswitch
  endfor
  if (! any (strcmp (distr, {"normal", "binomial", "poisson"})))
    error ("fitglme: Distribution must be 'normal', 'binomial', or 'poisson'.");
  endif
  if (isempty (link))
    link = canonical_link (distr);
  endif
  mkey = lower (method);
  if (! any (strcmp (mkey, {"mpl", "rempl", "laplace", "approximatelaplace"})))
    error (strcat ("fitglme: FitMethod must be 'MPL', 'REMPL', 'Laplace',", ...
                   " or 'ApproximateLaplace'."));
  endif

  ## --- decompose the formula ---
  S = parseWilkinsonFormula (formula, "mixed");
  if (! S.HasRandom)
    error (strcat ("fitglme: FORMULA must contain a random-effects term", ...
                   " '(...|...)'; use fitglm for fixed-effects models."));
  endif
  if (isempty (S.Response))
    error ("fitglme: FORMULA must specify a response variable.");
  endif

  ## --- drop rows with missing values in any model variable ---
  vars = collect_vars (S);
  tvars = tbl.Properties.VariableNames;
  mask = true (height (tbl), 1);
  for i = 1:numel (vars)
    if (! ismember (vars{i}, tvars))
      error ("fitglme: variable '%s' is not in the table.", vars{i});
    endif
    col = tbl.(vars{i});
    if (isnumeric (col))
      mask = mask & all (! isnan (col), 2);
    endif
  endfor
  tbl = tbl(mask, :);

  ## --- fixed design (formula-term order) and random designs ---
  [X, y, fenames] = parseWilkinsonFormula (S.FixedFormula, "model_matrix", tbl);
  perm = formula_order (fenames, S.FixedTerms, S.FixedIntercept);
  X = X(:, perm);
  fenames = fenames(perm);

  nt = numel (S.Random);
  Z = cell (1, nt);  G = cell (1, nt);
  renames = cell (1, nt);  grnames = cell (1, nt);
  for k = 1:nt
    rhs = terms_to_rhs (S.Random(k).Terms, S.Random(k).Intercept);
    [Zk, ~, znames] = parseWilkinsonFormula (["~ ", rhs], "model_matrix", tbl);
    Z{k} = Zk;  renames{k} = znames;
    G{k} = combine_groups (tbl, S.Random(k).GroupVars);
    grnames{k} = S.Random(k).Group;
  endfor

  ## --- fit ---
  fit = __glmefit__ (X, y, Z, G, distr, link, mkey);
  info = fit;
  info.X = X;  info.y = y;
  info.CoefficientNames = fenames;
  info.GroupNames = grnames;
  info.REPred = renames;
  info.FitMethod = method;
  info.Formula = formula;
  info.ResponseName = S.Response;
  glme = GeneralizedLinearMixedModel (info);

endfunction

function lnk = canonical_link (distr)
  switch (distr)
    case "normal",   lnk = "identity";
    case "binomial", lnk = "logit";
    case "poisson",  lnk = "log";
  endswitch
endfunction

## --- formula helpers (shared conventions with fitlme) ---
function vars = collect_vars (S)
  vars = {};
  if (! isempty (S.Response)), vars{end+1} = strtrim (S.Response); endif
  vars = [vars, flatten_terms(S.FixedTerms)];
  for k = 1:numel (S.Random)
    vars = [vars, flatten_terms(S.Random(k).Terms), S.Random(k).GroupVars];
  endfor
  vars = unique (vars);
endfunction

function names = flatten_terms (terms)
  names = {};
  for i = 1:numel (terms), names = [names, terms{i}]; endfor
endfunction

function perm = formula_order (names, terms, intercept)
  assigned = false (size (names));
  perm = [];
  if (intercept)
    j = find (strcmp (names, "(Intercept)"), 1);
    if (! isempty (j)), perm(end+1) = j; assigned(j) = true; endif
  endif
  for t = 1:numel (terms)
    tvars = sort (terms{t});
    for j = 1:numel (names)
      if (assigned(j)), continue; endif
      if (isequal (sort (strsplit (names{j}, ":")), tvars))
        perm(end+1) = j;  assigned(j) = true;
      endif
    endfor
  endfor
  perm = [perm, find(! assigned)'];
endfunction

function rhs = terms_to_rhs (terms, intercept)
  parts = cell (1, numel (terms));
  for i = 1:numel (terms), parts{i} = strjoin (terms{i}, ":"); endfor
  if (intercept)
    if (isempty (parts)), rhs = "1"; else, rhs = strjoin (parts, " + "); endif
  else
    rhs = [strjoin(parts, " + "), " - 1"];
  endif
endfunction

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
%! ## Poisson mixed model with a random intercept per group.
%! g = reshape (repmat (1:6, 7, 1), [], 1);
%! x = randn (42, 1);
%! y = poissrnd (exp (0.3 + 0.5 * x + 0.2 * reshape (repmat (randn (1, 6), 7, 1), [], 1)));
%! tbl = table (y, x, g);
%! glme = fitglme (tbl, "y ~ x + (1 | g)", "Distribution", "poisson");
%! disp (glme.Coefficients);

## MATLAB-verified parity (fitglme R2026a) on the reference GLME data.
%!shared tbl
%! xL = [0.032760004 0.70410822 -0.8646718 -0.28869454 0.51276678 -1.4975462 ...
%!  -1.4527871 -0.80013541 -1.644209 1.5137701 0.72905543 0.20880758 1.0856145 ...
%!  0.62862577 -0.87409978 1.9178276 0.09748204 0.50697633 1.0247569 ...
%!  -0.92789896 -0.88921018 -0.98322849 -0.031378913 0.86875961 -0.91481141 ...
%!  0.034324163 -0.25025257 -1.0575644 -0.86131607 -0.35355444 0.82950729 ...
%!  -0.36874363 0.061580868 0.55803564 -0.1763803 1.0482413 1.0137831 ...
%!  -0.94876976 -0.010703972 -0.35149845 -1.6828735 -1.0493301]';
%! yBin = [0 0 0 0 1 0 0 1 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 1 1 ...
%!  0 1 1 0 1 1 0 0]';
%! yPois = [3 3 1 1 1 1 1 2 1 5 2 0 5 0 1 5 0 2 0 0 1 0 0 5 3 0 1 0 0 1 1 1 2 2 ...
%!  1 1 4 0 1 0 0 1]';
%! g = [1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 ...
%!  6 1 2 3 4 5 6]';
%! tbl = table (yBin, yPois, xL, g);

%!test  # binomial logit, MPL -- matches MATLAB fitglme
%! glme = fitglme (tbl, "yBin ~ xL + (1 | g)", "Distribution", "binomial");
%! assert (isa (glme, "GeneralizedLinearMixedModel"));
%! assert (glme.Coefficients.Estimate, [-0.55912; 0.76062], 1e-3);
%! assert (glme.Coefficients.SE, [0.33856; 0.39971], 1e-3);
%! assert (glme.LogLikelihood, -92.58872, 1e-2);

%!test  # poisson log, REMPL -- non-degenerate random-effect variance
%! glme = fitglme (tbl, "yPois ~ xL + (1 | g)", "Distribution", "poisson", ...
%!                 "FitMethod", "REMPL");
%! assert (glme.Coefficients.Estimate, [0.23092; 0.67809], 1e-3);
%! [psi, ~] = covarianceParameters (glme);
%! assert (psi{1}, 0.015918, 1e-3);
%! assert (glme.LogLikelihood, -56.90298, 1e-2);

%!test  # Laplace reports the marginal log-likelihood
%! glme = fitglme (tbl, "yBin ~ xL + (1 | g)", "Distribution", "binomial", ...
%!                 "FitMethod", "Laplace");
%! assert (glme.LogLikelihood, -25.36010, 1e-2);

%!test  # coefficient stats: DF = n - p, tStat = Estimate / SE
%! glme = fitglme (tbl, "yPois ~ xL + (1 | g)", "Distribution", "poisson");
%! C = glme.Coefficients;
%! assert (C.DF, [40; 40]);
%! assert (C.tStat, C.Estimate ./ C.SE, 1e-10);

%!test  # metadata
%! glme = fitglme (tbl, "yBin ~ xL + (1 | g)", "Distribution", "binomial");
%! assert (glme.Distribution, "binomial");
%! assert (glme.Link, "logit");
%! assert (glme.FitMethod, "MPL");
%! assert (glme.ResponseName, "yBin");

## Input validation
%!error <Invalid call> fitglme (table ())
%!error <TBL must be a table> fitglme (magic (3), "y ~ x + (1|g)")
%!error <Distribution must be> fitglme (table ((1:3)', "VariableNames", {"y"}), "y ~ (1|y)", "Distribution", "xxx")
%!error <must contain a random-effects term> fitglme (table ((1:3)', "VariableNames", {"y"}), "y ~ 1")
%!error <unknown option 'bogus'> fitglme (table ((1:3)', "VariableNames", {"y"}), "y ~ (1|y)", "bogus", 1)
