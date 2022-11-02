## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
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
## @deftypefn {Function File} @var{tab} = fitlm (@var{X}, @var{y})
## @deftypefnx {Function File} @var{tab} = fitlm (@var{X}, @var{y}, "name", @var{value})
## @deftypefnx {Function File} @var{tab} = fitlm (@var{X}, @var{y}, @var{modelspec})
## @deftypefnx {Function File} @var{tab} = fitlm (@var{X}, @var{y}, @var{modelspec}, "name", @var{value})
## @deftypefnx {Function File} [@var{tab}] = fitlm (...)
## @deftypefnx {Function File} [@var{tab}, @var{stats}] = fitlm (...)
## @deftypefnx {Function File} [@var{tab}, @var{stats}] = fitlm (...)
##
## Regress the continuous outcome (i.e. dependent variable) @var{y} on
## continuous or categorical predictors (i.e. independent variables) @var{X}
## by minimizing the sum-of-squared residuals. Unless requested otherwise,
## @qcode{fitlm} prints the model formula, the regression coefficients (i.e. 
## parameters/contrasts) and an ANOVA table. Note that unlike @qcode{anovan}, 
## @qcode{fitlm} treats all factors as continuous by default. 
##
## @var{X} must be a column major matrix or cell array consisting of the
## predictors. By default, there is a constant term in the model, unless you,
## explicitly remove it, so do not include a column of 1s in X. @var{y} must be
## a column vector corresponding to the outcome variable. @var{modelspec} can
## specified as one of the following:
##
## @itemize
## @item
## "constant" : model contains only a constant (intercept) term.
##
## @item
## "linear" (default) : model contains an intercept and linear term for each
## predictor.
##
## @item
## "interactions" : model contains an intercept, linear term for each predictor
## and all products of pairs of distinct predictors.
##
## @item
## "full" : model contains an intercept, linear term for each predictor and
## all combinations of the predictors.
##
## @item
## a matrix of term definitions : an t-by-(N+1) matrix specifying terms in
## a model, where t is the number of terms, N is the number of predictor 
## variables, and +1 accounts for the outcome variable. The outcome variable
## is the last column in the terms matrix and must be a column of zeros.
## An intercept must be specified in the first row of the terms matrix and
## must be a row of zeros.
## @end itemize
##
## @qcode{fitlm} can take a number of optional parameters as name-value pairs.
##
## @code{[@dots{}] = fitlm (..., "CategoricalVars", @var{categorical})}
##
## @itemize
## @item
## @var{categorical} is a vector of indices indicating which of the columns 
## (i.e. variables) in @var{X} should be treated as categorical predictors
## rather than as continuous predictors. 
## @end itemize
##
## @qcode{fitlm} also accepts optional @qcode{anovan} parameters as name-value
## pairs (except for the "model" parameter). The accepted parameter names from
## @qcode{anovan} and their default values in @qcode{fitlm} are:
##
## @itemize
## @item
## @var{CONTRASTS} : "treatment"
##
## @item
## @var{SSTYPE}: 2
##
## @item
## @var{ALPHA}: 0.05
##
## @item
## @var{DISPLAY}: "on"
##
## @item
## @var{WEIGHTS}: [] (empty)
##
## @item
## @var{RANDOM}: [] (empty)
##
## @item
## @var{CONTINUOUS}: [1:N]
##
## @item
## @var{VARNAMES}: [] (empty)
## @end itemize
##
## Type '@qcode{help anovan}' to find out more about what these options do.
##
## @qcode{fitlm} can return up to two output arguments:
##
## [@var{tab}] = fitlm (@dots{}) returns a cell array containing a
## table of model parameters
##
## [@var{tab}, @var{stats}] = fitlm (@dots{}) returns a structure
## containing additional statistics, including degrees of freedom and effect
## sizes for each term in the linear model, the design matrix, the
## variance-covariance matrix, (weighted) model residuals, and the mean squared
## error. The columns of @var{stats}.coeffs (from left-to-right) report the
## model coefficients, standard errors, lower and upper 100*(1-alpha)%
## confidence interval bounds, t-statistics, and p-values relating to the
## contrasts. The number appended to each term name in @var{stats}.coeffnames
## corresponds to the column number in the relevant contrast matrix for that
## factor. The @var{stats} structure can be used as input for @qcode{multcompare}.
## The @var{stats} structure is recognised by the functions @qcode{bootcoeff}
## and @qcode{bootemm} from the statistics-bootstrap package. Note that if the
## model contains a continuous variable and you wish to use the @var{STATS}
## output as input to @qcode{multcompare}, then the model needs to be refit 
## with the "contrast" parameter set to a sum-to-zero contrast coding scheme,
## e.g."simple". 
##
## @seealso{anovan, multcompare}
## @end deftypefn

function [T, STATS] = fitlm (X, y, varargin)

    ## Check input and output arguments
    if (nargin < 2)
      error (strcat (["fitlm usage: ""fitlm (X, y, varargin)""; "], ...
                      [" atleast 2 input arguments required"]));
    endif
    if (nargout > 3)
      error ("fitlm: invalid number of output arguments requested");
    endif

    ## Evaluate input data 
    [n, N] = size (X);

    ## Fetch anovan options
    options = varargin;
    if (isempty(options))
      options{1} = "linear";
    endif

    ## Check if MODELSPEC was provided. If not create it.
    if (ischar (options{1}))
      if (! ismember (lower (options{1}), {"sstype", "varnames", "contrasts", ...
                    "weights", "alpha", "display", "continuous", ...
                    "categorical", "categoricalvars", "random", "model"}))
        MODELSPEC = options{1};
        options(1) = [];
        CONTINUOUS = [];
      else
        ## If MODELSPEC is not provided, set it for an additive linear model  
        MODELSPEC = zeros (N + 1);
        MODELSPEC(2:N+1, 1:N) = eye (N);
      end
    else
      MODELSPEC = options{1};
      options(1) = [];
    endif

    ## Evaluate MODELSPEC
    if (ischar (MODELSPEC))
      MODELSPEC = lower (MODELSPEC);
      if (! isempty (regexp (MODELSPEC, "~")))
        error ("fitlm: model formulae are not a supported format for MODELSPEC")
      endif
      if (! ismember (MODELSPEC, {"constant", "linear", "interaction", ...
                                  "interactions", "full"}))
        error ("fitlm: character vector for model specification not recognised")
      endif
      if strcmp (MODELSPEC, "constant") 
        X = [];
        MODELSPEC = "linear";
        N = 0;
      endif
    else
      if (size (MODELSPEC, 1) < N + 1)
        error ("fitlm: number of rows in MODELSPEC must  1 + number of columns in X");
      endif
      if (size (MODELSPEC, 2) != N + 1)
        error ("fitlm: number of columns in MODELSPEC must = 1 + number of columns in X");
      endif
      if (! all (ismember (MODELSPEC(:), [0,1])))
        error (strcat (["fitlm: elements of the model terms matrix must be "], ...
                       [" either 0 or 1. Higher order terms are not supported"]));
      endif
      MODELSPEC = logical (MODELSPEC(2:N+1,1:N));
    endif

    ## Check for unsupported options used by anovan
    if (any (strcmpi ("MODEL", options)))
      error (strcat(["fitlm: modelspec should be specified in the third"], ...
                    [" input argument of fitlm (if at all)"]));
    endif

    ## Check and set variable types
    idx = find (any (cat (1, strcmpi ("categorical", options), ...
                     strcmpi ("categoricalvars", options))));
    if (! isempty (idx))
      CONTINUOUS = [1:N];
      CONTINUOUS(ismember(CONTINUOUS,options{idx+1})) = [];
      options(idx:idx+1) = [];
    else
      CONTINUOUS = [1:N];
    endif
    idx = find (strcmpi ("continuous", options));
    if (! isempty (idx))
      ## Note that setting continuous parameter will override settings made
      ## to "categorical"
      CONTINUOUS = options{idx+1};
    endif

    ## Check if anovan CONTRASTS option was used
    idx = find (strcmpi ("contrasts", options));
    if (isempty (idx))
      CONTRASTS = "treatment";
    else
      CONTRASTS = options{idx+1};
      if (ischar(CONTRASTS))
        contr_str = CONTRASTS;
        CONTRASTS = cell (1, N);
        CONTRASTS(:) = {contr_str};
      endif
      if (! iscell (CONTRASTS))
        CONTRASTS = {CONTRASTS};
      endif
      for i = 1:N
        if (! isnumeric(CONTRASTS{i}))
          if (! isempty (CONTRASTS{i}))
            if (! ismember (CONTRASTS{i}, ...
                            {"simple","poly","helmert","effect","treatment"}))
              error (strcat(["fitlm: the choices for built-in contrasts are"], ...
                     [" ""simple"", ""poly"", ""helmert"", ""effect"", or ""treatment"""]));
            endif
          endif
        endif
      endfor
    endif

    ## Check if anovan SSTYPE option was used
    idx = find (strcmpi ("sstype", options));
    if (isempty (idx))
      SSTYPE = 2;
    else
      SSTYPE = options{idx+1};
    endif

    ## Perform model fit and ANOVA
    [jnk, jnk, STATS] = anovan (y, X, options{:}, ...
                            "model", MODELSPEC, ...
                            "contrasts", CONTRASTS, ...
                            "continuous", CONTINUOUS, ...
                            "sstype", SSTYPE);

    ## Create table of regression coefficients
    ncoeff = sum (STATS.df);
    T = cell (2 + ncoeff, 7);
    T(1,:) = {"Parameter", "Estimate", "SE", "Lower.CI", "Upper.CI", "t", "Prob>|t|"};
    T(2:end,1) = STATS.coeffnames;
    T(2:end,2:7) = num2cell (STATS.coeffs);

    ## Update STATS structure
    STATS.source = "fitlm";

endfunction

%!demo
%! y =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! X = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%!
%! [TAB,STATS] = fitlm (X,y,"linear","CategoricalVars",1,"display","on");

%!demo
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! brands = {'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'};
%! popper = {'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; ...
%!           'air', 'air', 'air'; 'air', 'air', 'air'; 'air', 'air', 'air'};
%!
%! [TAB, STATS] = fitlm ({brands(:),popper(:)},popcorn(:),"interactions",...
%!                          "CategoricalVars",[1,2],"display","on");

%!test
%! y =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! X = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%! [TAB,STATS] = fitlm (X,y,"continuous",[],"display","off");
%! [TAB,STATS] = fitlm (X,y,"CategoricalVars",1,"display","off");
%! [TAB,STATS] = fitlm (X,y,"constant","categorical",1,"display","off");
%! [TAB,STATS] = fitlm (X,y,"linear","categorical",1,"display","off");
%! [TAB,STATS] = fitlm (X,y,[0,0;1,0],"categorical",1,"display","off");
%! assert (TAB{2,2}, 10, 1e-04);
%! assert (TAB{3,2}, 7.99999999999999, 1e-09);
%! assert (TAB{4,2}, 8.99999999999999, 1e-09);
%! assert (TAB{5,2}, 11.0001428571429, 1e-09);
%! assert (TAB{6,2}, 19.0001111111111, 1e-09);
%! assert (TAB{2,3}, 1.01775379540949, 1e-09);
%! assert (TAB{3,3}, 1.64107868458008, 1e-09);
%! assert (TAB{4,3}, 1.43932122062479, 1e-09);
%! assert (TAB{5,3}, 1.48983900477565, 1e-09);
%! assert (TAB{6,3}, 1.3987687997822, 1e-09);
%! assert (TAB{2,6}, 9.82555903510687, 1e-09);
%! assert (TAB{3,6}, 4.87484242844031, 1e-09);
%! assert (TAB{4,6}, 6.25294748040552, 1e-09);
%! assert (TAB{5,6}, 7.38344399756088, 1e-09);
%! assert (TAB{6,6}, 13.5834536158296, 1e-09);
%! assert (TAB{3,7}, 2.85812420217862e-05, 1e-12);
%! assert (TAB{4,7}, 5.22936741204002e-07, 1e-06);
%! assert (TAB{5,7}, 2.12794763209106e-08, 1e-07);
%! assert (TAB{6,7}, 7.82091664406755e-15, 1e-08);

%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! brands = bsxfun (@times, ones(6,1), [1,2,3]);
%! popper = bsxfun (@times, [1;1;1;2;2;2], ones(1,3));
%!
%! [TAB, STATS] = fitlm ({brands(:),popper(:)},popcorn(:),"interactions",...
%!                          "categoricalvars",[1,2],"display","off");
%! assert (TAB{2,2}, 5.66666666666667, 1e-09);
%! assert (TAB{3,2}, -1.33333333333333, 1e-09);
%! assert (TAB{4,2}, -2.16666666666667, 1e-09);
%! assert (TAB{5,2}, 1.16666666666667, 1e-09);
%! assert (TAB{6,2}, -0.333333333333334, 1e-09);
%! assert (TAB{7,2}, -0.166666666666667, 1e-09);
%! assert (TAB{2,3}, 0.215165741455965, 1e-09);
%! assert (TAB{3,3}, 0.304290309725089, 1e-09);
%! assert (TAB{4,3}, 0.304290309725089, 1e-09);
%! assert (TAB{5,3}, 0.304290309725089, 1e-09);
%! assert (TAB{6,3}, 0.43033148291193, 1e-09);
%! assert (TAB{7,3}, 0.43033148291193, 1e-09);
%! assert (TAB{2,6}, 26.3362867542108, 1e-09);
%! assert (TAB{3,6}, -4.38178046004138, 1e-09);
%! assert (TAB{4,6}, -7.12039324756724, 1e-09);
%! assert (TAB{5,6}, 3.83405790253621, 1e-09);
%! assert (TAB{6,6}, -0.774596669241495, 1e-09);
%! assert (TAB{7,6}, -0.387298334620748, 1e-09);
%! assert (TAB{2,7}, 5.49841502258254e-12, 1e-09);
%! assert (TAB{3,7}, 0.000893505495903642, 1e-09);
%! assert (TAB{4,7}, 1.21291454302428e-05, 1e-09);
%! assert (TAB{5,7}, 0.00237798044119407, 1e-09);
%! assert (TAB{6,7}, 0.453570536021938, 1e-09);
%! assert (TAB{7,7}, 0.705316781644046, 1e-09);
%! ## Test with string ids for categorical variables
%! brands = {'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'};
%! popper = {'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; ...
%!           'air', 'air', 'air'; 'air', 'air', 'air'; 'air', 'air', 'air'};
%! [TAB, STATS] = fitlm ({brands(:),popper(:)},popcorn(:),"interactions",...
%!                          "categoricalvars",[1,2],"display","off");

%!test 
%! load carsmall
%! X = [Weight,Horsepower,Acceleration];
%! [TAB, STATS] = fitlm (X, MPG,"constant","display","off");
%! [TAB, STATS] = fitlm (X, MPG,"linear","display","off");
%! assert (TAB{2,2}, 47.9767628118615, 1e-09);
%! assert (TAB{3,2}, -0.00654155878851796, 1e-09);
%! assert (TAB{4,2}, -0.0429433065881864, 1e-09);
%! assert (TAB{5,2}, -0.0115826516894871, 1e-09);
%! assert (TAB{2,3}, 3.87851641748551, 1e-09);
%! assert (TAB{3,3}, 0.00112741016370336, 1e-09);
%! assert (TAB{4,3}, 0.0243130608813806, 1e-09);
%! assert (TAB{5,3}, 0.193325043113178, 1e-09);
%! assert (TAB{2,6}, 12.369874881944, 1e-09);
%! assert (TAB{3,6}, -5.80228828790225, 1e-09);
%! assert (TAB{4,6}, -1.76626492228599, 1e-09);
%! assert (TAB{5,6}, -0.0599128364487485, 1e-09);
%! assert (TAB{2,7}, 4.89570341688996e-21, 1e-09);
%! assert (TAB{3,7}, 9.87424814144e-08, 1e-09);
%! assert (TAB{4,7}, 0.0807803098213114, 1e-09);
%! assert (TAB{5,7}, 0.952359384151778, 1e-09);
