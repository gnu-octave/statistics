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
## @deftp {statistics} {} GeneralizedLinearModel
##
## Generalized linear regression model class.
##
## A @code{GeneralizedLinearModel} object encapsulates a generalized linear
## model (GLM) of a response on one or more predictors, fitted by iteratively
## reweighted least squares.  It is the GLM counterpart of @code{LinearModel}
## and is normally created with the @code{fitglm} function.
##
## The response is modelled through a distribution from the exponential family
## (@qcode{'normal'}, @qcode{'binomial'}, @qcode{'poisson'}, @qcode{'gamma'}, or
## @qcode{'inverse gaussian'}) and a link function @math{g} relating the mean
## @math{mu} to the linear predictor @math{eta = g (mu)}.
##
## The most useful properties are @code{Coefficients} (a table of estimates,
## standard errors, @math{t}-statistics and p-values), @code{Deviance},
## @code{Dispersion}, @code{Residuals}, @code{Fitted}, @code{Diagnostics},
## @code{Distribution}, and @code{Link}.  @code{ObservationInfo} records which
## rows were weighted, excluded, or missing, and @code{Variables} holds the data
## the model was built from.  Fitted models support the @code{predict} and
## @code{feval} methods for prediction.
##
## @code{Fitted}, @code{Residuals}, @code{Diagnostics}, and
## @code{ObservationInfo} have one row per @emph{input} observation, not per
## fitted observation.  Rows that were excluded with the @qcode{'Exclude'} pair
## still carry a fitted value and a residual, since the model can be evaluated
## there; rows dropped because a variable was missing carry @code{NaN}.
##
## For a binomial response given together with a @qcode{'BinomialSize'} vector
## @math{N}, the response is the @emph{proportion} of successes, as
## @code{fitglm} documents.  @code{Fitted.Response} is then the fitted count
## @math{N p} and @code{Residuals.Raw} is on that same count scale.
##
## A categorical predictor expands to indicator columns, one per level bar the
## reference level, which the intercept carries.  When the model has no
## intercept, the @emph{first} categorical predictor is given an indicator for
## every one of its levels instead, so that its coefficients are the group
## means; any further categorical predictor stays reference coded, which keeps
## the design full rank.  This differs from MATLAB, which omits the reference
## level whether or not an intercept is present and so cannot fit the reference
## group at all -- for a three-level grouping variable @code{g}, MATLAB fits
## @code{y ~ g - 1} with two coefficients, predicts exactly 0 for every
## observation in the omitted group, and reports a negative @math{R^2}.  This
## implementation returns three coefficients, one per group.
##
## @seealso{fitglm, LinearModel, glmfit, glmval}
## @end deftp

classdef GeneralizedLinearModel

  properties (GetAccess = public, SetAccess = protected)

    ## Table of coefficient estimates, standard errors, t-statistics, p-values.
    Coefficients = [];

    ## Cell array of coefficient names.
    CoefficientNames = {};

    ## Estimated covariance matrix of the coefficients.
    CoefficientCovariance = [];

    ## Number of coefficients in the model.
    NumCoefficients = [];

    ## Number of estimated (nonzero-freedom) coefficients.
    NumEstimatedCoefficients = [];

    ## Number of predictor variables.
    NumPredictors = [];

    ## Number of observations used in the fit (missing/excluded rows removed).
    NumObservations = [];

    ## Deviance of the fitted model.
    Deviance = [];

    ## Error (residual) degrees of freedom.
    DFE = [];

    ## Dispersion parameter (estimated or fixed at 1).
    Dispersion = [];

    ## True if the dispersion parameter was estimated from the data.
    DispersionEstimated = [];

    ## Structure describing the response distribution.
    Distribution = [];

    ## Structure describing the link function (Name, Link, Derivative, Inverse).
    Link = [];

    ## Table of fitted values (Response on the mean scale, LinearPredictor, and
    ## Probability for the binomial family).  One row per input observation.
    Fitted = [];

    ## Table of residuals (Raw, LinearPredictor, Pearson, Anscombe, Deviance).
    ## One row per input observation.
    Residuals = [];

    ## Table of per-observation diagnostics (Leverage, CooksDistance,
    ## HatMatrix).  One row per input observation.
    Diagnostics = [];

    ## Log-likelihood of the fitted model.
    LogLikelihood = [];

    ## Structure of information criteria (AIC, AICc, BIC, CAIC).
    ModelCriterion = [];

    ## Structure of R-squared measures (Ordinary, Adjusted, Deviance,
    ## AdjGeneralized, LLR).
    Rsquared = [];

    ## Error (residual) sum of squares on the response scale.
    SSE = [];

    ## Regression sum of squares on the response scale.
    SSR = [];

    ## Total sum of squares on the response scale.
    SST = [];

    ## Offset vector added to the linear predictor, one element per input
    ## observation (all zero when no offset was given).
    Offset = [];

    ## Penalty applied to the likelihood; always @qcode{"none"} here.
    LikelihoodPenalty = [];

    ## Name of the response variable.
    ResponseName = 'y';

    ## Cell array of predictor variable names.
    PredictorNames = {};

    ## Cell array of all variable names.
    VariableNames = {};

    ## Number of variables (predictors and response).
    NumVariables = [];

    ## Table of per-variable information (class, range, in-model, categorical).
    VariableInfo = [];

    ## Table of the data the model was built from.
    Variables = [];

    ## Table recording, per input observation, its weight and whether it was
    ## excluded, missing, or used in the fit.
    ObservationInfo = [];

    ## Cell array of observation names; empty unless the data carried row names.
    ObservationNames = {};

    ## LinearFormula object describing the model formula.
    Formula = [];

    ## Structure describing the stepwise term-selection history.  Empty unless
    ## the model was built by @code{stepwiseglm}.
    Steps = [];

  endproperties

  properties (Access = private, Hidden)
    b_          = [];   # coefficient vector aligned with the design columns
    stats_      = [];   # glmfit stats struct (for prediction CIs)
    distr_      = '';   # distribution name
    linkarg_    = [];   # link specification for glmfit/glmval
    binomsize_  = [];   # BinomialSize (trials) for the binomial family
    terms_      = [];   # terms matrix (for rebuilding the design in predict)
    catinfo_    = [];   # categorical level info (names + levels)
    encnames_   = {};   # encoded predictor column names
    prednames_  = {};   # raw predictor names
    nulldev_    = [];   # deviance of the intercept-only model
    llnull_     = [];   # log-likelihood of the intercept-only model
    design_     = [];   # fitted design matrix (for diagnostics/plots)
    subset_     = [];   # logical mask of the rows used in the fit
    formulastr_ = '';   # rendered formula, e.g. 'log(y) ~ 1 + x1 + x2'
    xmeans_     = [];   # mean of each raw predictor (for slice/effect plots)
  endproperties

  methods (Hidden)

    ## Custom display of the object name.
    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ("%s =\n", in_name);
      endif
      disp (this);
    endfunction

    ## Custom display of the model summary.
    function disp (this)
      fprintf ("\n  Generalized linear regression model:\n");
      if (! isempty (this.formulastr_))
        fprintf ("      %s\n", this.formulastr_);
      endif
      if (! isempty (this.Distribution) && isstruct (this.Distribution))
        fprintf ("      Distribution = %s,  Link = %s\n", ...
                 this.Distribution.Name, this.Link.Name);
      endif
      if (! isempty (this.Coefficients))
        fprintf ("\n  Coefficients:\n\n");
        disp (this.Coefficients);
      endif
      fprintf ("\n");
      if (! isempty (this.NumObservations) && ! isempty (this.DFE))
        fprintf (strcat ("Number of observations: %d,", ...
                         " Error degrees of freedom: %d\n"), ...
                 this.NumObservations, this.DFE);
      endif
      if (! isempty (this.Dispersion))
        fprintf ("Dispersion: %g\n", this.Dispersion);
      endif
      if (! isempty (this.Deviance))
        fprintf ("Deviance: %g\n", this.Deviance);
      endif
      if (! isempty (this.NumCoefficients) && this.NumCoefficients > 1)
        df1  = this.NumCoefficients - 1;
        drop = this.nulldev_ - this.Deviance;
        if (this.DispersionEstimated)
          Fstat = (drop / df1) / this.Dispersion;
          pval  = 1 - fcdf (Fstat, df1, this.DFE);
          fprintf ("F-statistic vs. constant model: %g, p-value = %g\n", ...
                   Fstat, pval);
        else
          pval = 1 - chi2cdf (drop, df1);
          fprintf (strcat ("Chi^2-statistic vs. constant model: %g,", ...
                           " p-value = %g\n"), drop, pval);
        endif
      endif
    endfunction

    ## Class specific subscripted reference.
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("GeneralizedLinearModel: () indexing is not", ...
                         " supported.  Use dot notation for properties."));
        case '{}'
          error (strcat ("GeneralizedLinearModel: {} indexing is not", ...
                         " supported.  Use dot notation for properties."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("GeneralizedLinearModel.subsref: property name", ...
                           " must be a character vector."));
          endif
          if (ismethod (this, s.subs))
            [varargout{1:nargout}] = builtin ('subsref', this, [s, chain_s]);
            return;
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("GeneralizedLinearModel.subsref: unknown", ...
                           " property '%s'."), s.subs);
          end_try_catch
      endswitch
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

    ## Attach a stepwise-selection history structure.  Used by @code{stepwiseglm}
    ## to record the term-selection trace on the returned object; not intended
    ## for direct use.
    function this = setSteps (this, steps)
      this.Steps = steps;
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{mdl} =} GeneralizedLinearModel (@var{data}, @var{resp}, @var{modelspec})
    ## @deftypefnx {GeneralizedLinearModel} {@var{mdl} =} GeneralizedLinearModel (@dots{}, @var{Name}, @var{Value})
    ##
    ## Fit a generalized linear model.  Prefer the @code{fitglm} function, which
    ## documents the accepted inputs and @var{Name}/@var{Value} pairs.
    ##
    ## @end deftypefn
    function this = GeneralizedLinearModel (data, resp, modelspec, varargin)

      if (nargin == 0)
        return;   # empty object
      endif
      if (nargin < 3)
        error ("GeneralizedLinearModel: DATA, RESP, and MODELSPEC are required.");
      endif

      opts       = glm_parse_nv (varargin);
      is_formula = ischar (modelspec) && any (modelspec == '~');
      tbl        = [];   # defined only for table input; passed through helpers

      ## ------------------------------------------------------------------ ##
      ## Intake: resolve predictors, response, names, and the response vector.
      ## ------------------------------------------------------------------ ##
      if (! istable (data))
        if (! (isnumeric (data) && isreal (data) && ismatrix (data)))
          error ("GeneralizedLinearModel: X must be a real matrix.");
        endif
        if (! (isnumeric (resp) && isreal (resp) && isvector (resp)))
          error ("GeneralizedLinearModel: Y must be a real vector.");
        endif
        X_raw   = double (data);
        n_total = rows (X_raw);
        p_raw   = columns (X_raw);
        y_full  = double (resp(:));
        if (rows (X_raw) != numel (y_full))
          error (strcat ("GeneralizedLinearModel: X and Y must have the", ...
                         " same number of observations."));
        endif
        if (! isempty (opts.VarNames))
          if (numel (opts.VarNames) != p_raw + 1)
            error ("GeneralizedLinearModel: VarNames must have %d elements.", ...
                   p_raw + 1);
          endif
          pred_names = opts.VarNames(1:p_raw)(:)';
          resp_name  = opts.VarNames{end};
        else
          pred_names = arrayfun (@(k) sprintf ("x%d", k), 1:p_raw, ...
                                 'UniformOutput', false);
          resp_name  = 'y';
        endif
        if (! isempty (opts.ResponseVar))
          resp_name = opts.ResponseVar;
        endif
        var_names_all = [pred_names, {resp_name}];
      else
        tbl           = data;
        X_raw         = [];   # numeric matrix built below via raw_to_codes
        col_names     = tbl.Properties.VariableNames;
        n_total       = height (tbl);
        var_names_all = col_names;
        if (ischar (resp) && ! isempty (resp))
          resp_name = resp;
        elseif (isnumeric (resp) && ! isempty (resp))
          resp_name = 'y';
          if (! isempty (opts.ResponseVar))
            resp_name = opts.ResponseVar;
          endif
          y_ext = double (resp(:));
        elseif (is_formula)
          tparts    = strsplit (modelspec, '~');
          resp_name = strtrim (tparts{1});
        else
          resp_name = col_names{end};
        endif
        if (! isempty (opts.PredictorVars))
          pred_names = opts.PredictorVars;
        elseif (is_formula)
          ## Only the variables the formula names take part in the model; any
          ## other column of the table is carried but is not a predictor.
          used       = formula_var_names (modelspec);
          in_formula = ismember (col_names, used) ...
                       & ! strcmp (col_names, resp_name);
          pred_names = col_names(in_formula);
        else
          pred_names = col_names(! strcmp (col_names, resp_name));
        endif
        p_raw = numel (pred_names);
        if (exist ('y_ext', 'var'))
          y_full = y_ext;
        else
          y_full = double (tbl.(resp_name)(:));
        endif
      endif

      ## Categorical predictor flags.
      cat_logical = false (1, p_raw);
      if (! isempty (opts.CategoricalVars))
        cv = opts.CategoricalVars;
        if (islogical (cv))
          n_cv = min (numel (cv), p_raw);
          cat_logical(1:n_cv) = cv(1:n_cv);
        elseif (isnumeric (cv))
          cat_logical(cv(cv > 0 & cv <= p_raw)) = true;
        elseif (iscell (cv))
          for i = 1:numel (cv)
            cat_logical(strcmp (pred_names, cv{i})) = true;
          endfor
        endif
      endif
      if (istable (data))
        for j = 1:p_raw
          col = tbl.(pred_names{j});
          ## A logical or string column groups its observations just as a cell
          ## or categorical one does, and is coded the same way.
          if (iscell (col) || isa (col, 'categorical') ...
              || islogical (col) || isa (col, 'string'))
            cat_logical(j) = true;
          endif
        endfor
      endif

      ## Missing/excluded masks and the fitting subset.  Only the variables the
      ## model actually uses can make a row missing; an unused table column with
      ## a gap in it does not cost an observation.
      if (! istable (data))
        missing_mask = any (isnan (X_raw), 2) | isnan (y_full);
      else
        used_cols    = [pred_names(:)', {resp_name}];
        used_cols    = used_cols(ismember (used_cols, col_names));
        missing_mask = any (ismissing (tbl(:, used_cols)), 2);
      endif
      missing_mask  = missing_mask(:);
      excluded_mask = false (n_total, 1);
      if (! isempty (opts.Exclude))
        ex = opts.Exclude(:);
        if (islogical (ex))
          excluded_mask(1:numel (ex)) = ex;
        else
          excluded_mask(ex) = true;
        endif
      endif
      avail_mask  = ! missing_mask;
      subset_mask = avail_mask & ! excluded_mask;
      n_obs       = sum (subset_mask);
      if (n_obs < 1)
        error (strcat ("GeneralizedLinearModel: no observations remain after", ...
                       " removing missing/excluded rows."));
      endif

      ## Weights, offset, and binomial trial counts.  Each is kept at full
      ## length so that the per-observation tables can span the input rows.
      has_weights = ! isempty (opts.Weights);
      w_full      = ones (n_total, 1);
      if (has_weights)
        w_full = expand_to_rows (double (opts.Weights(:)), n_total);
      endif
      w_sub = w_full(subset_mask);

      has_offset = ! isempty (opts.Offset);
      off_full   = zeros (n_total, 1);
      if (has_offset)
        off_full = expand_to_rows (double (opts.Offset(:)), n_total);
      endif
      off_sub = off_full(subset_mask);

      distr  = opts.Distribution;
      N_full = [];
      N_sub  = [];
      if (strcmp (distr, 'binomial') && ! isempty (opts.BinomialSize))
        N_full = expand_to_rows (double (opts.BinomialSize(:)), n_total);
        N_sub  = N_full(subset_mask);
      endif

      y_sub = y_full(subset_mask);

      ## Numeric predictor matrix (categorical columns as 1-based codes), used
      ## for encoding and for the slice/effect plots.
      [X_num_full, cat_levels] = raw_to_codes (data, X_raw, tbl, ...
                                   pred_names, cat_logical, n_total);

      ## ------------------------------------------------------------------ ##
      ## Design matrix (intercept included as a column when present).  It is
      ## built over every input row so that the fitted values, residuals, and
      ## diagnostics can be reported for rows kept out of the fit; the fit
      ## itself sees only the subset.
      ## ------------------------------------------------------------------ ##
      if (is_formula)
        if (! istable (data))
          tbl_all = array2table ([X_raw, y_full], ...
                                 'VariableNames', var_names_all);
        else
          tbl_all = tbl;
        endif
        ## Hand the parser only the rows it would keep anyway, so that the
        ## returned design lines up row for row with AVAIL_MASK.
        tbl_avail = tbl_all(avail_mask, :);
        [X_avail, ~, coef_names] = parseWilkinsonFormula ( ...
          modelspec, 'model_matrix', tbl_avail, pred_names(cat_logical));
        coef_names    = coef_names(:)';
        has_intercept = any (strcmp (coef_names, '(Intercept)'));
        enc_names     = coef_names(! strcmp (coef_names, '(Intercept)'));
        X_design_all  = NaN (n_total, columns (X_avail));
        X_design_all(avail_mask, :) = X_avail;
        [terms, cat_info, term_cols] = terms_from_coefnames (coef_names, ...
                              pred_names, cat_logical, data, tbl_avail);
      else
        ## Whether a categorical is given all its indicator columns depends
        ## on the intercept, which has to be settled before encoding.
        [X_enc_all, enc_names, cat_info] = encode_categorical ( ...
          X_num_full, cat_logical, pred_names, cat_levels, ...
          modelspec_has_intercept (modelspec, opts.Intercept));
        [terms, has_intercept, coef_names, emsg] = parse_modelspec ( ...
          modelspec, enc_names, columns (X_enc_all), opts.Intercept);
        if (! isempty (emsg))
          error ("GeneralizedLinearModel: %s", emsg);
        endif
        X_design_all = build_design (terms, X_enc_all);
        term_cols    = [enc_names, {''}];
        ## A missing categorical code encodes as an all-zero dummy row rather
        ## than NaN, so mark the missing rows explicitly.
        X_design_all(missing_mask, :) = NaN;
      endif
      X_design = X_design_all(subset_mask, :);

      ## ------------------------------------------------------------------ ##
      ## Fit the design via glmfit (the intercept is a design column already).
      ## ------------------------------------------------------------------ ##
      if (! isempty (opts.Link))
        linkspec = opts.Link;
      else
        linkspec = default_link_spec (distr);
      endif
      linkname = link_name (linkspec);

      gargs = {'link', linkspec, 'constant', 'off'};
      if (has_weights)
        gargs = [gargs, {'weights', w_sub}];
      endif
      if (has_offset)
        gargs = [gargs, {'offset', off_sub}];
      endif
      if (! isempty (opts.DispersionFlag))
        gargs = [gargs, {'estdisp', ternary(opts.DispersionFlag, 'on', 'off')}];
      endif

      yfit = y_sub;
      if (strcmp (distr, 'binomial') && ! isempty (N_sub))
        yfit = [y_sub .* N_sub, N_sub];
      endif

      [b, dev, stats] = glmfit (X_design, yfit, distr, gargs{:});

      ## ------------------------------------------------------------------ ##
      ## Assemble the object.
      ## ------------------------------------------------------------------ ##
      [flink, dlink, ilink] = getlinkfunctions (linkspec);
      this.Link = struct ('Name', linkname, 'Link', flink, ...
                          'Derivative', dlink, 'Inverse', ilink);
      [devfun, varfun] = glm_family_functions (distr);
      this.Distribution = struct ('Name', distribution_name (distr), ...
                                  'DevianceFunction', devfun, ...
                                  'VarianceFunction', varfun);

      this.Coefficients = table (b(:), stats.se(:), stats.t(:), stats.p(:), ...
        'VariableNames', {'Estimate', 'SE', 'tStat', 'pValue'}, ...
        'RowNames', coef_names(:));

      ## Fitted values over every input row.  Rows kept out of the fit by
      ## 'Exclude' still get a prediction; rows dropped as missing give NaN.
      eta_all     = X_design_all * b + off_full;
      mu_prob_all = ilink (eta_all);   # probability (binomial) or mean (others)
      if (isempty (N_full))
        N_all = ones (n_total, 1);
      else
        N_all = N_full;
      endif
      if (strcmp (distr, 'binomial'))
        mu_resp_all = N_all .* mu_prob_all;
        this.Fitted = table (mu_resp_all, eta_all, mu_prob_all, ...
          'VariableNames', {'Response', 'LinearPredictor', 'Probability'});
      else
        mu_resp_all = mu_prob_all;
        this.Fitted = table (mu_resp_all, eta_all, ...
          'VariableNames', {'Response', 'LinearPredictor'});
      endif

      [r_raw, r_pear, r_ansc, r_dev] = glm_residuals (distr, y_full, ...
                                                      mu_prob_all, N_full);
      r_lin = (y_full - mu_prob_all) .* dlink (mu_prob_all);
      this.Residuals = table (r_raw, r_lin, r_pear, r_ansc, r_dev, ...
        'VariableNames', {'Raw', 'LinearPredictor', 'Pearson', 'Anscombe', ...
                          'Deviance'});

      mu_prob = mu_prob_all(subset_mask);
      mu_resp = mu_resp_all(subset_mask);

      ## Log-likelihood, information criteria, and R-squared measures.  These
      ## use the intercept-only ("null") model as the baseline.
      w_ll = w_sub;
      [bn, nulldev, sn] = glmfit (ones (n_obs, 1), yfit, distr, gargs{:});
      eta0 = bn * ones (n_obs, 1);
      if (! isempty (off_sub))
        eta0 = eta0 + off_sub;
      endif
      mu0 = ilink (eta0);
      LL      = glm_loglik (distr, y_sub, mu_prob, N_sub, w_ll, stats.s);
      LL_null = glm_loglik (distr, y_sub, mu0, N_sub, w_ll, sn.s);
      this.LogLikelihood = LL;
      this.nulldev_ = nulldev;
      this.llnull_  = LL_null;

      k = numel (b);   # MATLAB counts the coefficients only, not the dispersion
      aic = -2 * LL + 2 * k;
      if (n_obs - k - 1 > 0)
        aicc = aic + 2 * k * (k + 1) / (n_obs - k - 1);
      else
        aicc = Inf;
      endif
      this.ModelCriterion = struct ('AIC', aic, ...
        'AICc', aicc, 'BIC', -2 * LL + k * log (n_obs), ...
        'CAIC', -2 * LL + k * (log (n_obs) + 1));

      ybar_w = sum (w_ll .* y_sub) / sum (w_ll);
      sse = sum (w_ll .* (y_sub - mu_resp) .^ 2);
      ssr = sum (w_ll .* (mu_resp - ybar_w) .^ 2);
      sst = sum (w_ll .* (y_sub - ybar_w) .^ 2);
      r2_ord = 1 - sse / max (sst, eps);
      dfe = stats.dfe;
      if (dfe > 0)
        r2_adj = 1 - (1 - r2_ord) * (n_obs - 1) / dfe;
      else
        r2_adj = NaN;
      endif
      r2_dev = 1 - dev / max (nulldev, eps);
      r2_llr = 1 - LL / LL_null;
      r2_gen = 1 - exp (2 * (LL_null - LL) / n_obs);
      r2_gen_max = 1 - exp (2 * LL_null / n_obs);
      this.Rsquared = struct ('Ordinary', r2_ord, 'Adjusted', r2_adj, ...
        'Deviance', r2_dev, 'AdjGeneralized', r2_gen / r2_gen_max, ...
        'LLR', r2_llr);

      ## Model formula.  TERMS is expressed over the encoded design columns; the
      ## formula is expressed over the model's variables, so a categorical's
      ## indicator columns have to be folded back onto the variable they came
      ## from before the term names can be built.
      n_vars  = numel (var_names_all);
      var_idx = zeros (1, p_raw);
      for j = 1:p_raw
        k = find (strcmp (var_names_all, pred_names{j}), 1);
        if (! isempty (k))
          var_idx(j) = k;
        endif
      endfor
      ## Indexed by the columns of TERMS, which are not the coefficient
      ## names: a factor appearing only inside an interaction or a power has a
      ## column without ever being a coefficient.
      [enc2raw, col_pow] = encodednames_to_row (term_cols(1:end-1), pred_names, ...
                                             cat_info);
      terms_var = variable_level_terms (terms(:, 1:end-1), enc2raw, var_idx, ...
                                   n_vars, col_pow);

      this.Formula = LinearFormula (terms_var, var_names_all(:)', ...
                                    'ResponseName', resp_name, ...
                                    'Link', linkspec);
      this.formulastr_ = char (this.Formula);
      in_model         = this.Formula.InModel;

      ## Per-variable information, over the full data rather than the fitting
      ## subset: a range is a property of the variable, not of the fit.
      vi_class   = cell (n_vars, 1);
      vi_range   = cell (n_vars, 1);
      vi_inmodel = in_model(:);
      vi_iscat   = false (n_vars, 1);
      vi_iscat(var_idx(var_idx > 0)) = cat_logical(var_idx > 0);
      for j = 1:n_vars
        if (istable (data))
          col = tbl.(var_names_all{j});
        elseif (j <= p_raw)
          col = X_raw(:, j);
        else
          col = y_full;
        endif
        [vi_class{j}, vi_range{j}] = variable_class_and_range (col);
      endfor
      this.VariableInfo = table (vi_class, vi_range, vi_inmodel, vi_iscat, ...
        'VariableNames', {'Class', 'Range', 'InModel', 'IsCategorical'}, ...
        'RowNames', var_names_all(:));

      if (istable (data))
        this.Variables = tbl;
      else
        this.Variables = array2table ([X_raw, y_full], ...
                                      'VariableNames', var_names_all);
      endif
      this.ObservationInfo = table (w_full, excluded_mask, missing_mask, ...
        subset_mask, ...
        'VariableNames', {'Weights', 'Excluded', 'Missing', 'Subset'});

      ## Leverage, hat matrix, and Cook's distance.  The hat matrix is the
      ## asymmetric IRLS form H = X inv(X'WX) X' W, whose diagonal is the
      ## leverage; rows outside the fit contribute nothing to it.
      dmu_deta = 1 ./ dlink (mu_prob);
      vwt      = glm_varfun (distr, mu_prob, N_sub);
      w_irls   = (dmu_deta .^ 2) ./ max (vwt, realmin) .* w_ll;
      XtW      = (X_design .* w_irls)';
      Hs       = X_design * ((XtW * X_design) \ XtW);
      lev      = diag (Hs);
      pear     = r_pear(subset_mask);
      phi_c    = ternary (logical (stats.estdisp), stats.s, 1);
      cooksd   = w_ll .* pear .^ 2 .* lev ...
                 ./ (max (1 - lev, eps) .^ 2 * numel (b) * phi_c);

      lev_all = zeros (n_total, 1);
      lev_all(subset_mask) = lev;
      cd_all = NaN (n_total, 1);
      cd_all(subset_mask) = cooksd;
      hat_all = zeros (n_total, n_total);
      hat_all(subset_mask, subset_mask) = Hs;
      this.Diagnostics = table (lev_all, cd_all, hat_all, ...
        'VariableNames', {'Leverage', 'CooksDistance', 'HatMatrix'});

      this.b_          = b;
      this.stats_      = stats;
      this.distr_      = distr;
      this.linkarg_    = linkspec;
      this.binomsize_  = opts.BinomialSize;
      this.terms_      = terms;
      this.catinfo_    = cat_info;
      this.encnames_   = enc_names;
      this.prednames_  = pred_names;
      this.design_     = X_design;
      this.subset_     = subset_mask;
      this.xmeans_     = mean (X_num_full(subset_mask,:), 1);

      this.CoefficientNames         = coef_names;
      this.CoefficientCovariance    = stats.covb;
      this.NumCoefficients          = numel (b);
      this.NumEstimatedCoefficients = numel (b);
      this.NumPredictors            = p_raw;
      this.NumObservations          = n_obs;
      this.NumVariables             = n_vars;
      this.Deviance                 = dev;
      this.DFE                      = stats.dfe;
      this.Dispersion               = stats.s;
      this.DispersionEstimated      = logical (stats.estdisp);
      this.SSE                      = sse;
      this.SSR                      = ssr;
      this.SST                      = sst;
      this.Offset                   = off_full;
      this.LikelihoodPenalty        = string ("none");
      this.ResponseName             = resp_name;
      this.PredictorNames           = pred_names(:);
      this.VariableNames            = var_names_all(:);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{yhat} =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {GeneralizedLinearModel} {[@var{yhat}, @var{yci}] =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {GeneralizedLinearModel} {[@dots{}] =} predict (@dots{}, @var{Name}, @var{Value})
    ##
    ## Predict the response of the model @var{mdl} at the new predictor data
    ## @var{Xnew} (a numeric matrix or a table).  Predictions are on the mean
    ## (response) scale.  With two outputs, @var{yci} is an @math{m}-by-2 matrix
    ## of confidence intervals.  The @qcode{'Alpha'} pair sets the confidence
    ## level to @math{100 (1 - @var{Alpha})%} (default 0.05).
    ##
    ## @end deftypefn
    function [yhat, yci] = predict (mdl, Xnew, varargin)
      if (nargin < 2)
        error ("GeneralizedLinearModel.predict: Xnew is required.");
      endif
      alpha  = 0.05;
      offnew = [];
      for k = 1:2:numel (varargin)
        switch (lower (varargin{k}))
          case 'alpha';   alpha  = varargin{k+1};
          case 'offset';  offnew = varargin{k+1};
          otherwise
            error (strcat ("GeneralizedLinearModel.predict: unknown", ...
                           " parameter '%s'."), varargin{k});
        endswitch
      endfor

      ## Marshal Xnew into a raw numeric predictor matrix (encode categoricals).
      p_raw = mdl.NumPredictors;
      if (istable (Xnew))
        n_new = height (Xnew);
        X_raw = zeros (n_new, p_raw);
        for j = 1:p_raw
          col = Xnew.(mdl.prednames_{j});
          cidx = [];
          if (! isempty (mdl.catinfo_.names))
            cidx = find (strcmp (mdl.catinfo_.names, mdl.prednames_{j}));
          endif
          if (! isempty (cidx) && iscell (col))
            levels_j = mdl.catinfo_.levels{cidx};
            codes    = zeros (n_new, 1);
            for L = 1:numel (levels_j)
              codes(strcmp (col, levels_j{L})) = L;
            endfor
            X_raw(:, j) = codes;
          else
            X_raw(:, j) = double (col);
          endif
        endfor
      else
        X_raw = double (Xnew);
        if (columns (X_raw) != p_raw)
          error ("GeneralizedLinearModel.predict: Xnew must have %d columns.", ...
                 p_raw);
        endif
      endif

      X_enc    = reencode_predictors (X_raw, mdl.prednames_, mdl.catinfo_, ...
                                      mdl.encnames_);
      X_design = build_design (mdl.terms_, X_enc);

      valargs = {'constant', 'off', 'confidence', 1 - alpha};
      if (! isempty (offnew))
        valargs = [valargs, {'offset', offnew(:)}];
      endif
      if (! isempty (mdl.binomsize_))
        N = mdl.binomsize_(:);
        if (isscalar (N))
          N = N * ones (rows (X_design), 1);
        endif
        valargs = [valargs, {'size', N}];
      endif

      if (nargout > 1)
        [yhat, ylo, yhi] = glmval (mdl.b_, X_design, mdl.linkarg_, ...
                                   mdl.stats_, valargs{:});
        yci = [yhat - ylo, yhat + yhi];
      else
        yhat = glmval (mdl.b_, X_design, mdl.linkarg_, valargs{:});
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{yhat} =} feval (@var{mdl}, @var{x1}, @var{x2}, @dots{})
    ##
    ## Predict the response by passing each predictor as a separate argument (a
    ## scalar or column vector), returning point predictions on the mean scale.
    ## Equivalent to @code{predict (@var{mdl}, [@var{x1}, @var{x2}, @dots{}])}.
    ##
    ## @end deftypefn
    function yhat = feval (mdl, varargin)
      if (numel (varargin) == 1 && size (varargin{1}, 2) == mdl.NumPredictors)
        Xnew = varargin{1};
      else
        cols = cellfun (@(c) c(:), varargin, 'UniformOutput', false);
        Xnew = [cols{:}];
      endif
      yhat = predict (mdl, Xnew);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{ci} =} coefCI (@var{mdl})
    ## @deftypefnx {GeneralizedLinearModel} {@var{ci} =} coefCI (@var{mdl}, @var{alpha})
    ##
    ## Confidence intervals for the coefficient estimates.  @var{ci} is a
    ## @math{k}-by-2 matrix of @math{100 (1 - @var{alpha})%} intervals (default
    ## @var{alpha} = 0.05).  The @math{t} distribution is used when the
    ## dispersion was estimated, the normal distribution otherwise.
    ##
    ## @end deftypefn
    function ci = coefCI (mdl, alpha)
      if (nargin < 2)
        alpha = 0.05;
      endif
      if (! (isscalar (alpha) && isnumeric (alpha) && alpha >= 0 && alpha <= 1))
        error ("GeneralizedLinearModel.coefCI: ALPHA must be in [0, 1].");
      endif
      b  = mdl.Coefficients.Estimate;
      se = mdl.Coefficients.SE;
      crit = tinv (1 - alpha / 2, mdl.DFE);
      ci = [b - crit .* se, b + crit .* se];
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{p} =} coefTest (@var{mdl})
    ## @deftypefnx {GeneralizedLinearModel} {[@var{p}, @var{stat}, @var{df}] =} coefTest (@var{mdl}, @var{H})
    ##
    ## Wald test of the linear hypothesis @math{H b = 0} on the coefficients.
    ## @var{H} is an @math{m}-by-@math{k} contrast matrix; when omitted it tests
    ## that all coefficients except the intercept are zero (the model versus the
    ## constant model).  Returns the p-value @var{p}, and optionally the test
    ## statistic @var{stat} and its numerator degrees of freedom @var{df}.
    ##
    ## @end deftypefn
    function [p, stat, df] = coefTest (mdl, H)
      k = mdl.NumCoefficients;
      if (nargin < 2)
        ipos = find (strcmp (mdl.CoefficientNames, '(Intercept)'));
        rows_h = setdiff (1:k, ipos);
        H = zeros (numel (rows_h), k);
        for i = 1:numel (rows_h)
          H(i, rows_h(i)) = 1;
        endfor
      endif
      b  = mdl.Coefficients.Estimate;
      df = rows (H);
      Hb = H * b;
      HVH = H * mdl.CoefficientCovariance * H';
      ## Wald F statistic (MATLAB uses the F distribution for coefTest).
      stat = (Hb' * (HVH \ Hb)) / df;
      p = 1 - fcdf (stat, df, mdl.DFE);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{tbl} =} devianceTest (@var{mdl})
    ##
    ## Likelihood-ratio (deviance) test of the fitted model against the
    ## intercept-only model.  Returns a table with the deviance, degrees of
    ## freedom, and p-value of each model, the last row giving the chi-square
    ## statistic (the drop in deviance) and its p-value.
    ##
    ## @end deftypefn
    function tbl = devianceTest (mdl)
      dev_full = mdl.Deviance;
      dev_null = mdl.nulldev_;
      df_diff  = mdl.NumCoefficients - 1;
      chi2stat = dev_null - dev_full;
      pval     = 1 - chi2cdf (chi2stat, df_diff);
      Deviance   = [dev_null; dev_full];
      DFE        = [mdl.DFE + df_diff; mdl.DFE];
      chi2Stat   = [NaN; chi2stat];
      pValue     = [NaN; pval];
      tbl = table (Deviance, DFE, chi2Stat, pValue, ...
        'VariableNames', {'Deviance', 'DFE', 'chi2Stat', 'pValue'}, ...
        'RowNames', {'(Constant)', mdl.formulastr_});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{ysim} =} random (@var{mdl})
    ## @deftypefnx {GeneralizedLinearModel} {@var{ysim} =} random (@var{mdl}, @var{Xnew})
    ##
    ## Simulate responses from the fitted model.  With one argument the fitted
    ## values are used; otherwise the mean is predicted at the new predictor
    ## data @var{Xnew}.  A random draw from the response distribution about that
    ## mean is returned.
    ##
    ## @end deftypefn
    function ysim = random (mdl, Xnew)
      if (nargin < 2)
        mu = mdl.Fitted.Response;
      else
        mu = predict (mdl, Xnew);
      endif
      switch (mdl.distr_)
        case 'normal'
          ysim = normrnd (mu, sqrt (mdl.Dispersion));
        case 'poisson'
          ysim = poissrnd (mu);
        case 'binomial'
          if (isempty (mdl.binomsize_))
            ysim = binornd (1, mu);
          else
            N = mdl.binomsize_(:);
            if (isscalar (N))
              N = N * ones (size (mu));
            endif
            ysim = binornd (N, mu ./ N);
          endif
        case 'gamma'
          ysim = gamrnd (1 / mdl.Dispersion, mu * mdl.Dispersion);
        case 'inverse gaussian'
          error (strcat ("GeneralizedLinearModel.random: simulation for the", ...
                         " 'inverse gaussian' distribution is not supported."));
      endswitch
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{h} =} plotResiduals (@var{mdl})
    ## @deftypefnx {GeneralizedLinearModel} {@var{h} =} plotResiduals (@var{mdl}, @var{plottype})
    ## @deftypefnx {GeneralizedLinearModel} {@var{h} =} plotResiduals (@dots{}, @qcode{'ResidualType'}, @var{rt})
    ##
    ## Plot the model residuals.  @var{plottype} is one of @qcode{'histogram'}
    ## (default), @qcode{'caseorder'}, @qcode{'fitted'}, @qcode{'lagged'}, or
    ## @qcode{'probability'}.  @qcode{'ResidualType'} picks the residual column
    ## (@qcode{'Raw'} default, @qcode{'Pearson'}, @qcode{'Deviance'},
    ## @qcode{'Anscombe'}).  Returns the graphics handle.
    ##
    ## @end deftypefn
    function h = plotResiduals (mdl, varargin)
      ptype = 'histogram';
      rtype = 'Raw';
      k = 1;
      if (numel (varargin) >= 1 && ischar (varargin{1}) ...
          && ! strcmpi (varargin{1}, 'ResidualType'))
        ptype = varargin{1};
        k = 2;
      endif
      for i = k:2:numel (varargin)
        if (strcmpi (varargin{i}, 'ResidualType'))
          rtype = varargin{i+1};
        endif
      endfor
      r = mdl.Residuals.(rtype);
      switch (lower (ptype))
        case 'histogram'
          h = hist (r);
          xlabel ("Residuals");  ylabel ("Frequency");
        case 'caseorder'
          h = plot (1:numel (r), r, 'x');
          xlabel ("Row number");  ylabel (sprintf ("%s residuals", rtype));
        case 'fitted'
          h = plot (mdl.Fitted.Response, r, 'x');
          xlabel ("Fitted values");  ylabel (sprintf ("%s residuals", rtype));
        case 'lagged'
          h = plot (r(1:end-1), r(2:end), 'x');
          xlabel ("Residual (t-1)");  ylabel ("Residual (t)");
        case 'probability'
          [rs, idx] = sort (r);
          n = numel (rs);
          pp = ((1:n)' - 0.5) / n;
          h = plot (rs, norminv (pp), 'x');
          xlabel (sprintf ("%s residuals", rtype));
          ylabel ("Standard normal quantiles");
        otherwise
          error ("GeneralizedLinearModel.plotResiduals: bad plot type '%s'.", ...
                 ptype);
      endswitch
      title (sprintf ("Residuals: %s", ptype));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{h} =} plotDiagnostics (@var{mdl})
    ## @deftypefnx {GeneralizedLinearModel} {@var{h} =} plotDiagnostics (@var{mdl}, @var{plottype})
    ##
    ## Plot observation diagnostics.  @var{plottype} is @qcode{'leverage'}
    ## (default) or @qcode{'cookd'} (Cook's distance).  A reference line marks
    ## the usual threshold.  Returns the graphics handle.
    ##
    ## @end deftypefn
    function h = plotDiagnostics (mdl, plottype)
      if (nargin < 2)
        plottype = 'leverage';
      endif
      n = mdl.NumObservations;
      idx = find (mdl.subset_);
      switch (lower (plottype))
        case 'leverage'
          lev = mdl.Diagnostics.Leverage(idx);
          h = stem (idx, lev, 'Marker', 'x');
          ref = 2 * mdl.NumCoefficients / n;
          ylabel ("Leverage");
        case 'cookd'
          cd = mdl.Diagnostics.CooksDistance(idx);
          h = stem (idx, cd, 'Marker', 'x');
          ref = 3 * mean (cd);
          ylabel ("Cook's distance");
        otherwise
          error (strcat ("GeneralizedLinearModel.plotDiagnostics: bad plot", ...
                         " type '%s'."), plottype);
      endswitch
      hold on;
      plot ([idx(1), idx(end)], [ref, ref], 'r--');
      hold off;
      xlabel ("Row number");
      title (sprintf ("Diagnostics: %s", plottype));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{h} =} plotEffects (@var{mdl})
    ##
    ## Main-effects plot: for each predictor, the change in the fitted mean
    ## response as that predictor sweeps its observed range while the others are
    ## held at their means.  Returns the graphics handle.
    ##
    ## @end deftypefn
    function h = plotEffects (mdl)
      p    = mdl.NumPredictors;
      base = mdl.xmeans_;
      eff  = zeros (p, 1);
      for j = 1:p
        a = base;  b = base;
        a(j) = base(j) - 1;
        b(j) = base(j) + 1;
        eff(j) = predict (mdl, b) - predict (mdl, a);
      endfor
      h = plot (eff, 1:p, 'o', 'MarkerFaceColor', 'b');
      hold on;
      for j = 1:p
        plot ([0, eff(j)], [j, j], 'b-');
      endfor
      plot ([0, 0], [0.5, p + 0.5], 'k:');
      hold off;
      ylim ([0.5, p + 0.5]);
      set (gca, 'YTick', 1:p, 'YTickLabel', mdl.PredictorNames);
      xlabel ("Effect on fitted response (two-unit change)");
      title ("Main effects");
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{h} =} plotAdjustedResponse (@var{mdl}, @var{var})
    ##
    ## Adjusted-response plot for the predictor @var{var} (a name or index): the
    ## fitted mean response as @var{var} sweeps its observed range with other
    ## predictors held at their means, overlaid on the partial residuals.
    ## Returns the graphics handle.
    ##
    ## @end deftypefn
    function h = plotAdjustedResponse (mdl, var)
      j = resolve_predictor (mdl, var);
      xj = linspace (mdl.xmeans_(j) - 2, mdl.xmeans_(j) + 2, 50)';
      Xg = repmat (mdl.xmeans_, numel (xj), 1);
      Xg(:,j) = xj;
      yg = predict (mdl, Xg);
      h = plot (xj, yg, 'b-', 'LineWidth', 1.5);
      xlabel (mdl.PredictorNames{j});
      ylabel (sprintf ("Adjusted %s", mdl.ResponseName));
      title (sprintf ("Adjusted response for %s", mdl.PredictorNames{j}));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{h} =} plotAdded (@var{mdl}, @var{var})
    ##
    ## Added-variable (partial-regression) plot for the predictor @var{var} (a
    ## name or index): the response residuals from the model without @var{var}
    ## against the residuals of @var{var} regressed on the remaining predictors.
    ## Returns the graphics handle.
    ##
    ## @end deftypefn
    function h = plotAdded (mdl, var)
      j = resolve_predictor (mdl, var);
      X = mdl.design_;
      ## Locate the design column of the requested predictor.
      cidx = find (strcmp (mdl.CoefficientNames, mdl.PredictorNames{j}), 1);
      if (isempty (cidx))
        error (strcat ("GeneralizedLinearModel.plotAdded: predictor is not a", ...
                       " single design column."));
      endif
      others = setdiff (1:columns (X), cidx);
      Xo = X(:,others);
      ## The design holds the fitted rows only, so take the residuals of the
      ## same rows.
      ry = mdl.Residuals.Raw(mdl.subset_);
      xj = X(:,cidx);
      bx = Xo \ xj;
      rx = xj - Xo * bx;
      h = plot (rx, ry, 'x');
      hold on;
      bb = rx \ ry;
      xr = [min(rx), max(rx)];
      plot (xr, bb * xr, 'r-');
      hold off;
      xlabel (sprintf ("%s (adjusted)", mdl.PredictorNames{j}));
      ylabel (sprintf ("%s (adjusted)", mdl.ResponseName));
      title (sprintf ("Added variable plot for %s", mdl.PredictorNames{j}));
    endfunction

  endmethods

endclassdef

## Parse GeneralizedLinearModel name/value options into a structure.
function opts = glm_parse_nv (nv)
  opts = struct ('Distribution', 'normal', 'Link', [], 'Weights', [], ...
                 'Offset', [], 'BinomialSize', [], 'Intercept', true, ...
                 'DispersionFlag', [], 'VarNames', {{}}, 'ResponseVar', '', ...
                 'PredictorVars', {{}}, 'CategoricalVars', [], 'Exclude', []);
  if (mod (numel (nv), 2) != 0)
    error (strcat ("GeneralizedLinearModel: optional arguments must be", ...
                   " Name-Value pairs."));
  endif
  for k = 1:2:numel (nv)
    switch (lower (nv{k}))
      case 'distribution';   opts.Distribution   = lower (nv{k+1});
      case 'link';           opts.Link           = nv{k+1};
      case 'weights';        opts.Weights        = nv{k+1};
      case 'offset';         opts.Offset         = nv{k+1};
      case 'binomialsize';   opts.BinomialSize   = nv{k+1};
      case 'intercept';      opts.Intercept      = logical (nv{k+1});
      case 'dispersionflag'; opts.DispersionFlag = logical (nv{k+1});
      case 'varnames';       opts.VarNames       = nv{k+1};
      case 'responsevar';    opts.ResponseVar    = nv{k+1};
      case 'predictorvars';  opts.PredictorVars  = nv{k+1};
      case 'categoricalvars'; opts.CategoricalVars = nv{k+1};
      case 'exclude';        opts.Exclude        = nv{k+1};
      otherwise
        error (strcat ("GeneralizedLinearModel: unknown parameter", ...
                       " name '%s'."), nv{k});
    endswitch
  endfor
  if (! any (strcmp (opts.Distribution, {'normal', 'binomial', 'poisson', ...
                                         'gamma', 'inverse gaussian'})))
    error ("GeneralizedLinearModel: unknown distribution '%s'.", ...
           opts.Distribution);
  endif
endfunction

## Convert raw predictor data (numeric matrix or table) to numeric level codes
## for categorical columns; return the codes and the per-column level labels.
function [X_num, cat_levels] = raw_to_codes (data, X_raw, tbl, pred_names, ...
                                             cat_logical, n_total)
  p = numel (pred_names);
  X_num      = zeros (n_total, p);
  cat_levels = cell (1, p);
  for j = 1:p
    if (istable (data))
      col = tbl.(pred_names{j});
      if (iscell (col))
        ## Appearance order, so that the omitted reference level is the one the
        ## data shows first; see parseWilkinsonFormula for the formula path.
        [cat_levels{j}, ~, ic] = unique (col, 'stable');
        X_num(:, j) = ic;
      elseif (isa (col, 'categorical'))
        cat_levels{j} = categories (col);
        [~, ic] = ismember (cellstr (col), cat_levels{j});
        X_num(:, j) = ic;
      else
        X_num(:, j) = double (col(:));
        cat_levels{j} = {};
      endif
    else
      if (cat_logical(j))
        uvals = sort (unique (X_raw(isfinite (X_raw(:,j)), j)));
        cat_levels{j} = strtrim (cellstr (num2str (uvals(:))));
        [~, ic] = ismember (X_raw(:,j), uvals);
        X_num(:, j) = ic;
      else
        X_num(:, j) = X_raw(:, j);
        cat_levels{j} = {};
      endif
    endif
  endfor
endfunction

## Broadcast a scalar to N rows; otherwise pass the vector through.
function v = expand_to_rows (v, n)
  if (isscalar (v))
    v = v * ones (n, 1);
  endif
endfunction

## Base names of the variables a Wilkinson formula refers to (response
## included), with any '^k' power suffix stripped.
function vnames = formula_var_names (modelspec)
  schema = parseWilkinsonFormula (modelspec, 'matrix');
  vnames = regexprep (schema.VariableNames, '\^\d+$', '');
  vnames = unique (vnames, 'stable');
endfunction

## Display name of a response distribution.
function name = distribution_name (distr)
  switch (distr)
    case 'normal';            name = 'Normal';
    case 'binomial';          name = 'Binomial';
    case 'poisson';           name = 'Poisson';
    case 'gamma';             name = 'Gamma';
    case 'inverse gaussian';  name = 'Inverse Gaussian';
  endswitch
endfunction

## Deviance and variance functions of a response distribution, as the handles
## reported by the Distribution property.
function [devfun, varfun] = glm_family_functions (distr)
  switch (distr)
    case 'normal'
      devfun = @(mu, y) (y - mu) .^ 2;
      varfun = @(mu) ones (size (mu));
    case 'binomial'
      devfun = @(mu, y, N) 2 * N .* (y .* log ((y + (y == 0)) ./ mu) ...
                           + (1 - y) .* log ((1 - y + (y == 1)) ./ (1 - mu)));
      varfun = @(mu, N) mu .* (1 - mu) ./ N;
    case 'poisson'
      devfun = @(mu, y) 2 * (y .* (log ((y + (y == 0)) ./ mu)) - (y - mu));
      varfun = @(mu) mu;
    case 'gamma'
      devfun = @(mu, y) 2 * (-log (y ./ mu) + (y - mu) ./ mu);
      varfun = @(mu) mu .^ 2;
    case 'inverse gaussian'
      devfun = @(mu, y) (((y - mu) ./ mu) .^ 2) ./ y;
      varfun = @(mu) mu .^ 3;
  endswitch
endfunction

## Raw, Pearson, Anscombe, and deviance residuals of a GLM fit, evaluated for
## every input row.  Y is the response (a proportion for the binomial family),
## MU the fitted probability/mean, and N the binomial trial counts (empty
## otherwise).  Matches the residuals @code{glmfit} returns for the fitted rows.
function [raw, pear, ansc, devr] = glm_residuals (distr, y, mu, N)

  [devfun, varfun] = glm_family_functions (distr);
  if (strcmp (distr, 'binomial'))
    if (isempty (N))
      N = ones (size (y));
    endif
    raw  = (y - mu) .* N;
    sd   = sqrt (varfun (mu, N));
    devn = devfun (mu, y, N);
  else
    raw  = y - mu;
    sd   = sqrt (varfun (mu));
    devn = devfun (mu, y);
  endif
  pear = (y - mu) ./ (sd + (y == mu));
  devr = sign (y - mu) .* sqrt (max (0, devn));

  switch (distr)
    case 'normal'
      ansc = y - mu;
    case 'binomial'
      ab   = 2 / 3;
      ansc = beta (ab, ab) * (betainc (y, ab, ab) - betainc (mu, ab, ab)) ...
             ./ ((mu .* (1 - mu)) .^ (1 / 6) ./ sqrt (N));
    case 'poisson'
      ansc = 1.5 * ((y .^ (2 / 3) - mu .^ (2 / 3)) ./ mu .^ (1 / 6));
    case 'gamma'
      pwr  = 1 / 3;
      ansc = 3 * (y .^ pwr - mu .^ pwr) ./ mu .^ pwr;
    case 'inverse gaussian'
      ansc = (log (y) - log (mu)) ./ mu;
  endswitch

endfunction

## Log-likelihood of a GLM fit.  Y is the response (proportion for binomial),
## MU the fitted probability/mean, N the binomial trials (empty otherwise), W
## the prior weights, and PHI the dispersion parameter.
function ll = glm_loglik (distr, y, mu, N, w, phi)
  rmin = realmin;
  switch (distr)
    case 'normal'
      ne  = sum (w);
      s2  = sum (w .* (y - mu) .^ 2) / ne;
      ll  = -0.5 * ne * (log (2 * pi * s2) + 1);
    case 'poisson'
      ll = sum (w .* (y .* log (max (mu, rmin)) - mu - gammaln (y + 1)));
    case 'binomial'
      if (isempty (N))
        N = ones (size (y));
      endif
      yc = y .* N;
      ll = sum (w .* (gammaln (N + 1) - gammaln (yc + 1) ...
           - gammaln (N - yc + 1) + yc .* log (max (mu, rmin)) ...
           + (N - yc) .* log (max (1 - mu, rmin))));
    case 'gamma'
      a  = 1 ./ phi;
      ll = sum (w .* (a .* log (a) - a .* log (mu) + (a - 1) .* log (y) ...
           - a .* y ./ mu - gammaln (a)));
    case 'inverse gaussian'
      ll = sum (w .* (-0.5 * (log (2 * pi * phi .* y .^ 3) ...
           + (y - mu) .^ 2 ./ (phi .* mu .^ 2 .* y))));
  endswitch
endfunction

## Resolve a predictor reference (name or index) to a column index.
function j = resolve_predictor (mdl, var)
  if (ischar (var))
    j = find (strcmp (mdl.PredictorNames, var), 1);
    if (isempty (j))
      error ("GeneralizedLinearModel: unknown predictor '%s'.", var);
    endif
  else
    j = var;
    if (! (isscalar (j) && j >= 1 && j <= mdl.NumPredictors))
      error ("GeneralizedLinearModel: predictor index out of range.");
    endif
  endif
endfunction

## Variance function V(mu) of the response distribution (for the IRLS weights
## used in leverage/Cook's-distance diagnostics).
function v = glm_varfun (distr, mu, N)
  [~, varfun] = glm_family_functions (distr);
  if (strcmp (distr, 'binomial'))
    if (isempty (N))
      N = ones (size (mu));
    endif
    v = varfun (mu, N);
  else
    v = varfun (mu);
  endif
endfunction

## Response distribution's canonical link specification.
function spec = default_link_spec (distr)
  switch (distr)
    case 'normal';            spec = 'identity';
    case 'binomial';          spec = 'logit';
    case 'poisson';           spec = 'log';
    case 'gamma';             spec = 'reciprocal';
    case 'inverse gaussian';  spec = -2;
  endswitch
endfunction

## Human-readable name of a link specification.
function name = link_name (linkarg)
  if (ischar (linkarg))
    name = linkarg;
  elseif (isnumeric (linkarg) && isscalar (linkarg))
    name = sprintf ('%g', linkarg);
  else
    name = 'custom';
  endif
endfunction

## Small inline conditional helper.
function out = ternary (cond, a, b)
  if (cond)
    out = a;
  else
    out = b;
  endif
endfunction


%!demo
%! ## Fit a Poisson regression and inspect the model object.
%! X = [0.1, 1.2; 0.4, 0.7; 1.1, 0.2; 1.5, 1.9; 0.3, 0.5; 1.8, 1.1; 0.9, 0.3];
%! y = [1; 0; 2; 3; 1; 4; 2];
%! mdl = fitglm (X, y, 'Distribution', 'poisson');
%! disp (mdl.Coefficients)
%! printf ("Deviance = %g,  AIC = %g\n", mdl.Deviance, mdl.ModelCriterion.AIC);

## Test direct construction and input validation
%!test
%! X = [1 2; 2 1; 3 4; 4 3; 5 6; 6 5; 1 3; 4 2];
%! y = [1; 0; 2; 3; 2; 4; 1; 3];
%! mdl = GeneralizedLinearModel (X, y, "linear", "Distribution", "poisson");
%! assert_equal (class (mdl), "GeneralizedLinearModel");
%! assert_equal (mdl.Distribution.Name, "Poisson");
%! assert_equal (mdl.Link.Name, "log");
%!error<DATA, RESP, and MODELSPEC are required> GeneralizedLinearModel (1)
%!error<X must be a real matrix.> ...
%! GeneralizedLinearModel ("a", [1;2], "linear")
%!error<Y must be a real vector.> ...
%! GeneralizedLinearModel ([1 2; 3 4], "a", "linear")
%!error<unknown distribution 'wibble'.> ...
%! GeneralizedLinearModel ([1 2; 3 4], [1;0], "linear", ...
%!                         "Distribution", "wibble")
%!error<GeneralizedLinearModel: \(\) indexing is not supported> ...
%! mdl = GeneralizedLinearModel ([1 2; 2 1; 3 4; 4 3; 5 6; 6 5], ...
%!                               [1;0;2;3;2;4], "linear", ...
%!                               "Distribution", "poisson"); mdl(1);

## Comprehensive property and method coverage
%!shared X, yp, yb, yn
%! X = [ 0.37,  0.06,  1.76; -0.76, -1.52,  0.84;  0.76, -0.19, -0.47; ...
%!      -0.80, -2.74, -0.90;  0.08,  0.39,  1.05; -0.41, -0.03,  0.74; ...
%!       0.23,  1.21,  0.35;  0.66,  0.94,  0.13;  0.66, -0.12, -0.06; ...
%!       2.09,  1.33, -0.71;  1.50,  0.08, -0.52;  0.59,  0.07, -1.13; ...
%!      -1.17, -0.35, -1.28;  0.68,  0.63, -0.80; -0.69,  0.08,  0.41; ...
%!       2.04,  0.96, -0.56];
%! yp = [5 2 0 3 1 1 0 1 2 1 3 0 0 1 1 3]';
%! yb = [1 1 1 0 0 1 1 1 1 1 1 0 0 0 0 1]';
%! yn = [2.1 -0.3 1.2 -1.1 0.8 0.4 1.5 1.1 0.6 2.9 2.0 0.7 -1.3 0.9 -0.2 2.5]';

%!test  # a normal-distribution GLM with identity link reproduces OLS exactly
%! mdl = fitglm (X, yn, "Distribution", "normal");
%! b_ols = [ones(16, 1), X] \ yn;
%! assert_equal (mdl.Coefficients.Estimate, b_ols, 1e-10);
%! assert_equal (mdl.Fitted.Response, [ones(16, 1), X] * b_ols, 1e-10);
%! assert_equal (mdl.Residuals.Raw, yn - [ones(16, 1), X] * b_ols, 1e-10);
%! assert_equal (mdl.Deviance, sum ((yn - [ones(16, 1), X] * b_ols) .^ 2), 1e-10);
%! assert_equal (mdl.Link.Name, "identity");

%!test  # the normal-GLM standard errors match those from regress
%! mdl = fitglm (X, yn, "Distribution", "normal");
%! [~, bint] = regress (yn, [ones(16, 1), X], 0.05);
%! se_reg = (bint(:,1) - bint(:,2)) / 2 / tinv (0.025, 12);
%! assert_equal (mdl.Coefficients.SE, se_reg, 1e-9);

%!test  # normal-GLM dispersion is estimated as SSE/DFE; log-likelihood closed form
%! mdl = fitglm (X, yn, "Distribution", "normal");
%! rss = mdl.Deviance;
%! assert_equal (mdl.DispersionEstimated, true);
%! assert_equal (mdl.Dispersion, rss / mdl.DFE, 1e-12);
%! assert_equal (mdl.LogLikelihood, -8 * (log (2 * pi * rss / 16) + 1), 1e-6);

%!test  # Poisson coefficients and fit statistics (verified against MATLAB)
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [-0.3420955; 1.2804868; -1.0743272; 0.8395779], 1e-6);
%! assert_equal (mdl.Deviance, 7.403008, 1e-5);
%! assert_equal (mdl.LogLikelihood, -18.543280, 1e-5);
%! assert_equal (mdl.ModelCriterion.AIC, 45.086559, 1e-5);
%! assert_equal (mdl.Rsquared.Deviance, 0.6627677, 1e-6);

%!test  # scalar count/size properties of the Poisson fit
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumEstimatedCoefficients, 4);
%! assert_equal (mdl.NumPredictors, 3);
%! assert_equal (mdl.NumObservations, 16);
%! assert_equal (mdl.DFE, 12);
%! assert_equal (mdl.ResponseName, "y");
%! assert_equal (mdl.CoefficientNames, {'(Intercept)', 'x1', 'x2', 'x3'});

%!test  # Poisson has a fixed unit dispersion (not estimated)
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.Dispersion, 1);
%! assert_equal (mdl.DispersionEstimated, false);
%! assert_equal (mdl.Distribution.Name, "Poisson");
%! assert_equal (mdl.Link.Name, "log");

%!test  # the coefficient covariance is symmetric with SE^2 on its diagonal
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! C = mdl.CoefficientCovariance;
%! assert_equal (size (C), [4, 4]);
%! assert_equal (C, C', 1e-14);
%! assert_equal (diag (C), mdl.Coefficients.SE .^ 2, 1e-12);

%!test  # predict at the training data reproduces the fitted response
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! [yhat, yci] = predict (mdl, X);
%! assert_equal (yhat, mdl.Fitted.Response, 1e-10);
%! assert_equal (size (yci), [16, 2]);
%! assert_equal (all (yci(:,1) <= yhat & yhat <= yci(:,2)), true);

%!test  # feval evaluates the model and agrees with predict
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (feval (mdl, X(:,1), X(:,2), X(:,3)), predict (mdl, X), 1e-12);

%!test  # coefCI matches the t-interval and honours a custom alpha
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! b = mdl.Coefficients.Estimate;  se = mdl.Coefficients.SE;
%! t95 = tinv (0.975, mdl.DFE);
%! assert_equal (coefCI (mdl), [b - t95 * se, b + t95 * se], 1e-12);
%! t90 = tinv (0.95, mdl.DFE);
%! assert_equal (coefCI (mdl, 0.10), [b - t90 * se, b + t90 * se], 1e-12);

%!test  # coefTest gives the Wald F statistic against the constant model
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! [p, F, df] = coefTest (mdl);
%! assert_equal (F, 3.685312, 1e-5);
%! assert_equal (p, 0.04331745, 1e-7);
%! assert_equal (df, 3);

%!test  # devianceTest chi-square equals the drop from the null deviance
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! dt = devianceTest (mdl);
%! assert_equal (class (dt), "table");
%! assert_equal (dt.chi2Stat(2), dt.Deviance(1) - dt.Deviance(2), 1e-10);

%!test  # information criteria satisfy their defining identities
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! k = mdl.NumEstimatedCoefficients;  ll = mdl.LogLikelihood;
%! assert_equal (mdl.ModelCriterion.AIC, -2 * ll + 2 * k, 1e-9);
%! assert_equal (mdl.ModelCriterion.BIC, -2 * ll + k * log (16), 1e-9);

%!test  # raw residuals are response minus fit; random draws match the response size
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.Residuals.Raw, yp - mdl.Fitted.Response, 1e-12);
%! ysim = random (mdl);
%! assert_equal (size (ysim), [16, 1]);
%! assert_equal (all (ysim == round (ysim) & ysim >= 0), true);

%!test  # binomial/logistic fit: coefficients agree with the glmfit engine
%! mdl = fitglm (X, yb, "Distribution", "binomial");
%! assert_equal (mdl.Coefficients.Estimate, glmfit (X, yb, "binomial"), 1e-8);
%! assert_equal (mdl.Link.Name, "logit");
%! assert_equal (all (mdl.Fitted.Response >= 0 & mdl.Fitted.Response <= 1), true);
%! assert_equal (mdl.Deviance, 10.997099, 1e-5);

%!test  # an interaction model adds the cross term and one coefficient
%! mdl = fitglm (X, yp, "interactions", "Distribution", "poisson");
%! assert_equal (any (strcmp (mdl.CoefficientNames, "x1:x2")), true);
%! assert_equal (mdl.NumCoefficients, 7);

%!test  # an offset is stored and applied
%! mdl = fitglm (X, yp, "Distribution", "poisson", "Offset", log (2 * ones (16, 1)));
%! assert_equal (numel (mdl.Offset), 16);

%!test  # disp prints the model header and the coefficient table
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! s = evalc ("disp (mdl)");
%! assert_equal (isempty (strfind (s, "Generalized linear regression model")), false);
%! assert_equal (isempty (strfind (s, "Estimate")), false);

%!test  # chained subsref reaches property -> table column -> element
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (numel (mdl.Coefficients.Estimate), 4);
%! assert_equal (mdl.Coefficients.Estimate(1), -0.3420955, 1e-6);
%! assert_equal (mdl.Coefficients.Estimate(2), mdl.Coefficients{2, "Estimate"}, 1e-12);

%!test  # the diagnostic and effect plots run without error
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! hf = figure ("visible", "off");
%! unwind_protect
%!   plotResiduals (mdl);
%!   plotResiduals (mdl, "fitted", "ResidualType", "Pearson");
%!   plotDiagnostics (mdl);
%!   plotDiagnostics (mdl, "cookd");
%!   plotEffects (mdl);
%!   plotAdjustedResponse (mdl, 1);
%!   plotAdded (mdl, "x2");
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## The property surface below is checked against MATLAB R2024a on the same data.

%!test  # each family reports its MATLAB display name and its own functions
%! names = {'Normal', 'Binomial', 'Poisson', 'Gamma', 'Inverse Gaussian'};
%! dists = {'normal', 'binomial', 'poisson', 'gamma', 'inverse gaussian'};
%! resp  = {yn, yb, yp, abs(yn) + 1, abs(yn) + 1};
%! for k = 1:numel (dists)
%!   mdl = fitglm (X, resp{k}, "Distribution", dists{k});
%!   assert_equal (mdl.Distribution.Name, names{k});
%!   dev = mdl.Distribution.DevianceFunction;
%!   var = mdl.Distribution.VarianceFunction;
%!   assert_equal (is_function_handle (dev), true);
%!   assert_equal (is_function_handle (var), true);
%! endfor

%!test  # the variance function is the family's variance, not its scale
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! v = mdl.Distribution.VarianceFunction;
%! assert_equal (v ([1; 4; 9]), [1; 4; 9]);
%! mdl = fitglm (X, abs (yn) + 1, "Distribution", "gamma");
%! v = mdl.Distribution.VarianceFunction;
%! assert_equal (v ([1; 2; 3]), [1; 4; 9]);

%!test  # the deviance function reproduces the reported deviance
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! d = mdl.Distribution.DevianceFunction;
%! assert_equal (sum (d (mdl.Fitted.Response, yp)), mdl.Deviance, 1e-10);

%!test  # sums of squares match MATLAB and reproduce Rsquared.Ordinary
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.SSE, 5.878748587925331, 1e-12);
%! assert_equal (mdl.SSR, 22.693546014008586, 1e-11);
%! assert_equal (mdl.SST, 30, 1e-12);
%! assert_equal (mdl.Rsquared.Ordinary, 1 - mdl.SSE / mdl.SST, 1e-14);

%!test  # SSE + SSR closes to SST only for the identity link
%! mdl = fitglm (X, yn);
%! assert_equal (mdl.SSE, 1.007522513007451, 1e-12);
%! assert_equal (mdl.SSR, 20.549977486992557, 1e-11);
%! assert_equal (mdl.SST, 21.557500000000005, 1e-12);
%! assert_equal (mdl.SSE + mdl.SSR, mdl.SST, 1e-12);
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.SSE + mdl.SSR < mdl.SST, true);

%!test  # binomial and gamma sums of squares
%! mdl = fitglm (X, yb, "Distribution", "binomial");
%! assert_equal ([mdl.SSE, mdl.SSR, mdl.SST], ...
%!   [1.529775270428043, 2.016466154627350, 3.75], 1e-11);
%! mdl = fitglm (X, abs (yn) + 1, "Distribution", "gamma");
%! assert_equal ([mdl.SSE, mdl.SSR, mdl.SST], ...
%!   [3.537122914639216, 5.173700107311616, 9.45], 1e-9);

%!test  # counts, penalty, and observation names
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.NumVariables, 4);
%! assert_equal (char (mdl.LikelihoodPenalty), "none");
%! assert_equal (mdl.ObservationNames, {});
%! assert_equal (size (mdl.Steps), [0, 0]);

%!test  # Variables holds the data the model was built from
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (class (mdl.Variables), "table");
%! assert_equal (size (mdl.Variables), [16, 4]);
%! assert_equal (mdl.Variables.Properties.VariableNames, ...
%!               {'x1', 'x2', 'x3', 'y'});
%! assert_equal (mdl.Variables{:, 'x2'}, X(:,2));
%! assert_equal (mdl.Variables{:, 'y'}, yp);

%!test  # ObservationInfo spans the input rows and records why each was used
%! mdl = fitglm (X, yn);
%! assert_equal (size (mdl.ObservationInfo), [16, 4]);
%! assert_equal (mdl.ObservationInfo.Properties.VariableNames, ...
%!               {'Weights', 'Excluded', 'Missing', 'Subset'});
%! assert_equal (mdl.ObservationInfo.Weights, ones (16, 1));
%! assert_equal (any (mdl.ObservationInfo.Excluded), false);
%! assert_equal (all (mdl.ObservationInfo.Subset), true);

%!test  # excluded rows keep their weight and drop out of the fit
%! mdl = fitglm (X, yn, "Exclude", [2 5], "Weights", (1:16)' / 16);
%! assert_equal (mdl.NumObservations, 14);
%! assert_equal (mdl.DFE, 10);
%! assert_equal (mdl.ObservationInfo.Weights, (1:16)' / 16);
%! assert_equal (find (mdl.ObservationInfo.Excluded), [2; 5]);
%! assert_equal (any (mdl.ObservationInfo.Missing), false);
%! assert_equal (find (! mdl.ObservationInfo.Subset), [2; 5]);
%! assert_equal (mdl.SSE, 0.376329223083291, 1e-11);
%! assert_equal (mdl.SST, 12.408120155038763, 1e-10);
%! assert_equal (mdl.Dispersion, 0.037632922308329, 1e-12);

%!test  # an excluded row still gets a fitted value; a missing row does not
%! mdl = fitglm (X, yn, "Exclude", [2 5], "Weights", (1:16)' / 16);
%! assert_equal (size (mdl.Fitted), [16, 2]);
%! assert_equal (mdl.Fitted.Response(1:3), ...
%!   [1.452778236704902; -0.409857651323262; 1.021905387197724], 1e-11);
%! assert_equal (mdl.Residuals.Raw(2), 0.109857651323262, 1e-11);
%! Xm = X;  Xm(3,2) = NaN;
%! mdl = fitglm (Xm, yn);
%! assert_equal (mdl.NumObservations, 15);
%! assert_equal (find (mdl.ObservationInfo.Missing), 3);
%! assert_equal (isnan (mdl.Fitted.Response(3)), true);
%! assert_equal (isnan (mdl.Residuals.Raw(3)), true);
%! assert_equal (mdl.Fitted.Response(1), 1.675087141815260, 1e-11);
%! assert_equal (mdl.SSE, 0.990693655817047, 1e-11);

%!test  # the binomial family adds a Probability column to Fitted
%! mdl = fitglm (X, yb, "Distribution", "binomial");
%! assert_equal (mdl.Fitted.Properties.VariableNames, ...
%!               {'Response', 'LinearPredictor', 'Probability'});
%! assert_equal (mdl.Fitted.Probability(1), 0.998111791212415, 1e-9);
%! assert_equal (mdl.Fitted.Response, mdl.Fitted.Probability, 1e-14);
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.Fitted.Properties.VariableNames, ...
%!               {'Response', 'LinearPredictor'});

%!test  # with BinomialSize the response is a proportion and the fit a count
%! N = 5 * ones (16, 1);
%! y = [3 4 5 1 0 4 3 5 4 5 5 1 0 2 1 5]' ./ N;
%! mdl = fitglm (X, y, "Distribution", "binomial", "BinomialSize", N);
%! assert_equal (mdl.Fitted.Response, N .* mdl.Fitted.Probability, 1e-12);
%! assert_equal (mdl.Residuals.Raw, (y - mdl.Fitted.Probability) .* N, 1e-12);

%!test  # residuals gain the linear-predictor column, in MATLAB's order
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.Residuals.Properties.VariableNames, ...
%!   {'Raw', 'LinearPredictor', 'Pearson', 'Anscombe', 'Deviance'});
%! assert_equal (mdl.Residuals.LinearPredictor(1), 0.066685156329760, 1e-12);
%! ## For the log link the working residual is the raw one over the mean.
%! assert_equal (mdl.Residuals.LinearPredictor, ...
%!               mdl.Residuals.Raw ./ mdl.Fitted.Response, 1e-12);

%!test  # the identity link leaves the linear-predictor residual raw
%! mdl = fitglm (X, yn);
%! assert_equal (mdl.Residuals.LinearPredictor, mdl.Residuals.Raw, 1e-14);
%! mdl = fitglm (X, yb, "Distribution", "binomial");
%! assert_equal (mdl.Residuals.LinearPredictor(1), 1.001891780864838, 1e-7);

%!test  # leverage, Cook's distance, and the hat matrix
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.Diagnostics.Properties.VariableNames, ...
%!               {'Leverage', 'CooksDistance', 'HatMatrix'});
%! assert_equal (mdl.Diagnostics.Leverage(1), 0.768034312650850, 1e-7);
%! assert_equal (mdl.Diagnostics.CooksDistance(1), 0.074381550609974, 1e-7);
%! assert_equal (size (mdl.Diagnostics.HatMatrix), [16, 16]);
%! assert_equal (diag (mdl.Diagnostics.HatMatrix), mdl.Diagnostics.Leverage, ...
%!               1e-14);
%! assert_equal (sum (mdl.Diagnostics.Leverage), mdl.NumCoefficients, 1e-8);

%!test  # the GLM hat matrix is the asymmetric IRLS form, so H(i,j) != H(j,i)
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! H = mdl.Diagnostics.HatMatrix;
%! assert_equal (H(1,2), 0.198729950284820, 1e-7);
%! assert_equal (H(2,1), 0.334913252272241, 1e-7);
%! ## The identity link with unit weights makes it symmetric again.
%! H = fitglm (X, yn).Diagnostics.HatMatrix;
%! assert_equal (H, H', 1e-12);

%!test  # normal-family diagnostics match MATLAB
%! mdl = fitglm (X, yn);
%! assert_equal (mdl.Diagnostics.Leverage(1:3), ...
%!   [0.392783227704423; 0.308688028422869; 0.103894702431155], 1e-12);
%! assert_equal (mdl.Diagnostics.CooksDistance(1), 0.562165139840681, 1e-11);

%!test  # rows outside the fit have no leverage and no Cook's distance
%! mdl = fitglm (X, yn, "Exclude", [2 5], "Weights", (1:16)' / 16);
%! assert_equal (mdl.Diagnostics.Leverage([2, 5]), [0; 0]);
%! assert_equal (isnan (mdl.Diagnostics.CooksDistance([2, 5])), [true; true]);
%! assert_equal (mdl.Diagnostics.Leverage(1), 0.116937877466335, 1e-11);
%! assert_equal (mdl.Diagnostics.CooksDistance(1), 0.026081406006346, 1e-10);
%! assert_equal (mdl.Diagnostics.HatMatrix([2, 5], :), zeros (2, 16));

%!test  # VariableInfo carries the class and range of every variable
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.VariableInfo.Properties.VariableNames, ...
%!               {'Class', 'Range', 'InModel', 'IsCategorical'});
%! assert_equal (size (mdl.VariableInfo), [4, 4]);
%! assert_equal (mdl.VariableInfo.Range{1}, [-1.17, 2.09], 1e-14);
%! assert_equal (mdl.VariableInfo.Range{2}, [-2.74, 1.33], 1e-14);
%! assert_equal (mdl.VariableInfo.Range{4}, [0, 5]);
%! assert_equal (mdl.VariableInfo.InModel, [true; true; true; false]);
%! assert_equal (mdl.VariableInfo.IsCategorical, false (4, 1));
%! assert_equal (mdl.VariableInfo.Class, repmat ({'double'}, 4, 1));

%!test  # a range spans the whole variable, not just the fitted rows
%! mdl = fitglm (X, yn, "Exclude", [1, 10]);
%! assert_equal (mdl.VariableInfo.Range{1}, [-1.17, 2.09], 1e-14);
%! assert_equal (mdl.VariableInfo.Range{4}, [-1.3, 2.9], 1e-14);

%!test  # the offset spans the input rows and is zero when none was given
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! assert_equal (mdl.Offset, zeros (16, 1));
%! mdl = fitglm (X, yp, "Distribution", "poisson", ...
%!               "Offset", log (2 * ones (16, 1)));
%! assert_equal (mdl.Offset, log (2) * ones (16, 1), 1e-14);

%!test  # the formula is a structure describing the model over its variables
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! f = mdl.Formula;
%! assert_equal (f.ResponseName, "y");
%! assert_equal (f.LinearPredictor, "1 + x1 + x2 + x3");
%! assert_equal (f.PredictorNames, {'x1', 'x2', 'x3'});
%! assert_equal (f.VariableNames, {'x1', 'x2', 'x3', 'y'});
%! assert_equal (f.TermNames, {'(Intercept)'; 'x1'; 'x2'; 'x3'});
%! assert_equal (f.Terms, [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0]);
%! assert_equal (f.Link, "log");
%! assert_equal (f.InModel, [true, true, true, false]);
%! assert_equal (f.HasIntercept, true);
%! assert_equal ([f.NTerms, f.NVars, f.NPredictors], [4, 4, 3]);
%! assert_equal (f.FunctionCalls, cell (1, 0));
%! assert_equal (f.ModelFun ([1; 2], [1, 3]), 7);

%!test  # a term whose parts are all present is written as a product
%! f = fitglm (X, yn, "interactions").Formula;
%! assert_equal (f.LinearPredictor, "1 + x1*x2 + x1*x3 + x2*x3");
%! assert_equal (f.NTerms, 7);
%! assert_equal (f.TermNames, ...
%!   {'(Intercept)'; 'x1'; 'x2'; 'x3'; 'x1:x2'; 'x1:x3'; 'x2:x3'});
%! f = fitglm (X(:,1:2), yn, "quadratic").Formula;
%! assert_equal (f.LinearPredictor, "1 + x1*x2 + x1^2 + x2^2");
%! assert_equal (f.TermNames, ...
%!   {'(Intercept)'; 'x1'; 'x2'; 'x1:x2'; 'x1^2'; 'x2^2'});
%! assert_equal (f.Terms, [0 0 0; 1 0 0; 0 1 0; 1 1 0; 2 0 0; 0 2 0]);

%!test  # only pairs collapse into a product, never a three-way interaction
%! f = fitglm (X, yn, "full").Formula;
%! assert_equal (f.LinearPredictor, "1 + x1*x2 + x1*x3 + x2*x3 + x1:x2:x3");

%!test  # dropping the intercept drops the leading 1
%! f = fitglm (X, yn, "Intercept", false).Formula;
%! assert_equal (f.LinearPredictor, "x1 + x2 + x3");
%! assert_equal (f.HasIntercept, false);
%! assert_equal (f.Terms, [1 0 0 0; 0 1 0 0; 0 0 1 0]);

%!test  # disp names the link alongside the response
%! s = evalc ("disp (fitglm (X, yp, 'Distribution', 'poisson'))");
%! assert_equal (isempty (strfind (s, "log(y) ~ 1 + x1 + x2 + x3")), false);
%! s = evalc ("disp (fitglm (X, yb, 'Distribution', 'binomial'))");
%! assert_equal (isempty (strfind (s, "logit(y) ~ 1 + x1 + x2 + x3")), false);
%! s = evalc ("disp (fitglm (X, yn))");
%! assert_equal (isempty (strfind (s, "y ~ 1 + x1 + x2 + x3")), false);

%!test  # a numeric link is named by its exponent and shown as a power
%! mdl = fitglm (X, abs (yn) + 1, "Distribution", "inverse gaussian");
%! assert_equal (mdl.Link.Name, "-2");
%! s = evalc ("disp (mdl)");
%! assert_equal (isempty (strfind (s, "power(y,-2) ~ 1 + x1")), false);

%!test  # a table model keeps every column but only models what it names
%! tbl = array2table ([X, yn], "VariableNames", {'a', 'b', 'c', 'resp'});
%! mdl = fitglm (tbl, "resp ~ 1 + a + b");
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.NumVariables, 4);
%! assert_equal (mdl.PredictorNames, {'a'; 'b'});
%! assert_equal (mdl.VariableNames, {'a'; 'b'; 'c'; 'resp'});
%! assert_equal (mdl.Formula.LinearPredictor, "1 + a + b");
%! assert_equal (mdl.Formula.InModel, [true, true, false, false]);
%! assert_equal (mdl.VariableInfo.InModel, [true; true; false; false]);
%! assert_equal (size (mdl.Variables), [16, 4]);

%!test  # a column the model does not use cannot cost an observation
%! tbl = array2table ([X, yn], "VariableNames", {'a', 'b', 'c', 'resp'});
%! tbl.c(4) = NaN;
%! mdl = fitglm (tbl, "resp ~ 1 + a + b");
%! assert_equal (mdl.NumObservations, 16);
%! assert_equal (any (mdl.ObservationInfo.Missing), false);

%!test  # a grouping column is reported by its own class and its levels
%! g = {'lo'; 'hi'; 'lo'; 'hi'; 'lo'; 'hi'; 'lo'; 'hi'; ...
%!      'lo'; 'hi'; 'lo'; 'hi'; 'lo'; 'hi'; 'lo'; 'hi'};
%! tbl = table (X(:,1), g, yn, "VariableNames", {'u', 'grp', 'resp'});
%! mdl = fitglm (tbl, "resp ~ 1 + u + grp");
%! assert_equal (mdl.VariableInfo.Class, {'double'; 'cell'; 'double'});
%! assert_equal (mdl.VariableInfo.IsCategorical, [false; true; false]);
%! ## The levels are listed as the design codes them: first seen first.
%! assert_equal (mdl.VariableInfo.Range{2}, {'lo', 'hi'});
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.NumVariables, 3);
%! ## The indicator columns fold back onto the variable they came from.
%! assert_equal (mdl.Formula.NTerms, 3);
%! assert_equal (sort (mdl.Formula.TermNames), {'(Intercept)'; 'grp'; 'u'});

%!test  # devianceTest labels the alternative with the rendered formula
%! mdl = fitglm (X, yp, "Distribution", "poisson");
%! dt = devianceTest (mdl);
%! assert_equal (dt.Properties.RowNames, ...
%!               {'(Constant)'; 'log(y) ~ 1 + x1 + x2 + x3'});

%!test  # a power keeps its exponent in the variable's own column
%! u = [1;2;3;4;5;6;7;8;9;10;11;12];
%! v = [2;1;4;3;6;5;8;7;10;9;12;11];
%! cnt = [2;5;3;7;4;6;8;3;5;9;4;6];
%! tbl = table (u, v, cnt);
%! f = fitglm (tbl, 'cnt ~ 1 + u*v + u^2', 'Distribution', 'poisson').Formula;
%! assert_equal (char (f), 'log(cnt) ~ 1 + u*v + u^2');
%! assert_equal (f.Terms, [0 0 0; 1 0 0; 0 1 0; 1 1 0; 2 0 0]);
%! assert_equal (f.TermNames, {'(Intercept)'; 'u'; 'v'; 'u:v'; 'u^2'});

%!test  # an interaction survives when one factor is not a main effect
%! u = [1;2;3;4;5;6;7;8;9;10;11;12];
%! v = [2;1;4;3;6;5;8;7;10;9;12;11];
%! cnt = [2;5;3;7;4;6;8;3;5;9;4;6];
%! tbl = table (u, v, cnt);
%! f = fitglm (tbl, 'cnt ~ 1 + u + u:v', 'Distribution', 'poisson').Formula;
%! assert_equal (char (f), 'log(cnt) ~ 1 + u + u:v');
%! assert_equal (f.Terms, [0 0 0; 1 0 0; 1 1 0]);

%!test  # dropping the intercept codes the first categorical in full, both paths
%! u = [1;2;3;4;5;6;7;8;9;10;11;12];
%! h2 = {'b';'c';'a';'b';'c';'a';'b';'c';'a';'b';'c';'a'};
%! cnt = [2;5;3;7;4;6;8;3;5;9;4;6];
%! tbl = table (u, h2, cnt);
%! m = fitglm (tbl, 'cnt ~ h2 - 1', 'Distribution', 'poisson');
%! assert_equal (m.CoefficientNames, {'h2_b', 'h2_c', 'h2_a'});
%! m = fitglm (tbl, 'linear', 'Distribution', 'poisson', 'Intercept', false);
%! assert_equal (m.CoefficientNames, {'u', 'h2_b', 'h2_c', 'h2_a'});
