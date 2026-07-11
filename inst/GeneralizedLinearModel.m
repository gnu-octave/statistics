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
## @code{Dispersion}, @code{Residuals}, @code{Fitted}, @code{Distribution}, and
## @code{Link}.  Fitted models support the @code{predict} and @code{feval}
## methods for prediction.
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

    ## Structure of fitted values (Response on the mean scale, LinearPredictor).
    Fitted = [];

    ## Table of residuals (Raw, Pearson, Deviance, Anscombe).
    Residuals = [];

    ## Offset vector added to the linear predictor (empty if none).
    Offset = [];

    ## Name of the response variable.
    ResponseName = 'y';

    ## Cell array of predictor variable names.
    PredictorNames = {};

    ## Cell array of all variable names (predictors then response).
    VariableNames = {};

    ## Table of per-variable information (class, range, in-model, categorical).
    VariableInfo = [];

    ## Character-vector representation of the model formula.
    Formula = '';

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
      if (! isempty (this.Formula))
        fprintf ("      %s\n", this.Formula);
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
          if (iscell (col) || isa (col, 'categorical'))
            cat_logical(j) = true;
          endif
        endfor
      endif

      ## Missing/excluded masks and the fitting subset.
      if (! istable (data))
        missing_mask = any (isnan (X_raw), 2) | isnan (y_full);
      else
        missing_mask = any (ismissing (tbl), 2);
      endif
      excluded_mask = false (n_total, 1);
      if (! isempty (opts.Exclude))
        ex = opts.Exclude(:);
        if (islogical (ex))
          excluded_mask(1:numel (ex)) = ex;
        else
          excluded_mask(ex) = true;
        endif
      endif
      subset_mask = ! missing_mask & ! excluded_mask;
      n_obs       = sum (subset_mask);
      if (n_obs < 1)
        error (strcat ("GeneralizedLinearModel: no observations remain after", ...
                       " removing missing/excluded rows."));
      endif

      ## Weights, offset, and binomial trial counts, subset to the fitting rows.
      if (isempty (opts.Weights))
        w_sub = [];
      else
        w_full = double (opts.Weights(:));
        w_sub  = w_full(subset_mask);
      endif
      off_sub = [];
      if (! isempty (opts.Offset))
        off_full = double (opts.Offset(:));
        off_sub  = off_full(subset_mask);
      endif
      distr = opts.Distribution;
      N_sub = [];
      if (strcmp (distr, 'binomial') && ! isempty (opts.BinomialSize))
        N = opts.BinomialSize(:);
        if (isscalar (N))
          N = N * ones (n_total, 1);
        endif
        N_sub = N(subset_mask);
      endif

      y_sub = y_full(subset_mask);

      ## ------------------------------------------------------------------ ##
      ## Design matrix (intercept included as a column when present).
      ## ------------------------------------------------------------------ ##
      if (is_formula)
        if (! istable (data))
          tbl_sub = array2table ([X_raw(subset_mask,:), y_sub], ...
                                 'VariableNames', var_names_all);
        else
          tbl_sub = tbl(subset_mask, :);
        endif
        [X_design, ~, coef_names] = parseWilkinsonFormula ( ...
          modelspec, 'model_matrix', tbl_sub);
        coef_names    = coef_names(:)';
        has_intercept = any (strcmp (coef_names, '(Intercept)'));
        enc_names     = coef_names(! strcmp (coef_names, '(Intercept)'));
        [terms, cat_info] = terms_from_coefnames (coef_names, pred_names, ...
                                                  cat_logical, data, tbl_sub);
      else
        [X_num_full, cat_levels] = raw_to_codes (data, X_raw, tbl, ...
                                     pred_names, cat_logical, n_total);
        X_num_sub = X_num_full(subset_mask, :);
        [X_enc_sub, enc_names, cat_info] = encode_categorical ( ...
          X_num_sub, cat_logical, pred_names, cat_levels);
        [terms, has_intercept, coef_names, emsg] = parse_modelspec ( ...
          modelspec, enc_names, columns (X_enc_sub), opts.Intercept);
        if (! isempty (emsg))
          error ("GeneralizedLinearModel: %s", emsg);
        endif
        X_design = build_design (terms, X_enc_sub);
      endif

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
      if (! isempty (w_sub))
        gargs = [gargs, {'weights', w_sub}];
      endif
      if (! isempty (off_sub))
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
      this.Distribution = struct ('Name', distr, ...
                                  'DevianceFunction', [], ...
                                  'VarianceFunction', []);

      this.Coefficients = table (b(:), stats.se(:), stats.t(:), stats.p(:), ...
        'VariableNames', {'Estimate', 'SE', 'tStat', 'pValue'}, ...
        'RowNames', coef_names(:));
      this.Residuals = table (stats.resid(:), stats.residp(:), ...
        stats.residd(:), stats.resida(:), ...
        'VariableNames', {'Raw', 'Pearson', 'Deviance', 'Anscombe'});

      eta = X_design * b;
      if (! isempty (off_sub))
        eta = eta + off_sub;
      endif
      mu = ilink (eta);
      if (strcmp (distr, 'binomial') && ! isempty (N_sub))
        mu = N_sub .* mu;
      endif
      this.Fitted = struct ('Response', mu, 'LinearPredictor', eta);

      ## Model formula string.
      if (is_formula)
        this.Formula = modelspec;
      else
        this.Formula = formula_string (resp_name, coef_names, has_intercept);
      endif

      ## VariableInfo (predictors and response).
      vi_iscat = [cat_logical, false](:);
      vi_class = repmat ({'double'}, numel (var_names_all), 1);
      vi_class(vi_iscat) = {'categorical'};
      vi_inmodel = [true(1, p_raw), false](:);
      this.VariableInfo = table (vi_class, vi_iscat, vi_inmodel, ...
        'VariableNames', {'Class', 'IsCategorical', 'InModel'}, ...
        'RowNames', var_names_all(:));

      this.b_          = b;
      this.stats_      = stats;
      this.distr_      = distr;
      this.linkarg_    = linkspec;
      this.binomsize_  = opts.BinomialSize;
      this.terms_      = terms;
      this.catinfo_    = cat_info;
      this.encnames_   = enc_names;
      this.prednames_  = pred_names;

      this.CoefficientNames         = coef_names;
      this.CoefficientCovariance    = stats.covb;
      this.NumCoefficients          = numel (b);
      this.NumEstimatedCoefficients = numel (b);
      this.NumPredictors            = p_raw;
      this.NumObservations          = n_obs;
      this.Deviance                 = dev;
      this.DFE                      = stats.dfe;
      this.Dispersion               = stats.s;
      this.DispersionEstimated      = logical (stats.estdisp);
      this.Offset                   = off_sub;
      this.ResponseName             = resp_name;
      this.PredictorNames           = pred_names;
      this.VariableNames            = var_names_all;

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
        [cat_levels{j}, ~, ic] = unique (col);
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
        cat_levels{j} = cellstr (num2str (uvals(:)));
        [~, ic] = ismember (X_raw(:,j), uvals);
        X_num(:, j) = ic;
      else
        X_num(:, j) = X_raw(:, j);
        cat_levels{j} = {};
      endif
    endif
  endfor
endfunction

## Reconstruct a terms matrix and categorical info from formula coefficient
## names (PATH A), mirroring the encoding used for the design.
function [terms, cat_info] = terms_from_coefnames (coef_names, pred_names, ...
                                                   cat_logical, data, tbl_sub)
  cat_info.names  = {};
  cat_info.levels = {};
  if (istable (data))
    for j = 1:numel (pred_names)
      if (cat_logical(j))
        col = tbl_sub.(pred_names{j});
        if (iscell (col))
          levels_j = unique (col);
        elseif (isa (col, 'categorical'))
          levels_j = categories (col);
        else
          levels_j = {};
        endif
        cat_info.names{end+1}  = pred_names{j};
        cat_info.levels{end+1} = levels_j;
      endif
    endfor
  endif

  nc = numel (coef_names);
  atomic = {};
  for t = 1:nc
    if (strcmp (coef_names{t}, '(Intercept)'))
      continue;
    endif
    for f = strsplit (coef_names{t}, ':')
      if (! any (strcmp (atomic, f{1})))
        atomic{end+1} = f{1};
      endif
    endfor
  endfor
  terms = zeros (nc, numel (atomic) + 1);
  for t = 1:nc
    if (strcmp (coef_names{t}, '(Intercept)'))
      continue;
    endif
    for f = strsplit (coef_names{t}, ':')
      terms(t, strcmp (atomic, f{1})) = 1;
    endfor
  endfor
endfunction

## Build a Wilkinson-style formula string from coefficient names.
function s = formula_string (resp_name, coef_names, has_intercept)
  terms = coef_names(! strcmp (coef_names, '(Intercept)'));
  if (has_intercept)
    rhs = strjoin ([{'1'}, terms], ' + ');
  elseif (isempty (terms))
    rhs = '1';
  else
    rhs = strjoin (terms, ' + ');
  endif
  s = sprintf ('%s ~ %s', resp_name, rhs);
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
    name = sprintf ('power(%g)', linkarg);
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
