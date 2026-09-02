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
## @deftp {statistics} {} CoxModel
##
## Cox proportional hazards regression model class.
##
## A @code{CoxModel} object encapsulates a Cox proportional hazards model of a
## survival time on one or more predictors, fitted by maximizing the Cox partial
## likelihood.  It is the object counterpart of @code{coxphfit} and is normally
## created with the @code{fitcox} function.
##
## The model states that an observation with predictor values @var{x} has hazard
##
## @tex
## $$ h(x, t) = h_0(t)\exp\left(\sum_{j=1}^{p} x_{j} b_j\right) $$
## @end tex
## @ifnottex
## @math{h(x, t) = h_0(t) exp (x' b)}
## @end ifnottex
##
## where @math{h_0(t)} is an unspecified baseline hazard.  The model carries no
## constant term: any constant is absorbed into that baseline.
##
## The most useful properties are @code{Coefficients} (a table of estimates,
## standard errors, @math{z}-statistics and p-values), @code{Hazard} (the
## estimated baseline cumulative hazard), @code{LogLikelihood},
## @code{Residuals}, and the three p-values @code{LikelihoodRatioTestPValue},
## @code{ProportionalHazardsPValue} and
## @code{ProportionalHazardsPValueGlobal}.  Fitted models support the
## @code{survival}, @code{hazardratio}, @code{coefci}, @code{linhyptest},
## @code{plotSurvival} and @code{discardResiduals} methods.
##
## A categorical predictor expands to indicator columns, one per level bar the
## first, which the baseline hazard carries; the indicator columns are named
## @qcode{@var{name}_@var{level}} and enter the default baseline as zero, while
## a numeric predictor enters it as its mean.
##
## @code{ProportionalHazardsPValue} is a Grambsch-Therneau test of each
## coefficient against the mid-ranks of the event times, and
## @code{ProportionalHazardsPValueGlobal} the same test taken over the whole
## model.  A small p-value is evidence that the hazard ratio moves with time,
## which is what proportionality denies.
##
## @strong{Deviations from MATLAB, all in naming.}  MATLAB derives the names
## reported by a fitted model from three different places and they need not
## agree with one another: with default predictor names its @code{Formula}
## reads @qcode{'y ~ x1 + x2'} in lower case while @code{PredictorNames} holds
## @qcode{'X1'} and @qcode{'X2'}, and supplying @qcode{'PredictorNames'}
## changes @code{ResponseName} from @qcode{'y'} to the name of the variable
## passed as the response.  Here the names are consistent by construction:
## @code{ResponseName} is @qcode{'y'} unless the data came from a table, the
## @code{Formula} is built from @code{PredictorNames} and
## @code{ResponseName}, and neither depends on which optional arguments were
## given.  Every fitted quantity agrees with MATLAB.
##
## @seealso{fitcox, coxphfit, GeneralizedLinearModel, LinearModel}
## @end deftp

classdef CoxModel

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} Coefficients
    ##
    ## Coefficient estimates and their statistics
    ##
    ## A table with one row per encoded predictor column, its row names the
    ## encoded column names, and the variables @qcode{Beta}, @qcode{SE},
    ## @qcode{zStat} and @qcode{pValue}.  This property is read-only.
    ##
    ## @end deftp
    Coefficients = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer counting the predictor variables of the model,
    ## before any categorical predictor is encoded.  This property is
    ## read-only.
    ##
    ## @end deftp
    NumPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} LogLikelihood
    ##
    ## Log-likelihood of the fitted model
    ##
    ## A scalar, the maximised Cox partial log-likelihood.  This property is
    ## read-only.
    ##
    ## @end deftp
    LogLikelihood = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} Hazard
    ##
    ## Estimated baseline cumulative hazard
    ##
    ## A numeric matrix of event times in its first column and the cumulative
    ## hazard at them in its second.  A stratified model adds a third column
    ## holding the stratum each row belongs to.  This property is read-only.
    ##
    ## @end deftp
    Hazard = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} PredictorNames
    ##
    ## Names of the predictor variables
    ##
    ## A cell array of character vectors with one name per predictor.  This
    ## property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector naming the response, which is the survival time.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName = 'y';

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} Formula
    ##
    ## Model formula
    ##
    ## A @code{LinearFormula} object describing the terms of the model.  This
    ## property is read-only.
    ##
    ## @end deftp
    Formula = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} Baseline
    ##
    ## Predictor values the baseline hazard is evaluated at
    ##
    ## The baseline the fit used, one row per stratum.  It is reported as it
    ## was given when @code{fitcox} was given one, a scalar staying a scalar,
    ## and otherwise holds the rows the fit was centred on.  This property is
    ## read-only.
    ##
    ## @end deftp
    Baseline = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} Stratification
    ##
    ## Stratification levels used in the fit
    ##
    ## The distinct levels of the stratification variable.  It is empty when
    ## the model is not stratified.  This property is read-only.
    ##
    ## @end deftp
    Stratification = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} CoefficientCovariance
    ##
    ## Estimated covariance of the coefficients
    ##
    ## A square numeric matrix, one row and column per encoded predictor
    ## column, holding the estimated covariance of the estimates in
    ## @qcode{Coefficients}.  This property is read-only.
    ##
    ## @end deftp
    CoefficientCovariance = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} StandardError
    ##
    ## Standard errors of the coefficients
    ##
    ## A numeric column vector, the square roots of the diagonal of
    ## @qcode{CoefficientCovariance}, which is column @qcode{SE} of
    ## @qcode{Coefficients}.  This property is read-only.
    ##
    ## @end deftp
    StandardError = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} Residuals
    ##
    ## Residuals of the fitted model
    ##
    ## A table with one row per observation and the variables
    ## @qcode{CoxSnell}, @qcode{Deviance}, @qcode{Martingale},
    ## @qcode{Schoenfeld}, @qcode{ScaledSchoenfeld}, @qcode{Score} and
    ## @qcode{ScaledScore}.  This property is read-only.
    ##
    ## @end deftp
    Residuals = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} ProportionalHazardsPValue
    ##
    ## Proportional hazards test, one predictor at a time
    ##
    ## A numeric vector with one @math{p} value per predictor, testing whether
    ## that predictor's effect is constant over time.  A small value is
    ## evidence against the proportional hazards assumption.  This property is
    ## read-only.
    ##
    ## @end deftp
    ProportionalHazardsPValue = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} ProportionalHazardsPValueGlobal
    ##
    ## Proportional hazards test over the whole model
    ##
    ## A scalar @math{p} value testing the proportional hazards assumption for
    ## all predictors at once.  This property is read-only.
    ##
    ## @end deftp
    ProportionalHazardsPValueGlobal = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} LikelihoodRatioTestPValue
    ##
    ## Likelihood ratio test against the null model
    ##
    ## A scalar @math{p} value comparing the fitted model with the model that
    ## carries no predictors.  This property is read-only.
    ##
    ## @end deftp
    LikelihoodRatioTestPValue = [];

    ## -*- texinfo -*-
    ## @deftp {CoxModel} {property} VariableInfo
    ##
    ## Information about the variables
    ##
    ## A table with one row per variable, its row names the variable names,
    ## and the variables @qcode{Class}, @qcode{Range}, @qcode{InModel} and
    ## @qcode{IsCategorical}.  This property is read-only.
    ##
    ## @end deftp
    VariableInfo = [];

  endproperties

  properties (Access = private, Hidden)
    b_          = [];   # coefficient vector aligned with the encoded columns
    X_          = [];   # encoded predictor matrix used in the fit
    T_          = [];   # response as passed to coxphfit
    cens_       = [];   # censoring indicator
    freq_       = [];   # frequency weights
    strata_     = [];   # stratum label per observation ([] when unstratified)
    baserow_    = [];   # baseline actually used, one row per stratum
    ties_       = 'breslow';
    catinfo_    = [];   # categorical encoding, for new data
    catcols_    = [];   # logical mask of categorical predictors
    catbase_    = [];   # logical mask of encoded columns that are indicators
    formulastr_ = '';
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
      fprintf ("\n  Cox proportional hazards regression model:\n");
      if (! isempty (this.formulastr_))
        fprintf ("      %s\n", this.formulastr_);
      endif
      if (! isempty (this.Coefficients))
        fprintf ("\n  Coefficients:\n\n");
        disp (this.Coefficients);
      endif
      fprintf ("\n");
      if (! isempty (this.LogLikelihood))
        fprintf ("Log-likelihood: %g\n", this.LogLikelihood);
      endif
      if (! isempty (this.LikelihoodRatioTestPValue))
        fprintf ("Likelihood ratio test vs. constant model: p-value = %g\n", ...
                 this.LikelihoodRatioTestPValue);
      endif
    endfunction

    ## Class specific subscripted reference.
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("CoxModel: () indexing is not supported.", ...
                         "  Use dot notation for properties."));
        case '{}'
          error (strcat ("CoxModel: {} indexing is not supported.", ...
                         "  Use dot notation for properties."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("CoxModel.subsref: property name must be a", ...
                           " character vector."));
          endif
          if (ismethod (this, s.subs))
            [varargout{1:nargout}] = builtin ('subsref', this, [s, chain_s]);
            return;
          endif
          try
            out = this.(s.subs);
          catch
            error ("CoxModel.subsref: unknown property '%s'.", s.subs);
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
    ## @deftypefn  {CoxModel} {@var{mdl} =} CoxModel (@var{X}, @var{T})
    ## @deftypefnx {CoxModel} {@var{mdl} =} CoxModel (@var{tbl}, @var{respvar})
    ## @deftypefnx {CoxModel} {@var{mdl} =} CoxModel (@dots{}, @var{Name}, @var{Value})
    ##
    ## Fit a Cox proportional hazards regression model.
    ##
    ## @code{@var{mdl} = CoxModel (@var{X}, @var{T})} fits the model to the
    ## @math{n}-by-@math{p} numeric predictor matrix @var{X} and the
    ## @math{n}-by-1 vector of event times @var{T}.  @var{T} may instead be an
    ## @math{n}-by-2 matrix giving a @math{(start, stop]} interval of exposure,
    ## the counting process form.
    ##
    ## @code{@var{mdl} = CoxModel (@var{tbl}, @var{respvar})} takes the data
    ## from the table @var{tbl}, using the variable named @var{respvar} as the
    ## response and every other variable as a predictor.  A @code{categorical}
    ## variable is encoded as indicator columns.
    ##
    ## @var{X} must not contain a constant column: the model has no constant
    ## term, since any constant is absorbed into the baseline hazard.
    ##
    ## The following @var{Name}/@var{Value} pairs are accepted:
    ##
    ## @multitable @columnfractions 0.25 0.75
    ## @headitem Name @tab Value
    ## @item @qcode{"Baseline"} @tab The @var{X} values at which the baseline
    ## hazard is computed, either a scalar or a 1-by-@math{p} vector.  The
    ## default is the mean of each numeric predictor and zero for each indicator
    ## column of a categorical predictor, taken within each stratum.
    ## @item @qcode{"Beta"} @tab The starting value of the iteration, a vector
    ## of length @math{p}.  The default is @code{0.01 ./ std (@var{X})}.
    ## @item @qcode{"CategoricalPredictors"} @tab The predictors to treat as
    ## categorical, given as column indices, a logical vector, or a cell array
    ## of predictor names.  Table variables of class @code{categorical} are
    ## detected without this argument.
    ## @item @qcode{"Censoring"} @tab A logical or 0/1 vector of length
    ## @math{n}, where 1 marks an observation right-censored at its recorded
    ## time.  The default is a vector of zeros.
    ## @item @qcode{"Frequency"} @tab A vector of length @math{n} of
    ## non-negative values giving the number of observations each row
    ## represents, or a weight.  The default is a vector of ones.
    ## @item @qcode{"OptimizationOptions"} @tab A structure of iteration
    ## settings, as built by @code{statset ("fitcox")}.  The fields used are
    ## @qcode{"MaxIter"}, @qcode{"TolX"} and @qcode{"Display"}.
    ## @item @qcode{"PredictorNames"} @tab A cell array of @math{p} predictor
    ## names.  The default is @qcode{"X1"}, @qcode{"X2"}, and so on, or the
    ## table variable names.
    ## @item @qcode{"Stratification"} @tab A vector of length @math{n} of
    ## stratum labels.  Each stratum carries its own baseline hazard and its own
    ## risk sets, while the coefficients are shared across all of them.
    ## @item @qcode{"TieBreakMethod"} @tab The method of handling tied event
    ## times, either @qcode{"breslow"} (default) or @qcode{"efron"}.
    ## @end multitable
    ##
    ## @end deftypefn

    function this = CoxModel (X, T, varargin)

      if (nargin < 2)
        error ("CoxModel: too few input arguments.");
      endif

      ## --- the data ----------------------------------------------------
      istbl = isa (X, 'table');
      if (istbl)
        if (! ((ischar (T) && isrow (T))
               || (isa (T, 'string') && isscalar (T))))
          error (strcat ("CoxModel: RESPVAR must be a character vector", ...
                         " naming a variable of the data table."));
        endif
        resp_name = char (T);
        all_names = X.Properties.VariableNames;
        r_idx     = find (strcmp (all_names, resp_name), 1);
        if (isempty (r_idx))
          error ("CoxModel: '%s' is not a variable of the data table.", ...
                 resp_name);
        endif
        T_val      = X.(resp_name);
        pred_names = all_names;
        pred_names(r_idx) = [];
        tbl = X;
      else
        if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
          error ("CoxModel: X must be a real numeric matrix.");
        endif
        T_val      = T;
        resp_name  = 'y';
        pred_names = arrayfun (@(k) sprintf ("X%d", k), 1:columns (X), ...
                               'UniformOutput', false);
        tbl = [];
      endif
      p_raw = numel (pred_names);
      if (p_raw < 1)
        error ("CoxModel: the model needs at least one predictor.");
      endif

      ## --- optional arguments ------------------------------------------
      Baseline  = [];
      Beta0     = [];
      CatPreds  = [];
      Censoring = [];
      Frequency = [];
      Options   = [];
      PredNames = {};
      Strata    = [];
      Ties      = 'breslow';

      if (mod (numel (varargin), 2) != 0)
        error ("CoxModel: optional arguments must be name-value pairs.");
      endif
      for i = 1:2:numel (varargin)
        name = varargin{i};
        if (! (ischar (name) && isrow (name)))
          error ("CoxModel: parameter name must be a character vector.");
        endif
        switch (lower (name))
          case 'baseline'
            Baseline = varargin{i+1};
          case 'beta'
            Beta0 = varargin{i+1};
          case 'categoricalpredictors'
            CatPreds = varargin{i+1};
          case 'censoring'
            Censoring = varargin{i+1};
          case 'frequency'
            Frequency = varargin{i+1};
          case 'optimizationoptions'
            Options = varargin{i+1};
          case 'predictornames'
            PredNames = varargin{i+1};
          case 'stratification'
            Strata = varargin{i+1};
          case 'tiebreakmethod'
            Ties = varargin{i+1};
          otherwise
            error ("CoxModel: unknown parameter name '%s'.", name);
        endswitch
      endfor

      if (! isempty (PredNames))
        if (! (iscellstr (PredNames) && numel (PredNames) == p_raw))
          error (strcat ("CoxModel: 'PredictorNames' must be a cell array", ...
                         " of one character vector per predictor."));
        endif
        pred_names = PredNames(:)';
      endif
      if (! (ischar (Ties) && any (strcmpi (Ties, {'breslow', 'efron'}))))
        error (strcat ("CoxModel: 'TieBreakMethod' must be either", ...
                       " 'breslow' or 'efron'."));
      endif
      Ties = lower (Ties);

      ## --- the predictor data, under the names the model will report ----
      ## A renamed table variable is renamed once, here, so that every later
      ## lookup goes through PredictorNames alone.
      if (istbl)
        X_raw   = [];
        tbl_p   = tbl(:, setdiff (1:numel (all_names), r_idx));
        tbl_p.Properties.VariableNames = pred_names;
        n_total = rows (tbl);
        data_in = tbl_p;
      else
        X_raw   = X;
        tbl_p   = [];
        n_total = rows (X);
        data_in = X;
      endif

      ## --- which predictors are categorical ----------------------------
      cat_cols = false (1, p_raw);
      if (istbl)
        for j = 1:p_raw
          col = tbl_p.(pred_names{j});
          cat_cols(j) = isa (col, 'categorical') || iscell (col);
        endfor
      endif
      if (! isempty (CatPreds))
        cat_cols = cat_cols | categorical_mask (CatPreds, pred_names);
      endif

      ## --- numeric codes, then indicator columns -----------------------
      [X_num, cat_levels] = raw_to_codes (data_in, X_raw, tbl_p, pred_names, ...
                                          cat_cols, n_total);

      ## The baseline hazard plays the part of an intercept -- it is what the
      ## reference level of a categorical predictor is folded into -- so the
      ## encoding is reference coded even though the model has no constant
      ## term, which is what MATLAB does here too.
      [X_enc, enc_names, cat_info] = encode_categorical (X_num, cat_cols, ...
                                                         pred_names, ...
                                                         cat_levels, true);
      p = columns (X_enc);

      ## Mark the encoded columns that came from a categorical predictor: they
      ## enter the default baseline as zero rather than as their mean.
      ## encode_categorical gives a categorical predictor one column per level
      ## bar the reference and a numeric one a single column, so the encoded
      ## columns are walked in the same order to find which came from which.
      cat_enc = false (1, p);
      k = 0;
      for j = 1:p_raw
        if (cat_cols(j))
          nlev = numel (cat_levels{j}) - 1;
          cat_enc(k+1:k+nlev) = true;
          k += nlev;
        else
          k += 1;
        endif
      endfor

      ## --- the fit ------------------------------------------------------
      args = {};
      if (! isempty (Censoring))
        args = [args, {'Censoring', Censoring}];
      endif
      if (! isempty (Frequency))
        args = [args, {'Frequency', Frequency}];
      endif
      if (! isempty (Beta0))
        args = [args, {'B0', Beta0}];
      endif
      if (! isempty (Options))
        args = [args, {'Options', Options}];
      endif
      if (! isempty (Strata))
        args = [args, {'Strata', Strata}];
      endif
      args = [args, {'Ties', Ties}];

      ## The baseline each stratum is centred on.  A numeric predictor enters
      ## it as its frequency-weighted mean and a categorical one as zero, the
      ## reference level being what the baseline hazard already carries.
      w = ones (n_total, 1);
      if (! isempty (Frequency))
        w = Frequency(:);
      endif
      base_given = ! isempty (Baseline);
      if (isempty (Strata))
        slev = [];
        nS   = 1;
        smask = {true(n_total, 1)};
      else
        slev  = unique (Strata(:));
        nS    = numel (slev);
        smask = arrayfun (@(v) Strata(:) == v, slev, 'UniformOutput', false);
      endif
      brow = zeros (nS, p);   # baseline reported and predicted against
      bdef = zeros (nS, p);   # baseline coxphfit centres on by default
      for s = 1:nS
        m         = smask{s};
        bdef(s,:) = sum (w(m) .* X_enc(m,:), 1) / sum (w(m));
        if (base_given)
          bl = Baseline;
          if (isscalar (bl))
            bl = repmat (bl, 1, p);
          endif
          brow(s,:) = bl(:)';
        else
          brow(s,:) = bdef(s,:);
          brow(s,cat_enc) = 0;
        endif
      endfor

      ## coxphfit takes one baseline for every stratum, so a per-stratum
      ## baseline is reached by fitting on its own default and rescaling the
      ## hazard afterwards -- exact, the hazard being scaled by exp (B b) and
      ## nothing else.  Only a stratified model with a categorical predictor
      ## needs it; anything else states its baseline outright.
      rescale = false;
      if (base_given)
        args = [args, {'Baseline', Baseline}];
      elseif (nS == 1)
        args = [args, {'Baseline', brow(1,:)}];
      elseif (any (cat_enc))
        rescale = true;
      endif

      [b, logl, H, stats] = coxphfit (X_enc, T_val, args{:});

      if (rescale)
        for s = 1:nS
          m      = H(:,3) == slev(s);
          H(m,2) = H(m,2) * exp ((brow(s,:) - bdef(s,:)) * b);
        endfor
      endif

      ## --- properties ---------------------------------------------------
      this.b_       = b;
      this.X_       = X_enc;
      this.T_       = T_val;
      this.cens_    = Censoring;
      this.freq_    = Frequency;
      this.strata_  = Strata;
      this.ties_    = Ties;
      this.catinfo_ = cat_info;
      this.catcols_ = cat_cols;
      this.catbase_ = cat_enc;

      this.NumPredictors         = p_raw;
      this.LogLikelihood         = logl;
      this.Hazard                = H;
      this.PredictorNames        = pred_names;
      this.ResponseName          = resp_name;
      this.CoefficientCovariance = stats.covb;
      this.StandardError         = stats.se(:);
      this.LikelihoodRatioTestPValue = stats.LikelihoodRatioTestP;

      this.Coefficients = table (b(:), stats.se(:), stats.z(:), stats.p(:), ...
        'VariableNames', {'Beta', 'SE', 'zStat', 'pValue'}, ...
        'RowNames', enc_names(:));

      this.Residuals = table (stats.csres, stats.devres, stats.martres, ...
        stats.schres, stats.sschres, stats.scores, stats.sscores, ...
        'VariableNames', {'CoxSnell', 'Deviance', 'Martingale', ...
                          'Schoenfeld', 'ScaledSchoenfeld', 'Score', ...
                          'ScaledScore'});

      ## The baseline actually used, one row per stratum, kept for prediction.
      ## Reported as given when it was given -- a scalar stays a scalar -- and
      ## otherwise as the rows the fit was centred on.
      this.Stratification = slev;
      this.baserow_       = brow;
      if (base_given)
        this.Baseline = Baseline;
      else
        this.Baseline = brow;
      endif

      ## --- formula and per-variable information -------------------------
      var_names_all = [pred_names, {resp_name}];
      terms         = [eye(p_raw), zeros(p_raw, 1)];
      this.Formula  = LinearFormula (terms, var_names_all, ...
                                     'ResponseName', resp_name);
      this.formulastr_ = char (this.Formula);

      vi_class   = cell (p_raw + 1, 1);
      vi_range   = cell (p_raw + 1, 1);
      vi_inmodel = [true(p_raw, 1); false];
      vi_iscat   = [cat_cols(:); false];
      for j = 1:p_raw
        if (istbl)
          col = tbl_p.(pred_names{j});
        else
          col = X_raw(:, j);
        endif
        [vi_class{j}, vi_range{j}] = variable_class_and_range (col);
      endfor
      [vi_class{end}, vi_range{end}] = variable_class_and_range (T_val(:,end));
      this.VariableInfo = table (vi_class, vi_range, vi_inmodel, vi_iscat, ...
        'VariableNames', {'Class', 'Range', 'InModel', 'IsCategorical'}, ...
        'RowNames', var_names_all(:));

      ## --- the proportional hazards assumption --------------------------
      [phz, phg] = ph_assumption_test (stats, T_val);
      this.ProportionalHazardsPValue       = phz;
      this.ProportionalHazardsPValueGlobal = phg;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CoxModel} {@var{ci} =} coefci (@var{obj})
    ## @deftypefnx {CoxModel} {@var{ci} =} coefci (@var{obj}, @var{level})
    ##
    ## Confidence intervals for the coefficients of a Cox model.
    ##
    ## @code{@var{ci} = coefci (@var{obj})} returns a two-column matrix with one
    ## row per coefficient, holding the 95% confidence interval of each.
    ##
    ## @code{@var{ci} = coefci (@var{obj}, @var{level})} uses a
    ## @code{100 (1 - @var{level})}% interval.  @var{level} must be a positive
    ## scalar smaller than 1; it is a significance level, not a coverage.
    ##
    ## @end deftypefn

    function ci = coefci (this, level)

      if (nargin < 2)
        level = 0.05;
      endif
      if (! (isnumeric (level) && isreal (level) && isscalar (level)
             && level > 0 && level < 1))
        error (strcat ("CoxModel.coefci: LEVEL must be a real scalar", ...
                       " greater than 0 and smaller than 1."));
      endif
      z  = norminv (1 - level / 2);
      se = this.StandardError(:);
      b  = this.b_(:);
      ci = [b - z * se, b + z * se];

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CoxModel} {@var{obj} =} discardResiduals (@var{obj})
    ##
    ## Drop the stored residuals of a Cox model.
    ##
    ## @code{@var{obj} = discardResiduals (@var{obj})} returns the model with an
    ## empty @code{Residuals} property.  The residual table holds one row per
    ## observation and is the largest thing a fitted model carries, so
    ## discarding it makes a model that is only going to be used for prediction
    ## considerably smaller.  Nothing else about the model changes.
    ##
    ## @end deftypefn

    function this = discardResiduals (this)
      this.Residuals = table ();
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CoxModel} {@var{hr} =} hazardratio (@var{obj}, @var{X})
    ## @deftypefnx {CoxModel} {@var{hr} =} hazardratio (@var{obj}, @var{X}, @var{S})
    ## @deftypefnx {CoxModel} {@var{hr} =} hazardratio (@dots{}, @qcode{"Baseline"}, @var{B})
    ##
    ## Hazard of a Cox model relative to its baseline.
    ##
    ## @code{@var{hr} = hazardratio (@var{obj}, @var{X})} returns the hazard at
    ## the predictor values @var{X} relative to the baseline the model was
    ## fitted with, @code{exp ((@var{X} - @var{B}) * @var{b})}.  @var{X} has one
    ## row per evaluation point and is a numeric matrix, or a table when the
    ## model was fitted from one.
    ##
    ## @code{@var{hr} = hazardratio (@var{obj}, @var{X}, @var{S})} gives the
    ## stratum of each row of @var{X}, and is required when the model is
    ## stratified, each stratum having its own baseline.
    ##
    ## @code{@var{hr} = hazardratio (@dots{}, "Baseline", @var{B})} evaluates
    ## the ratio against the baseline @var{B} instead, either a scalar or a row
    ## vector with one element per encoded predictor column.
    ##
    ## @end deftypefn

    function hr = hazardratio (this, varargin)

      if (numel (varargin) < 1)
        error ("CoxModel.hazardratio: X is required.");
      endif
      [Xq, Sq, opts] = split_prediction_args (this, 'hazardratio', varargin);
      base = prediction_baseline (this, Sq, opts.Baseline, 'hazardratio');
      hr   = exp (sum ((Xq - base) .* this.b_(:)', 2));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CoxModel} {@var{tbl} =} linhyptest (@var{obj})
    ##
    ## Sequential tests on the coefficients of a Cox model.
    ##
    ## @code{@var{tbl} = linhyptest (@var{obj})} returns a table with one row
    ## per predictor, whose @math{k}-th row tests the hypothesis that the
    ## coefficients of the @math{k}-th and every later predictor are jointly
    ## zero.  The @code{Predictor} column names the predictors the hypothesis
    ## leaves in the model, so its first row reads @qcode{"Empty Model"} and
    ## tests every coefficient at once, and its last row tests the last
    ## coefficient alone, reproducing that coefficient's own p-value.
    ##
    ## Each test is a Wald test on the fitted model, not a refit.
    ##
    ## @end deftypefn

    function tbl = linhyptest (this)

      p     = numel (this.b_);
      pvals = zeros (p, 1);
      names = cell (p, 1);
      for k = 1:p
        idx      = k:p;
        bk       = this.b_(idx);
        Vk       = this.CoefficientCovariance(idx, idx);
        stat     = bk(:)' * (Vk \ bk(:));
        pvals(k) = 1 - chi2cdf (stat, numel (idx));
        if (k == 1)
          names{k} = 'Empty Model';
        else
          names{k} = strjoin (this.Coefficients.Properties.RowNames(1:k-1)', ...
                              ', ');
        endif
      endfor
      tbl = table (names, pvals, ...
                   'VariableNames', {'Predictor', 'pValue'});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CoxModel} {} plotSurvival (@var{obj})
    ## @deftypefnx {CoxModel} {} plotSurvival (@var{obj}, @var{X})
    ## @deftypefnx {CoxModel} {} plotSurvival (@var{obj}, @var{X}, @var{S})
    ## @deftypefnx {CoxModel} {@var{h} =} plotSurvival (@dots{})
    ##
    ## Plot the survival function of a Cox model.
    ##
    ## @code{plotSurvival (@var{obj})} draws the survival function at the
    ## model's baseline as a stairstep plot.  @code{plotSurvival (@var{obj},
    ## @var{X})} draws it at the predictor values @var{X}, one curve per row,
    ## and @var{S} gives the stratum of each row when the model is stratified.
    ## A stratified model with no @var{X} draws one curve per stratum.
    ##
    ## @code{@var{h} = plotSurvival (@dots{})} returns the handles of the
    ## stairstep lines.
    ##
    ## @end deftypefn

    function varargout = plotSurvival (this, varargin)

      [s, tout] = survival (this, varargin{:});
      ## Unstratified curves come back as the columns of a matrix over one
      ## grid; stratified ones as a cell array, each on its own grid.
      if (iscell (s))
        curves = s;
        times  = tout;
      else
        curves = num2cell (s, 1);
        times  = repmat ({tout}, 1, columns (s));
      endif
      h = zeros (numel (curves), 1);
      hold_state = ishold ();
      for k = 1:numel (curves)
        h(k) = stairs (times{k}, curves{k});
        if (k == 1)
          hold on;
        endif
      endfor
      if (! hold_state)
        hold off;
      endif
      xlabel (this.ResponseName);
      ylabel ("Survival probability");
      title ("Cox proportional hazards model");
      if (nargout > 0)
        varargout{1} = h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CoxModel} {@var{s} =} survival (@var{obj})
    ## @deftypefnx {CoxModel} {@var{s} =} survival (@var{obj}, @var{X})
    ## @deftypefnx {CoxModel} {@var{s} =} survival (@var{obj}, @var{X}, @var{S})
    ## @deftypefnx {CoxModel} {@var{s} =} survival (@dots{}, @var{Name}, @var{Value})
    ## @deftypefnx {CoxModel} {[@var{s}, @var{T}] =} survival (@dots{})
    ##
    ## Survival function of a Cox model.
    ##
    ## @code{@var{s} = survival (@var{obj})} returns the survival probability at
    ## the model's baseline, evaluated at each row of the @code{Hazard}
    ## property.  @code{@var{s} = survival (@var{obj}, @var{X})} evaluates it at
    ## the predictor values @var{X}, and @var{S} gives the stratum of each row
    ## when the model is stratified.  For a stratified model @var{s} is a cell
    ## array holding one column vector per curve.
    ##
    ## @code{[@var{s}, @var{T}] = survival (@dots{})} also returns the times the
    ## probabilities refer to.
    ##
    ## The following @var{Name}/@var{Value} pairs are accepted:
    ##
    ## @multitable @columnfractions 0.3 0.7
    ## @headitem Name @tab Value
    ## @item @qcode{"Time"} @tab The times at which to evaluate the survival
    ## function.  The default is the model's own event times.  The baseline
    ## survival is interpolated linearly between them and raised to the hazard
    ## ratio of @var{X}.
    ## @item @qcode{"ExtrapolationMethod"} @tab How to evaluate a time outside
    ## the model's event times: @qcode{"nearest"} (default), @qcode{"linear"},
    ## @qcode{"next"}, @qcode{"previous"}, or @qcode{"none"}.  @qcode{"none"}
    ## returns @qcode{NaN} outside the range, as do @qcode{"next"} above it and
    ## @qcode{"previous"} below it.
    ## @end multitable
    ##
    ## @end deftypefn

    function [s, Tout] = survival (this, varargin)

      [Xq, Sq, opts] = split_prediction_args (this, 'survival', varargin);
      base = prediction_baseline (this, Sq, [], 'survival');
      if (isempty (Xq))
        hr = ones (max (1, numel (Sq)), 1);
      else
        hr = exp (sum ((Xq - base) .* this.b_(:)', 2));
      endif

      strat = this.Stratification;
      if (isempty (strat))
        grid_t = {this.Hazard(:,1)};
        grid_H = {this.Hazard(:,2)};
        which_grid = ones (numel (hr), 1);
      else
        nS     = numel (strat);
        grid_t = cell (nS, 1);
        grid_H = cell (nS, 1);
        for k = 1:nS
          m         = this.Hazard(:,3) == strat(k);
          grid_t{k} = this.Hazard(m,1);
          grid_H{k} = this.Hazard(m,2);
        endfor
        if (isempty (Sq))
          which_grid = (1:nS)';
          hr         = ones (nS, 1);
        else
          [~, which_grid] = ismember (Sq(:), strat);
        endif
      endif

      n_out = numel (which_grid);
      s     = cell (n_out, 1);
      Tout  = cell (n_out, 1);
      for k = 1:n_out
        g   = which_grid(k);
        tg  = grid_t{g};
        S0  = exp (-grid_H{g});
        if (isempty (opts.Time))
          s{k}    = S0 .^ hr(k);
          Tout{k} = tg;
        else
          tq      = opts.Time(:);
          s{k}    = interp_survival (tg, S0, tq, opts.Extrap) .^ hr(k);
          Tout{k} = tq;
        endif
      endfor

      ## Unstratified curves all share one time grid, so they come back as the
      ## columns of a matrix, one per row of X.  A stratified model gives each
      ## curve its own stratum's grid, so those stay a cell array.
      if (isempty (strat))
        s    = cell2mat (s(:)');
        Tout = Tout{1};
      else
        s    = s(:)';
        Tout = Tout(:)';
      endif

    endfunction

  endmethods

endclassdef

## Resolve a 'CategoricalPredictors' value to a logical mask over the
## predictors.
function mask = categorical_mask (spec, pred_names)

  p    = numel (pred_names);
  mask = false (1, p);
  if (islogical (spec))
    if (numel (spec) != p)
      error (strcat ("CoxModel: a logical 'CategoricalPredictors' must", ...
                     " have one element per predictor."));
    endif
    mask = spec(:)';
  elseif (isnumeric (spec))
    if (any (spec != fix (spec)) || any (spec < 1) || any (spec > p))
      error (strcat ("CoxModel: 'CategoricalPredictors' indices must be", ...
                     " integers between 1 and the number of predictors."));
    endif
    mask(spec) = true;
  elseif (iscellstr (spec) || ischar (spec))
    if (ischar (spec))
      spec = {spec};
    endif
    for k = 1:numel (spec)
      j = find (strcmp (pred_names, spec{k}), 1);
      if (isempty (j))
        error ("CoxModel: '%s' is not a predictor name.", spec{k});
      endif
      mask(j) = true;
    endfor
  else
    error (strcat ("CoxModel: 'CategoricalPredictors' must be indices, a", ...
                   " logical vector, or predictor names."));
  endif

endfunction

## Split the arguments of a prediction method into the predictor values, the
## stratum labels, and the name-value options.
function [Xq, Sq, opts] = split_prediction_args (this, caller, args)

  opts = struct ('Time', [], 'Extrap', 'nearest', 'Baseline', []);
  ## Each method takes only its own options: 'Baseline' belongs to
  ## hazardratio, the times to survival.  Accepting an option and ignoring it
  ## would claim a capability that is not there.
  if (strcmp (caller, 'hazardratio'))
    keys = {'baseline'};
  else
    keys = {'time', 'extrapolationmethod'};
  endif
  ## Every name is known to the scan that separates the positional arguments
  ## from the pairs, so that a name meant for the other method is reported as
  ## the wrong name rather than as a malformed argument list.
  scan_keys = {'time', 'extrapolationmethod', 'baseline'};

  ## Positional arguments run until the first name that starts a pair.
  npos = 0;
  while (npos < numel (args) && npos < 2)
    a = args{npos+1};
    if ((ischar (a) && isrow (a)) && any (strcmpi (a, scan_keys))
        && numel (args) > npos + 1)
      break;
    endif
    npos++;
  endwhile

  Xq = [];
  Sq = [];
  if (npos >= 1)
    Xq = args{1};
  endif
  if (npos >= 2)
    Sq = args{2};
  endif

  rest = args(npos+1:end);
  if (mod (numel (rest), 2) != 0)
    error ("CoxModel.%s: optional arguments must be name-value pairs.", caller);
  endif
  for i = 1:2:numel (rest)
    name = rest{i};
    if (! (ischar (name) && isrow (name)))
      error ("CoxModel.%s: parameter name must be a character vector.", caller);
    endif
    if (! any (strcmpi (name, keys)))
      error ("CoxModel.%s: unknown parameter name '%s'.", caller, name);
    endif
    switch (lower (name))
      case 'time'
        t = rest{i+1};
        if (! (isnumeric (t) && isreal (t) && isvector (t)))
          error ("CoxModel.%s: 'Time' must be a real vector.", caller);
        endif
        opts.Time = t;
      case 'extrapolationmethod'
        m = rest{i+1};
        valid = {'nearest', 'linear', 'next', 'previous', 'none'};
        if (! (ischar (m) && any (strcmpi (m, valid))))
          error (strcat ("CoxModel.%s: 'ExtrapolationMethod' must be one", ...
                         " of 'nearest', 'linear', 'next', 'previous', or", ...
                         " 'none'."), caller);
        endif
        opts.Extrap = lower (m);
      case 'baseline'
        opts.Baseline = rest{i+1};
      otherwise
        error ("CoxModel.%s: unknown parameter name '%s'.", caller, name);
    endswitch
  endfor

  ## Encode a table or raw predictor matrix the same way the fit did.
  if (! isempty (Xq))
    Xq = encode_new_data (this, Xq, caller);
  endif

  if (! isempty (this.Stratification) && ! isempty (Xq) && isempty (Sq))
    error (strcat ("CoxModel.%s: the model is stratified, so the stratum", ...
                   " of each row of X is required."), caller);
  endif
  if (! isempty (Sq) && ! isempty (Xq) && numel (Sq) != rows (Xq))
    error (strcat ("CoxModel.%s: S must have one element for each row", ...
                   " of X."), caller);
  endif

endfunction

## Encode new predictor values with the fitted model's own encoding.
function Xq = encode_new_data (this, Xq, caller)

  if (isa (Xq, 'table'))
    names = Xq.Properties.VariableNames;
    X_num = zeros (rows (Xq), numel (this.PredictorNames));
    for j = 1:numel (this.PredictorNames)
      k = find (strcmp (names, this.PredictorNames{j}), 1);
      if (isempty (k))
        error ("CoxModel.%s: X has no variable named '%s'.", caller, ...
               this.PredictorNames{j});
      endif
      col = Xq.(this.PredictorNames{j});
      if (isa (col, 'categorical') || iscell (col))
        levels = this.catinfo_.levels{j};
        [tf, ic] = ismember (cellstr (col), levels);
        if (! all (tf))
          error (strcat ("CoxModel.%s: X holds a level of '%s' that the", ...
                         " model was not fitted with."), caller, ...
                 this.PredictorNames{j});
        endif
        X_num(:,j) = ic;
      else
        X_num(:,j) = double (col(:));
      endif
    endfor
    Xq = encode_categorical (X_num, this.catcols_, this.PredictorNames, ...
                             this.catinfo_.levels, true);
  else
    if (! (isnumeric (Xq) && isreal (Xq)))
      error ("CoxModel.%s: X must be a real numeric matrix.", caller);
    endif
    if (columns (Xq) == numel (this.PredictorNames) && any (this.catcols_))
      Xq = encode_categorical (Xq, this.catcols_, this.PredictorNames, ...
                               this.catinfo_.levels, true);
    endif
  endif
  if (columns (Xq) != numel (this.b_))
    error (strcat ("CoxModel.%s: X must have one column for each", ...
                   " predictor."), caller);
  endif

endfunction

## The baseline each row of a prediction is measured against.
function base = prediction_baseline (this, Sq, given, caller)

  p = numel (this.b_);
  if (! isempty (given))
    if (! (isnumeric (given) && isreal (given)
           && (isscalar (given) || numel (given) == p)))
      error (strcat ("CoxModel.%s: 'Baseline' must be a scalar or a row", ...
                     " vector with one element per predictor column."), caller);
    endif
    if (isscalar (given))
      base = repmat (given, 1, p);
    else
      base = given(:)';
    endif
    return;
  endif

  if (isempty (this.Stratification) || isempty (Sq))
    base = this.baserow_(1,:);
  else
    [tf, idx] = ismember (Sq(:), this.Stratification);
    if (! all (tf))
      error (strcat ("CoxModel.%s: S holds a stratum the model was not", ...
                     " fitted with."), caller);
    endif
    base = this.baserow_(idx,:);
  endif

endfunction

## Interpolate a baseline survival curve, honouring the extrapolation rule.
function sq = interp_survival (tg, S0, tq, method)

  ## A repeated time carries two values -- the curve steps there -- and the
  ## lower one is what the interval to its right starts from.
  [tu, iu] = unique (tg(:), 'last');
  su       = S0(iu);

  below = tq < tu(1);
  above = tq > tu(end);
  if (numel (tu) < 2)
    ## One event time leaves nothing to interpolate between.
    sq = NaN (size (tq));
    sq(tq == tu(1)) = su(1);
  else
    sq = interp1 (tu, su, tq, 'linear');
  endif

  switch (method)
    case 'nearest'
      sq(below) = S0(1);
      sq(above) = su(end);
    case 'linear'
      sq(below) = S0(1);
      if (any (above) && numel (tu) > 1)
        slope = (su(end) - su(end-1)) / (tu(end) - tu(end-1));
        sq(above) = min (1, max (0, su(end) + slope * (tq(above) - tu(end))));
      endif
    case 'next'
      sq(below) = S0(1);
      sq(above) = NaN;
    case 'previous'
      sq(below) = NaN;
      sq(above) = su(end);
    case 'none'
      sq(below) = NaN;
      sq(above) = NaN;
  endswitch

endfunction

## Grambsch-Therneau test of the proportional hazards assumption, on event
## times replaced by their ranks.  The scaled Schoenfeld residual of a
## covariate is regressed on time; a slope that is not zero is a hazard ratio
## that moves with time, which is what proportionality denies.
##
## Tied event times take their mid-rank, the average of the ranks they span.
## Only ties tell the rankings apart -- dense ranks, ordinal ranks and
## mid-ranks agree on distinct times -- and the mid-rank is the one that
## reproduces MATLAB, measured under both tie-breaking methods.
function [pz, pg] = ph_assumption_test (stats, T)

  ev = ! any (isnan (stats.schres), 2);
  if (! any (ev))
    pz = NaN (1, columns (stats.schres));
    pg = NaN;
    return;
  endif
  s = stats.schres(ev,:);
  t = T(ev,end);
  d = sum (ev);
  V = stats.covb;

  [ut, ~, ic] = unique (t);
  [~, ord]    = sort (t);
  r           = zeros (numel (t), 1);
  r(ord)      = 1:numel (t);
  g           = zeros (numel (t), 1);
  for j = 1:numel (ut)
    m    = ic == j;
    g(m) = mean (r(m));
  endfor
  w    = g - mean (g);
  varx = sum (w .^ 2);
  if (varx == 0)
    pz = NaN (1, columns (s));
    pg = NaN;
    return;
  endif

  ## The scaled Schoenfeld residual is beta + d V s; the constant drops out of
  ## the centred regression, so d V s is what is correlated with time.
  u  = (w' * (s * V * d))';
  z  = (u .^ 2) ./ (varx * d * diag (V));
  pz = (1 - chi2cdf (z, 1))';
  pg = 1 - chi2cdf (u' * ((V * d * varx) \ u), columns (s));

endfunction

%!demo
%! ## Fit a Cox proportional hazards model and read its coefficients
%! X = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! T = [4; 6; 8; 11; 13; 16; 18; 21; 25; 30];
%! mdl = fitcox (X, T)

%!demo
%! ## The survival function at the model's baseline
%! X = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! T = [4; 6; 8; 11; 13; 16; 18; 21; 25; 30];
%! mdl = fitcox (X, T);
%! [s, t] = survival (mdl);
%! [t, s]

%!shared X, T, C, Tt, S, F, T2
%! X  = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! T  = [4; 6; 8; 11; 13; 16; 18; 21; 25; 30];
%! C  = [0; 0; 1; 0; 0; 1; 0; 0; 1; 0];
%! Tt = [4; 4; 6; 6; 8; 8; 11; 11; 13; 13];
%! S  = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];
%! F  = [1; 2; 1; 1; 3; 1; 1; 2; 1; 1];
%! T2 = [0 4; 0 6; 2 8; 0 11; 3 13; 0 16; 5 18; 0 21; 7 25; 0 30];

## The fitted coefficients are those of coxphfit
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (mdl.Coefficients.Beta, ...
%!               [-1.3886093196382836; 4.3814437183613322], 1e-8);
%! assert_equal (mdl.Coefficients.SE, ...
%!               [0.52737766743917369; 1.8537400210498443], 1e-8);
%! assert_equal (mdl.Coefficients.zStat, ...
%!               [-2.6330453589001133; 2.3635696853973904], 1e-8);
%! assert_equal (mdl.Coefficients.pValue, ...
%!               [0.0084623045168834322; 0.018099822319697333], 1e-10);

%!test
%! mdl = CoxModel (X, T);
%! assert_equal (mdl.LogLikelihood, -8.8069639381632356, 1e-10);
%! assert_equal (mdl.LikelihoodRatioTestPValue, 0.0018409958426933715, 1e-10);
%! assert_equal (mdl.NumPredictors, 2);

## Coefficient table shape and names
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (mdl.Coefficients.Properties.VariableNames, ...
%!               {'Beta', 'SE', 'zStat', 'pValue'});
%! assert_equal (mdl.Coefficients.Properties.RowNames, {'X1'; 'X2'});
%! assert_equal (size (mdl.Coefficients), [2, 4]);

## The covariance, standard errors, and baseline
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (mdl.CoefficientCovariance, ...
%!               [0.27812720411358371, -0.88932589928247008; ...
%!                -0.88932589928247008, 3.4363520656418767], 1e-8);
%! assert_equal (mdl.StandardError, ...
%!               [0.52737766743917369; 1.8537400210498443], 1e-8);
%! assert_equal (mdl.Baseline, [5.9, 0.5], 1e-12);

## The baseline cumulative hazard, as coxphfit returns it
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (size (mdl.Hazard), [11, 2]);
%! assert_equal (mdl.Hazard(1,:), [4, 0]);
%! assert_equal (mdl.Hazard(end,2), 39.969396549779539, 1e-6);

## Names and formula are consistent with one another
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (mdl.PredictorNames, {'X1', 'X2'});
%! assert_equal (mdl.ResponseName, 'y');
%! assert_equal (char (mdl.Formula), 'y ~ X1 + X2');

%!test
%! mdl = CoxModel (X, T, 'PredictorNames', {'age', 'trt'});
%! assert_equal (mdl.PredictorNames, {'age', 'trt'});
%! assert_equal (mdl.ResponseName, 'y');
%! assert_equal (char (mdl.Formula), 'y ~ age + trt');
%! assert_equal (mdl.Coefficients.Properties.RowNames, {'age'; 'trt'});

## Per-variable information
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (mdl.VariableInfo.Properties.VariableNames, ...
%!               {'Class', 'Range', 'InModel', 'IsCategorical'});
%! assert_equal (mdl.VariableInfo.Properties.RowNames, {'X1'; 'X2'; 'y'});
%! assert_equal (mdl.VariableInfo.InModel, [true; true; false]);
%! assert_equal (mdl.VariableInfo.IsCategorical, [false; false; false]);

## The residual table carries every type coxphfit computes
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (mdl.Residuals.Properties.VariableNames, ...
%!               {'CoxSnell', 'Deviance', 'Martingale', 'Schoenfeld', ...
%!                'ScaledSchoenfeld', 'Score', 'ScaledScore'});
%! assert_equal (size (mdl.Residuals), [10, 7]);
%! assert_equal (mdl.Residuals.Martingale(1), 0.6260391348376092, 1e-8);

%!test
%! mdl = CoxModel (X, T);
%! mdl = discardResiduals (mdl);
%! assert_equal (isempty (mdl.Residuals), true);
%! assert_equal (mdl.LogLikelihood, -8.8069639381632356, 1e-10);

## The proportional hazards assumption test
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (mdl.ProportionalHazardsPValue, ...
%!               [0.43865495271858179, 0.71413497767000234], 1e-8);
%! assert_equal (mdl.ProportionalHazardsPValueGlobal, ...
%!               0.53179255135825187, 1e-8);

%!test
%! mdl = CoxModel (X, T, 'Censoring', C);
%! assert_equal (mdl.ProportionalHazardsPValue, ...
%!               [0.38377818764906124, 0.62418297490071406], 1e-8);
%! assert_equal (mdl.ProportionalHazardsPValueGlobal, ...
%!               0.56302705208766479, 1e-8);

## Tied event times take their mid-rank; only a tie tells the rankings apart
%!test
%! mdl = CoxModel (X, Tt, 'Censoring', C);
%! assert_equal (mdl.ProportionalHazardsPValue, ...
%!               [0.21978449883707185, 0.47185588628217012], 1e-8);
%! assert_equal (mdl.ProportionalHazardsPValueGlobal, ...
%!               0.33816199611424813, 1e-8);

%!test
%! mdl = CoxModel (X, Tt, 'Censoring', C, 'TieBreakMethod', 'efron');
%! assert_equal (mdl.Coefficients.SE, ...
%!               [0.40807722610877606; 1.7349017938498392], 1e-8);
%! assert_equal (mdl.ProportionalHazardsPValue, ...
%!               [0.1946475346918406, 0.44447162515215921], 1e-8);
%! assert_equal (mdl.ProportionalHazardsPValueGlobal, ...
%!               0.29473648771539396, 1e-8);

## Confidence intervals
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (coefci (mdl), ...
%!               [-2.422250554069806, -0.35496808520676093; ...
%!                0.74818004040311514, 8.0147073963195492], 1e-8);
%! assert_equal (coefci (mdl, 0.01), ...
%!               [-2.747044169465374, -0.030174469811192983; ...
%!                -0.39347414902021249, 9.1563615857428768], 1e-8);

## Sequential Wald tests
%!test
%! mdl = CoxModel ([X, [1;3;2;5;4;6;8;7;9;10]], T);
%! tbl = linhyptest (mdl);
%! assert_equal (size (tbl), [3, 2]);
%! assert_equal (tbl.Predictor, {'Empty Model'; 'X1'; 'X1, X2'});
%! assert_equal (tbl.pValue, [0.11426498364159436; 0.066700738069398954; ...
%!                            0.040978941270458896], 1e-8);

## The last sequential test is the last coefficient's own p-value
%!test
%! mdl = CoxModel ([X, [1;3;2;5;4;6;8;7;9;10]], T);
%! tbl = linhyptest (mdl);
%! assert_equal (tbl.pValue(end), mdl.Coefficients.pValue(end), 1e-10);

## Hazard ratios
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (hazardratio (mdl, X(1,:)), 25.149914260347884, 1e-6);
%! assert_equal (hazardratio (mdl, X(1,:), 'Baseline', 0), ...
%!               0.062211299031685895, 1e-8);

%!test
%! mdl = CoxModel (X, T);
%! hr = hazardratio (mdl, X);
%! assert_equal (numel (hr), 10);
%! assert_equal (hr(end), 0.030119684486042686, 1e-8);

## The survival function
%!test
%! mdl = CoxModel (X, T);
%! [s, t] = survival (mdl);
%! assert_equal (numel (s), 11);
%! assert_equal (s(1), 1);
%! assert_equal (s(2), 0.98524073172292947, 1e-8);
%! assert_equal (t, mdl.Hazard(:,1));

%!test
%! mdl = CoxModel (X, T);
%! s = survival (mdl, X(1,:));
%! assert_equal (s(2), 0.6880038362205394, 1e-8);
%! assert_equal (s(3), 0.37858862123606984, 1e-8);

## Interpolation is linear on the survival scale, then raised to the ratio
%!test
%! mdl = CoxModel (X, T);
%! s = survival (mdl, 'Time', [4; 5; 6; 7]);
%! assert_equal (s, [0.98524073172292947; 0.97367819309820003; ...
%!                   0.96211565447347058; 0.91995051472240608], 1e-8);

%!test
%! mdl = CoxModel (X, T);
%! s = survival (mdl, X(1,:), 'Time', [5; 12; 20]);
%! assert_equal (s, [0.51126892454645156; 9.5014819976515837e-06; ...
%!                   1.7030500478184053e-37], 1e-8);

## Outside the event times, the extrapolation rule decides
%!test
%! mdl = CoxModel (X, T);
%! tq = [1; 2; 3; 31; 40];
%! assert_equal (survival (mdl, 'Time', tq), ...
%!               [1; 1; 1; 4.3803784471438163e-18; ...
%!                4.3803784471438163e-18], 1e-8);
%! assert_equal (survival (mdl, 'Time', tq, 'ExtrapolationMethod', 'none'), ...
%!               [NaN; NaN; NaN; NaN; NaN]);

%!test
%! mdl = CoxModel (X, T);
%! tq = [1; 2; 3; 31; 40];
%! s = survival (mdl, 'Time', tq, 'ExtrapolationMethod', 'previous');
%! assert_equal (isnan (s(1:3)), [true; true; true]);
%! assert_equal (s(4), 4.3803784471438163e-18, 1e-8);
%! s = survival (mdl, 'Time', tq, 'ExtrapolationMethod', 'next');
%! assert_equal (s(1:3), [1; 1; 1]);
%! assert_equal (isnan (s(4:5)), [true; true]);

## Stratification
%!test
%! mdl = CoxModel (X, T, 'Censoring', C, 'Stratification', S);
%! assert_equal (mdl.Stratification, [1; 2]);
%! assert_equal (size (mdl.Hazard), [9, 3]);
%! assert_equal (mdl.Coefficients.Beta, ...
%!               [-0.7705046389; 3.1181961544], 1e-6);
%! assert_equal (mdl.ProportionalHazardsPValue, ...
%!               [0.59926200661233664, 0.7801291832899343], 1e-6);

%!test
%! mdl = CoxModel (X, T, 'Censoring', C, 'Stratification', S);
%! s = survival (mdl);
%! assert_equal (iscell (s), true);
%! assert_equal (numel (s), 2);
%! assert_equal (numel (s{1}), 5);
%! assert_equal (numel (s{2}), 4);

%!test
%! mdl = CoxModel (X, T, 'Censoring', C, 'Stratification', S);
%! assert_equal (hazardratio (mdl, X(1,:), 1), 1.825643762907144, 1e-6);

## Unstratified curves are the columns of a matrix, one per row of X
%!test
%! mdl = CoxModel (X, T);
%! [s, t] = survival (mdl, X(1:3,:));
%! assert_equal (size (s), [11, 3]);
%! assert_equal (size (t), [11, 1]);
%! assert_equal (s(2,1), 0.6880038362205394, 1e-8);
%! assert_equal (s(2,2), 0.62879787425614297, 1e-8);
%! assert_equal (s(3,3), 0.78484832904143043, 1e-8);

## A stratified curve follows the stratum of its own row of X
%!test
%! mdl = CoxModel (X, T, 'Censoring', C, 'Stratification', S);
%! s = survival (mdl, X(1:2,:), [2; 1]);
%! assert_equal (iscell (s), true);
%! assert_equal (numel (s{1}), 4);
%! assert_equal (numel (s{2}), 5);

## A baseline given to the method overrides the fitted one
%!test
%! mdl = CoxModel (X, T);
%! assert_equal (hazardratio (mdl, X(1,:), 'Baseline', [1 1]), ...
%!               0.0031195920528550632, 1e-8);

## Discarding the residuals leaves prediction intact
%!test
%! mdl = discardResiduals (CoxModel (X, T));
%! s = survival (mdl, X(1,:));
%! assert_equal (s(2), 0.6880038362205394, 1e-8);

## Frequency weights may be given as a row
%!test
%! mdl = CoxModel (X, T, 'Frequency', F');
%! assert_equal (mdl.Coefficients.Beta, ...
%!               [-1.4447610177189139; 4.6300103378728581], 1e-6);

## plotSurvival draws one stairstep line per curve
%!test
%! f = figure ('visible', 'off');
%! unwind_protect
%!   mdl = CoxModel (X, T, 'Censoring', C, 'Stratification', S);
%!   h = plotSurvival (mdl);
%!   assert_equal (numel (h), 2);
%!   assert_equal (all (ishghandle (h)), true);
%! unwind_protect_cleanup
%!   close (f);
%! end_unwind_protect

## One line per curve, a curve per row of X when the model is not stratified
%!test
%! f = figure ('visible', 'off');
%! unwind_protect
%!   mdl = CoxModel (X, T);
%!   assert_equal (numel (plotSurvival (mdl)), 1);
%!   assert_equal (numel (plotSurvival (mdl, X(1:3,:))), 3);
%! unwind_protect_cleanup
%!   close (f);
%! end_unwind_protect

## The counting process form
%!test
%! mdl = CoxModel (X, T2, 'Censoring', C);
%! assert_equal (mdl.Coefficients.Beta, [-1.0104; 3.2523], 1e-4);

## Frequency weights
%!test
%! mdl = CoxModel (X, T, 'Frequency', F);
%! assert_equal (mdl.Coefficients.Beta, ...
%!               [-1.4447610177189139; 4.6300103378728581], 1e-6);

## A starting value changes nothing but the path taken to the answer
%!test
%! mdl = CoxModel (X, T, 'Beta', [0.1; -0.1]);
%! assert_equal (mdl.Coefficients.Beta, ...
%!               [-1.3886093196382836; 4.3814437183613322], 1e-6);

## An options structure is accepted
%!test
%! mdl = CoxModel (X, T, 'OptimizationOptions', statset ('fitcox'));
%! assert_equal (mdl.LogLikelihood, -8.8069639381632356, 1e-10);

## Errors
%!error<CoxModel: too few input arguments.> CoxModel (1)
%!error<CoxModel: X must be a real numeric matrix.> CoxModel ({1}, [1; 2])
%!error<CoxModel: optional arguments must be name-value pairs.> ...
%! CoxModel (X, T, 'Censoring')
%!error<CoxModel: unknown parameter name 'Ties'.> ...
%! CoxModel (X, T, 'Ties', 'efron')
%!error<CoxModel: 'TieBreakMethod' must be either 'breslow' or 'efron'.> ...
%! CoxModel (X, T, 'TieBreakMethod', 'exact')
%!error<CoxModel: 'PredictorNames' must be a cell array of one character vector per predictor.> ...
%! CoxModel (X, T, 'PredictorNames', {'only_one'})
%!error<CoxModel.coefci: LEVEL must be a real scalar greater than 0 and smaller than 1.> ...
%! coefci (CoxModel (X, T), 1.5)
%!error<CoxModel.hazardratio: X is required.> hazardratio (CoxModel (X, T))
%!error<CoxModel.survival: unknown parameter name 'Baseline'.> ...
%! survival (CoxModel (X, T), X(1,:), 'Baseline', 0)
%!error<CoxModel.hazardratio: unknown parameter name 'Time'.> ...
%! hazardratio (CoxModel (X, T), X(1,:), 'Time', 5)
%!error<CoxModel.survival: 'ExtrapolationMethod' must be one of 'nearest', 'linear', 'next', 'previous', or 'none'.> ...
%! survival (CoxModel (X, T), 'Time', 5, 'ExtrapolationMethod', 'cubic')
%!error<CoxModel.survival: the model is stratified, so the stratum of each row of X is required.> ...
%! survival (CoxModel (X, T, 'Stratification', S), X(1,:))
%!error<CoxModel: \(\) indexing is not supported.  Use dot notation for properties.> ...
%! CoxModel (X, T)(1)
%!error<CoxModel.subsref: unknown property 'nosuch'.> ...
%! CoxModel (X, T).nosuch
