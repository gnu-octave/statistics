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

    ## Number of observations used in the fit (missing rows excluded).
    NumObservations = [];

    ## Deviance of the fitted model.
    Deviance = [];

    ## Error (residual) degrees of freedom.
    DFE = [];

    ## Dispersion parameter (estimated or fixed at 1).
    Dispersion = [];

    ## True if the dispersion parameter was estimated from the data.
    DispersionEstimated = [];

    ## Structure describing the response distribution (Name and its deviance and
    ## variance functions).
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

    ## Character-vector representation of the model formula.
    Formula = '';

  endproperties

  properties (Access = private, Hidden)
    b_        = [];   # coefficient vector (intercept first when present)
    stats_    = [];   # glmfit stats struct (for prediction CIs)
    distr_    = '';   # distribution name
    linkarg_  = [];   # link specification passed to glmfit/glmval
    linkname_ = '';   # resolved link name (or 'power(p)')
    intercept_ = true;
    binomsize_ = [];  # BinomialSize (trials) for the binomial family
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
                         " supported.  Use dot notation to access properties."));
        case '{}'
          error (strcat ("GeneralizedLinearModel: {} indexing is not", ...
                         " supported.  Use dot notation to access properties."));
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
    ## @deftypefn {GeneralizedLinearModel} {@var{mdl} =} GeneralizedLinearModel (@var{X}, @var{y})
    ## @deftypefnx {GeneralizedLinearModel} {@var{mdl} =} GeneralizedLinearModel (@var{X}, @var{y}, @var{Name}, @var{Value})
    ##
    ## Fit a generalized linear model.  Prefer the @code{fitglm} function; see
    ## its help for the accepted @var{Name}/@var{Value} pairs.
    ##
    ## @end deftypefn
    function this = GeneralizedLinearModel (X, y, varargin)

      if (nargin == 0)
        return;   # empty object
      endif
      if (nargin < 2)
        error ("GeneralizedLinearModel: X and Y are required.");
      endif
      if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
        error ("GeneralizedLinearModel: X must be a real matrix.");
      endif
      if (! (isnumeric (y) && isreal (y) && isvector (y)))
        error ("GeneralizedLinearModel: Y must be a real vector.");
      endif
      y = y(:);
      if (rows (X) != numel (y))
        error (strcat ("GeneralizedLinearModel: X and Y must have the same", ...
                       " number of observations."));
      endif

      ## Defaults and Name-Value parsing.
      distr     = 'normal';
      linkarg   = [];
      weights   = [];
      offset    = [];
      binomsize = [];
      intercept = true;
      estdisp   = [];
      varnames  = {};
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("GeneralizedLinearModel: optional arguments must be", ...
                       " Name-Value pairs."));
      endif
      for k = 1:2:numel (varargin)
        switch (lower (varargin{k}))
          case 'distribution';  distr     = lower (varargin{k+1});
          case 'link';          linkarg   = varargin{k+1};
          case 'weights';       weights   = varargin{k+1};
          case 'offset';        offset    = varargin{k+1};
          case 'binomialsize';  binomsize = varargin{k+1};
          case 'intercept';     intercept = logical (varargin{k+1});
          case 'dispersionflag'; estdisp  = logical (varargin{k+1});
          case 'varnames';      varnames  = varargin{k+1};
          otherwise
            error (strcat ("GeneralizedLinearModel: unknown parameter", ...
                           " name '%s'."), varargin{k});
        endswitch
      endfor
      if (! any (strcmp (distr, {'normal', 'binomial', 'poisson', 'gamma', ...
                                 'inverse gaussian'})))
        error ("GeneralizedLinearModel: unknown distribution '%s'.", distr);
      endif

      p = columns (X);

      ## Assemble the glmfit call.
      args = {};
      if (! isempty (linkarg))
        args = [args, {'link', linkarg}];
      endif
      if (! isempty (weights))
        args = [args, {'weights', weights}];
      endif
      if (! isempty (offset))
        args = [args, {'offset', offset(:)}];
      endif
      args = [args, {'constant', ternary(intercept, 'on', 'off')}];
      if (! isempty (estdisp))
        args = [args, {'estdisp', ternary(estdisp, 'on', 'off')}];
      endif

      ## For the binomial family with a trials count, form the [successes trials]
      ## response that glmfit expects.
      yfit = y;
      if (strcmp (distr, 'binomial') && ! isempty (binomsize))
        N = binomsize(:);
        if (isscalar (N))
          N = N * ones (numel (y), 1);
        endif
        yfit = [y .* N, N];
      endif

      [b, dev, stats] = glmfit (X, yfit, distr, args{:});

      ## Coefficient names, formula, and variable names.
      if (isempty (varnames))
        prednames = arrayfun (@(j) sprintf ('x%d', j), 1:p, ...
                              'UniformOutput', false);
        respname = 'y';
      else
        prednames = varnames(1:p)(:)';
        respname  = varnames{p + 1};
      endif
      if (intercept)
        coefnames = [{'(Intercept)'}, prednames];
      else
        coefnames = prednames;
      endif
      terms = strjoin (prednames, ' + ');
      if (intercept)
        this.Formula = sprintf ('%s ~ 1 + %s', respname, terms);
      else
        this.Formula = sprintf ('%s ~ %s', respname, terms);
      endif

      ## Link and distribution descriptors.  LINKSPEC is what glmfit/glmval and
      ## getlinkfunctions understand (a name or a numeric power exponent);
      ## LINKNAME is the human-readable form for display.
      if (isempty (linkarg))
        linkspec = default_link_spec (distr);
      else
        linkspec = linkarg;
      endif
      linkname = link_name (linkspec);
      [flink, dlink, ilink] = getlinkfunctions (linkspec);
      this.Link = struct ('Name', linkname, 'Link', flink, ...
                          'Derivative', dlink, 'Inverse', ilink);
      this.Distribution = struct ('Name', distr, ...
                                  'DevianceFunction', [], ...
                                  'VarianceFunction', []);

      ## Coefficient and residual tables.
      this.Coefficients = table (b(:), stats.se(:), stats.t(:), stats.p(:), ...
        'VariableNames', {'Estimate', 'SE', 'tStat', 'pValue'}, ...
        'RowNames', coefnames(:));
      this.Residuals = table (stats.resid(:), stats.residp(:), ...
        stats.residd(:), stats.resida(:), ...
        'VariableNames', {'Raw', 'Pearson', 'Deviance', 'Anscombe'});

      ## Fitted values on both scales.
      valargs = {'constant', ternary(intercept, 'on', 'off')};
      if (! isempty (offset))
        valargs = [valargs, {'offset', offset(:)}];
      endif
      if (strcmp (distr, 'binomial') && ! isempty (binomsize))
        valargs = [valargs, {'size', N}];
      endif
      yhat = glmval (b, X, linkspec, valargs{:});
      eta = ilink_eta (X, b, intercept, offset);
      this.Fitted = struct ('Response', yhat, 'LinearPredictor', eta);

      ## Scalar summaries.
      this.b_         = b;
      this.stats_     = stats;
      this.distr_     = distr;
      this.linkarg_   = linkspec;
      this.linkname_  = linkname;
      this.intercept_ = intercept;
      this.binomsize_ = binomsize;

      this.CoefficientNames         = coefnames;
      this.CoefficientCovariance    = stats.covb;
      this.NumCoefficients          = numel (b);
      this.NumEstimatedCoefficients = numel (b);
      this.NumPredictors            = p;
      this.NumObservations          = numel (y) - sum (any (isnan ([X, y]), 2));
      this.Deviance                 = dev;
      this.DFE                      = stats.dfe;
      this.Dispersion               = stats.s;
      this.DispersionEstimated      = logical (stats.estdisp);
      this.Offset                   = offset;
      this.ResponseName             = respname;
      this.PredictorNames           = prednames;
      this.VariableNames            = [prednames, {respname}];

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearModel} {@var{yhat} =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {GeneralizedLinearModel} {[@var{yhat}, @var{yci}] =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {GeneralizedLinearModel} {[@dots{}] =} predict (@dots{}, @var{Name}, @var{Value})
    ##
    ## Predict the response of the model @var{mdl} at the new predictor data
    ## @var{Xnew} (an @math{m}-by-@math{p} numeric matrix).  Predictions are on
    ## the mean (response) scale.  With two outputs, @var{yci} is an
    ## @math{m}-by-2 matrix of confidence intervals.  The @qcode{'Alpha'} pair
    ## sets the confidence level to @math{100 (1 - @var{Alpha})%} (default 0.05).
    ##
    ## @end deftypefn
    function [yhat, yci] = predict (mdl, Xnew, varargin)
      if (nargin < 2)
        error ("GeneralizedLinearModel.predict: Xnew is required.");
      endif
      alpha = 0.05;
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
      valargs = {'constant', ternary(mdl.intercept_, 'on', 'off'), ...
                 'confidence', 1 - alpha};
      if (! isempty (offnew))
        valargs = [valargs, {'offset', offnew(:)}];
      endif
      if (! isempty (mdl.binomsize_))
        N = mdl.binomsize_(:);
        if (isscalar (N))
          N = N * ones (rows (Xnew), 1);
        endif
        valargs = [valargs, {'size', N}];
      endif
      if (nargout > 1)
        [yhat, ylo, yhi] = glmval (mdl.b_, Xnew, mdl.linkarg_, mdl.stats_, ...
                                   valargs{:});
        yci = [yhat - ylo, yhat + yhi];
      else
        yhat = glmval (mdl.b_, Xnew, mdl.linkarg_, valargs{:});
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

## Response distribution's canonical link specification (name, or numeric power
## exponent) as understood by glmfit/glmval/getlinkfunctions.
function spec = default_link_spec (distr)
  switch (distr)
    case 'normal';            spec = 'identity';
    case 'binomial';          spec = 'logit';
    case 'poisson';           spec = 'log';
    case 'gamma';             spec = 'reciprocal';
    case 'inverse gaussian';  spec = -2;
  endswitch
endfunction

## Human-readable name of a user-supplied link specification.
function name = link_name (linkarg)
  if (ischar (linkarg))
    name = linkarg;
  elseif (isnumeric (linkarg) && isscalar (linkarg))
    name = sprintf ('power(%g)', linkarg);
  else
    name = 'custom';
  endif
endfunction

## Linear predictor eta = [1 X] * b (+ offset).
function eta = ilink_eta (X, b, intercept, offset)
  if (intercept)
    eta = b(1) + X * b(2:end);
  else
    eta = X * b;
  endif
  if (! isempty (offset))
    eta = eta + offset(:);
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
