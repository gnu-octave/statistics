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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

classdef RepeatedMeasuresModel
  ## -*- texinfo -*-
  ## @deftp {statistics} RepeatedMeasuresModel
  ##
  ## Repeated measures model.
  ##
  ## A @qcode{RepeatedMeasuresModel} object holds a multivariate linear model
  ## of several responses measured on the same subjects, the repeated
  ## measures, on between-subject predictors.  Each row of the data is a
  ## subject and each response a measurement taken on it; the between-subject
  ## model describes how the mean of every response depends on the
  ## predictors, and the within-subject design describes how the responses
  ## relate to one another, as the levels of one or more within-subject
  ## factors.
  ##
  ## The between-subject coefficients are estimated by least squares, one
  ## column per response, with categorical predictors in effects coding: a
  ## predictor of @math{L} levels contributes @math{L - 1} columns named
  ## @qcode{@var{name}_@var{level}}, the last level coded -1 in all of them.
  ## @code{ranova} tests the within-subject effects, reporting the p-values
  ## corrected for departures from sphericity beside the uncorrected ones,
  ## @code{epsilon} gives those corrections and @code{mauchly} tests
  ## sphericity itself.
  ##
  ## Create a @qcode{RepeatedMeasuresModel} object with @code{fitrm}.
  ##
  ## @seealso{fitrm, anova2, manova1}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} BetweenDesign
    ##
    ## Design of the between-subject predictors
    ##
    ## The table the model was fitted to, one row per subject, responses
    ## included.  A subject missing a response or a predictor is kept here,
    ## although the fit leaves it out.  This property is read-only.
    ##
    ## @end deftp
    BetweenDesign       = [];

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} BetweenModel
    ##
    ## Between-subject model
    ##
    ## A character vector holding the right side of the model formula, with
    ## its intercept, as in @qcode{'1 + species'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    BetweenModel        = '';

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} BetweenFactorNames
    ##
    ## Names of the between-subject predictors
    ##
    ## A cell row of character vectors naming the variables the
    ## between-subject model draws on, empty for a model with an intercept
    ## alone.  This property is read-only.
    ##
    ## @end deftp
    BetweenFactorNames  = {};

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} ResponseNames
    ##
    ## Names of the responses
    ##
    ## A cell row of character vectors naming the repeated measures, in the
    ## order of the columns of @code{Coefficients}.  This property is
    ## read-only.
    ##
    ## @end deftp
    ResponseNames       = {};

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} WithinDesign
    ##
    ## Design of the within-subject factors
    ##
    ## A table with one row per response, named by the responses, and one
    ## variable per within-subject factor.  Without a design given at fitting
    ## it holds a single factor @qcode{Time} of the values 1 to @math{k} for
    ## @math{k} responses.  This property is read-only.
    ##
    ## @end deftp
    WithinDesign        = [];

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} WithinModel
    ##
    ## Within-subject model
    ##
    ## @qcode{'separatemeans'}, @qcode{'orthogonalcontrasts'}, a formula over
    ## the within-subject factors or a contrast matrix, as given at fitting.
    ## This property is read-only.
    ##
    ## @end deftp
    WithinModel         = 'separatemeans';

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} WithinFactorNames
    ##
    ## Names of the within-subject factors
    ##
    ## A cell row of character vectors, the variable names of
    ## @code{WithinDesign}.  This property is read-only.
    ##
    ## @end deftp
    WithinFactorNames   = {};

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} Coefficients
    ##
    ## Estimated between-subject coefficients
    ##
    ## A table with one row per coefficient, named after the terms of the
    ## between-subject model, and one variable per response.  This property is
    ## read-only.
    ##
    ## @end deftp
    Coefficients        = [];

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} Covariance
    ##
    ## Estimated covariance of the responses
    ##
    ## A table holding the covariance of the residuals of the responses,
    ## their cross products over @code{DFE}, with rows and variables named by
    ## the responses.  This property is read-only.
    ##
    ## @end deftp
    Covariance          = [];

    ## -*- texinfo -*-
    ## @deftp {RepeatedMeasuresModel} {property} DFE
    ##
    ## Error degrees of freedom
    ##
    ## The number of subjects the fit used less the number of between-subject
    ## coefficients.  This property is read-only.
    ##
    ## @end deftp
    DFE                 = [];

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    ## The between-subject design matrix, one row per subject of
    ## BetweenDesign, NaN on a row the fit left out.
    DesignMatrix        = [];
  endproperties

  properties (Access = private, Hidden)
    ## Complete subjects' design matrix and residuals, the coefficients as a
    ## matrix, for each between term its name, its columns and its variables
    ## (the intercept's first, and none), and the coding of each variable.
    X_          = [];
    R_          = [];
    B_          = [];
    TermNames_  = {};
    TermCols_   = {};
    TermVars_   = {};
    Coding_     = struct ();
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

    ## Custom display of the model, grouped as MATLAB groups it.
    function disp (this)
      fprintf ("\n  RepeatedMeasuresModel with properties:\n\n");
      fprintf ("   Between Subjects:\n");
      fprintf ("%+22s: [%dx%d table]\n", 'BetweenDesign', ...
               size (this.BetweenDesign));
      fprintf ("%+22s: %s\n", 'ResponseNames', cellrow (this.ResponseNames));
      fprintf ("%+22s: %s\n", 'BetweenFactorNames', ...
               cellrow (this.BetweenFactorNames));
      fprintf ("%+22s: '%s'\n\n", 'BetweenModel', this.BetweenModel);
      fprintf ("   Within Subjects:\n");
      fprintf ("%+22s: [%dx%d table]\n", 'WithinDesign', ...
               size (this.WithinDesign));
      fprintf ("%+22s: %s\n", 'WithinFactorNames', ...
               cellrow (this.WithinFactorNames));
      if (ischar (this.WithinModel))
        fprintf ("%+22s: '%s'\n\n", 'WithinModel', this.WithinModel);
      else
        fprintf ("%+22s: [%dx%d double]\n\n", 'WithinModel', ...
                 size (this.WithinModel));
      endif
      fprintf ("   Estimates:\n");
      fprintf ("%+22s: [%dx%d table]\n", 'Coefficients', ...
               size (this.Coefficients));
      fprintf ("%+22s: [%dx%d table]\n\n", 'Covariance', ...
               size (this.Covariance));
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{rm} =} RepeatedMeasuresModel (@var{t}, @var{modelspec})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{rm} =} RepeatedMeasuresModel (@var{t}, @var{modelspec}, @var{name}, @var{value})
    ##
    ## Fit a repeated measures model.
    ##
    ## @code{@var{rm} = RepeatedMeasuresModel (@var{t}, @var{modelspec})} fits
    ## the repeated measures in the table @var{t} on the between-subject
    ## model given by the formula @var{modelspec}, a character vector or a
    ## string scalar of the form @qcode{'@var{responses} ~ @var{terms}'}.
    ## The responses are a range of the table's variables, as in
    ## @qcode{'y1-y6'}, which takes every variable from @qcode{y1} to
    ## @qcode{y6} in the table's order, a comma list, as in
    ## @qcode{'y1,y2,y3'}, or both.  The terms are a Wilkinson formula over the
    ## other variables, as in @qcode{'species'} or @qcode{'g*x'}, or
    ## @qcode{'1'} for an intercept alone.  A @qcode{categorical}, logical,
    ## text or string variable is a categorical predictor and any other a
    ## continuous one.
    ##
    ## A subject missing a response or a predictor is left out of the fit.
    ##
    ## The following name-value arguments are accepted.
    ##
    ## @multitable @columnfractions 0.25 0.75
    ## @headitem Name @tab Value
    ## @item @qcode{'WithinDesign'} @tab The design of the within-subject
    ## factors: a table with one row per response and one variable per factor,
    ## or a numeric vector with one element per response, which becomes a
    ## single factor @qcode{Time}.  The default is @qcode{Time} holding the
    ## values 1 to @math{k} for @math{k} responses.
    ##
    ## @item @qcode{'WithinModel'} @tab The within-subject model:
    ## @qcode{'separatemeans'}, the default, which compares the means of the
    ## responses; @qcode{'orthogonalcontrasts'}, which tests the orthogonal
    ## polynomial trends over a single numeric within-subject factor; a
    ## formula over the within-subject factors, as in @qcode{'A*B'}; or a
    ## contrast matrix with one row per response.
    ## @end multitable
    ##
    ## @qcode{'orthogonalcontrasts'} is refused at fitting unless the
    ## within-subject design holds a single numeric factor.  MATLAB accepts
    ## the fit and refuses the model only when a test is asked of it.
    ##
    ## @seealso{fitrm}
    ## @end deftypefn
    function this = RepeatedMeasuresModel (t, modelspec, varargin)

      if (nargin < 2)
        error ("RepeatedMeasuresModel: too few input arguments.");
      endif
      if (! istable (t))
        error ("RepeatedMeasuresModel: T must be a table.");
      endif
      if (isa (modelspec, 'string') && isscalar (modelspec))
        modelspec = char (modelspec);
      endif
      if (! (ischar (modelspec) && isrow (modelspec)))
        error (strcat ("RepeatedMeasuresModel: MODELSPEC must be a", ...
                       " character vector or a string scalar."));
      endif

      ## Parse optional paired arguments.  An empty WithinDesign is the
      ## default Time factor.
      optNames = {'WithinDesign', 'WithinModel'};
      dfValues = {[], 'separatemeans'};
      [WD, WM, args] = parsePairedArguments (optNames, dfValues, varargin(:));
      if (! isempty (args))
        error ("RepeatedMeasuresModel: invalid optional paired argument.");
      endif

      ## Responses and between-subject terms from the formula
      vnames = t.Properties.VariableNames;
      parts = strsplit (modelspec, '~');
      if (numel (parts) != 2 || isempty (strtrim (parts{1})) ...
                             || isempty (strtrim (parts{2})))
        error (strcat ("RepeatedMeasuresModel: MODELSPEC must be of the", ...
                       " form 'responses ~ terms'."));
      endif
      [resp, errmsg] = responseNames (parts{1}, vnames);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel: %s", errmsg);
      endif
      rhs = strtrim (parts{2});
      [terms, intercept, errmsg] = betweenTerms (rhs, vnames, resp);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel: %s", errmsg);
      endif

      Y = zeros (rows (t), numel (resp));
      for j = 1:numel (resp)
        y = t.(resp{j});
        if (! (isnumeric (y) && isreal (y) && isvector (y)))
          error ("RepeatedMeasuresModel: response '%s' must be numeric.", ...
                 resp{j});
        endif
        Y(:,j) = double (y(:));
      endfor
      k = numel (resp);

      ## Effects-coded design over the complete subjects
      vars = unique ([terms{:}], 'stable');
      [X, cnames, tcols, tnames, bad, coding] = effectsDesign (t, terms, ...
                                                           intercept, vars);
      ok = ! bad & all (! isnan (Y), 2);
      ## A continuous predictor enters a marginal mean at its mean over the
      ## subjects the fit used, as in MATLAB
      for i = 1:numel (vars)
        if (! coding.(vars{i}).iscat)
          coding.(vars{i}).mean = mean (double (t.(vars{i})(ok)));
        endif
      endfor
      Xc = X(ok,:);
      Yc = Y(ok,:);
      if (rank (Xc) < columns (Xc) || rows (Xc) <= rank (Xc))
        error (strcat ("RepeatedMeasuresModel: the between-subject", ...
                       " design is rank deficient or leaves no error", ...
                       " degrees of freedom."));
      endif
      B = Xc \ Yc;
      R = Yc - Xc * B;
      dfe = rows (Xc) - rank (Xc);

      ## Within-subject design
      if (isempty (WD))
        WD = table ((1:k)', 'VariableNames', {'Time'}, 'RowNames', resp);
      elseif (isnumeric (WD) && isvector (WD))
        if (numel (WD) != k)
          error (strcat ("RepeatedMeasuresModel: 'WithinDesign' has %d", ...
                         " points, where %d are required."), numel (WD), k);
        endif
        WD = table (double (WD(:)), 'VariableNames', {'Time'}, ...
                    'RowNames', resp);
      elseif (istable (WD))
        if (rows (WD) != k)
          error (strcat ("RepeatedMeasuresModel: 'WithinDesign' has %d", ...
                         " points, where %d are required."), rows (WD), k);
        endif
        WD.Properties.RowNames = resp;
      else
        error (strcat ("RepeatedMeasuresModel: 'WithinDesign' must be", ...
                       " a table or a numeric vector."));
      endif
      wnames = WD.Properties.VariableNames;
      [~, errmsg] = withinTerms (WM, WD, wnames, k);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel: %s", errmsg);
      endif
      if (isa (WM, 'string'))
        WM = char (WM);
      endif

      ## Populate the properties
      this.BetweenDesign = t;
      this.ResponseNames = resp;
      this.BetweenFactorNames = cell (1, 0);
      if (! isempty (vars))
        this.BetweenFactorNames = vars;
      endif
      this.BetweenModel = modelText (rhs, intercept);
      this.WithinDesign = WD;
      this.WithinFactorNames = wnames;
      this.WithinModel = WM;
      this.Coefficients = array2table (B, 'VariableNames', resp, ...
                                       'RowNames', cnames);
      this.Covariance = array2table (R' * R / dfe, 'VariableNames', resp, ...
                                     'RowNames', resp);
      this.DFE = dfe;
      X(! ok,:) = NaN;
      this.DesignMatrix = X;
      this.X_ = Xc;
      this.R_ = R;
      this.B_ = B;
      this.TermNames_ = tnames;
      this.TermCols_ = tcols;
      this.TermVars_ = terms;
      this.Coding_ = coding;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{tbl} =} ranova (@var{rm})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{tbl} =} ranova (@var{rm}, @qcode{'WithinModel'}, @var{WM})
    ## @deftypefnx {RepeatedMeasuresModel} {[@var{tbl}, @var{A}, @var{C}, @var{D}] =} ranova (@dots{})
    ##
    ## Repeated measures analysis of variance.
    ##
    ## @code{@var{tbl} = ranova (@var{rm})} tests the within-subject effects of
    ## the repeated measures model @var{rm}: whether the means of the
    ## responses differ, and whether each between-subject term changes that
    ## difference.  @var{tbl} holds one row per test and one per error term,
    ## with the variables @qcode{SumSq}, @qcode{DF}, @qcode{MeanSq}, @qcode{F}
    ## and @qcode{pValue}, and beside them the p-values under the
    ## Greenhouse-Geisser, Huynh-Feldt and lower bound corrections for
    ## departures from sphericity, @qcode{pValueGG}, @qcode{pValueHF} and
    ## @qcode{pValueLB}.  The corrections are reported and never applied.
    ##
    ## @code{ranova (@var{rm}, @qcode{'WithinModel'}, @var{WM})} tests the
    ## within-subject model @var{WM} instead, which takes the forms the
    ## @qcode{'WithinModel'} argument of @code{fitrm} takes.  The default is
    ## @qcode{'separatemeans'} whatever the model was fitted with, as in
    ## MATLAB.  Under @qcode{'separatemeans'} or a contrast matrix the rows
    ## are named after the within-subject factor, or @qcode{Time} when there
    ## are several.  A formula or @qcode{'orthogonalcontrasts'} gives one
    ## block of rows per within-subject term, the constant term first, whose
    ## tests are the between-subject tests of the mean of the responses.
    ##
    ## @code{[@var{tbl}, @var{A}, @var{C}, @var{D}] = ranova (@dots{})} also
    ## returns the hypothesis of each test as @math{A B C = D}: @var{A} a cell
    ## column holding the between-subject hypothesis matrix of each term,
    ## @var{C} the within-subject contrast, a cell row of one per term when
    ## there are several, and @var{D} zero.
    ##
    ## @code{'orthogonalcontrasts'} works here where MATLAB R2024a and R2026a
    ## both fail inside @code{ranova}; its results agree with those of
    ## MATLAB's @code{anova} for the same contrasts.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.epsilon,
    ## RepeatedMeasuresModel.mauchly}
    ## @end deftypefn
    function [tbl, A, C, D] = ranova (this, varargin)

      [WM, args] = parsePairedArguments ({'WithinModel'}, {'separatemeans'}, ...
                                         varargin(:));
      if (! isempty (args))
        error ("RepeatedMeasuresModel.ranova: invalid optional paired argument.");
      endif
      k = numel (this.ResponseNames);
      [W, errmsg] = withinTerms (WM, this.WithinDesign, ...
                                 this.WithinFactorNames, k);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel.ranova: %s", errmsg);
      endif

      X = this.X_;
      E = this.R_' * this.R_;
      S = E / this.DFE;
      XtXi = inv (X' * X);
      nt = numel (this.TermNames_);
      A = cell (nt, 1);
      for i = 1:nt
        A{i} = eye (columns (X))(this.TermCols_{i},:);
      endfor

      names = {};
      vals = zeros (0, 8);
      for w = 1:numel (W.C)
        Q = orth (W.C{w});
        p = columns (Q);
        [gg, hf, lb] = sphericity (Q' * S * Q, p, this.DFE);
        SSE = trace (Q' * E * Q);
        dfE = p * this.DFE;
        for i = 1:nt
          AB = A{i} * this.B_;
          H = AB' * ((A{i} * XtXi * A{i}') \ AB);
          SSH = trace (Q' * H * Q);
          df = p * rows (A{i});
          F = (SSH / df) / (SSE / dfE);
          vals(end+1,:) = [SSH, df, SSH / df, F, ...
                           fcdf(F, df, dfE, 'upper'), ...
                           fcdf(F, df * gg, dfE * gg, 'upper'), ...
                           fcdf(F, df * hf, dfE * hf, 'upper'), ...
                           fcdf(F, df * lb, dfE * lb, 'upper')];
          if (isempty (W.names{w}))
            names{end+1} = this.TermNames_{i};
          else
            names{end+1} = [this.TermNames_{i}, ':', W.names{w}];
          endif
        endfor
        vals(end+1,:) = [SSE, dfE, SSE / dfE, NaN(1, 5)];
        if (isempty (W.names{w}))
          names{end+1} = 'Error';
        else
          names{end+1} = ['Error(', W.names{w}, ')'];
        endif
      endfor

      tbl = array2table (vals, 'VariableNames', {'SumSq', 'DF', 'MeanSq', ...
                         'F', 'pValue', 'pValueGG', 'pValueHF', 'pValueLB'}, ...
                         'RowNames', names);
      if (numel (W.C) == 1)
        C = W.C{1};
      else
        C = W.C;
      endif
      D = 0;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{tbl} =} epsilon (@var{rm})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{tbl} =} epsilon (@var{rm}, @var{C})
    ##
    ## Epsilon adjustments for repeated measures analysis of variance.
    ##
    ## @code{@var{tbl} = epsilon (@var{rm})} returns the corrections for
    ## departures from sphericity that @code{ranova} applies to its
    ## p-values, over the contrasts between successive responses: the
    ## variables @qcode{Uncorrected}, which is 1, @qcode{GreenhouseGeisser},
    ## @qcode{HuynhFeldt} and @qcode{LowerBound}.  The Huynh-Feldt value is
    ## Lecoutre's corrected form, capped at 1, and the lower bound is
    ## @math{1 / p} for @math{p} contrasts.
    ##
    ## @code{@var{tbl} = epsilon (@var{rm}, @var{C})} computes them over the
    ## contrast matrix @var{C}, which must have one row per response.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.ranova,
    ## RepeatedMeasuresModel.mauchly}
    ## @end deftypefn
    function tbl = epsilon (this, C)

      if (nargin < 2)
        C = [];
      endif
      Q = contrastBasis (this, C, 'epsilon');
      p = columns (Q);
      S = this.R_' * this.R_ / this.DFE;
      [gg, hf, lb] = sphericity (Q' * S * Q, p, this.DFE);
      tbl = array2table ([1, gg, hf, lb], 'VariableNames', {'Uncorrected', ...
                         'GreenhouseGeisser', 'HuynhFeldt', 'LowerBound'});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{tbl} =} mauchly (@var{rm})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{tbl} =} mauchly (@var{rm}, @var{C})
    ##
    ## Mauchly's test of sphericity.
    ##
    ## @code{@var{tbl} = mauchly (@var{rm})} tests whether the covariance of
    ## the contrasts between successive responses is a multiple of the
    ## identity, the sphericity that the uncorrected p-values of
    ## @code{ranova} assume.  @var{tbl} holds Mauchly's @qcode{W}, the
    ## chi-square statistic @qcode{ChiStat} with Bartlett's factor, its
    ## degrees of freedom @qcode{DF} and the @qcode{pValue}.  A singular
    ## covariance gives @var{W} 0, @var{ChiStat} @code{Inf} and a p-value of
    ## 0.
    ##
    ## @code{@var{tbl} = mauchly (@var{rm}, @var{C})} tests the contrasts of
    ## the matrix @var{C}, which must have one row per response.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.ranova,
    ## RepeatedMeasuresModel.epsilon}
    ## @end deftypefn
    function tbl = mauchly (this, C)

      if (nargin < 2)
        C = [];
      endif
      Q = contrastBasis (this, C, 'mauchly');
      p = columns (Q);
      S = Q' * (this.R_' * this.R_ / this.DFE) * Q;
      W = det (S) / (trace (S) / p) ^ p;
      df = p * (p + 1) / 2 - 1;
      ## Singularity is read from the conditioning, since rounding leaves the
      ## determinant of a singular matrix a small number of either sign
      if (rcond (S) < p * eps || ! (W > 0))
        W = 0;
        chi = Inf;
        pval = 0;
      else
        chi = -(this.DFE - (2 * p ^ 2 + p + 2) / (6 * p)) * log (W);
        pval = chi2cdf (chi, df, 'upper');
        if (df == 0)
          pval = 1;
        endif
      endif
      tbl = array2table ([W, chi, df, pval], 'VariableNames', ...
                         {'W', 'ChiStat', 'DF', 'pValue'});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{tbl} =} anova (@var{rm})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{tbl} =} anova (@var{rm}, @qcode{'WithinModel'}, @var{WM})
    ##
    ## Analysis of variance for between-subject effects.
    ##
    ## @code{@var{tbl} = anova (@var{rm})} tests the between-subject terms of
    ## the repeated measures model @var{rm} on the average of the repeated
    ## measures, one univariate analysis of variance.  @var{tbl} has the
    ## variables @qcode{Within}, the within-subject response analysed,
    ## @qcode{Between}, the between-subject term tested or @qcode{Error},
    ## @qcode{SumSq}, @qcode{DF}, @qcode{MeanSq}, @qcode{F} and
    ## @qcode{pValue}.  The intercept is named @qcode{constant}.
    ##
    ## @code{anova (@var{rm}, @qcode{'WithinModel'}, @var{WM})} analyses
    ## other responses built from the repeated measures, one block of rows
    ## each:
    ##
    ## @itemize
    ## @item @qcode{'separatemeans'}, the default: the average, named
    ## @qcode{Constant}.
    ## @item @qcode{'orthogonalcontrasts'}: the average and the orthogonal
    ## polynomial trends over a single numeric within-subject factor.
    ## @item A formula over the within-subject factors: each column of its
    ## effects-coded design, scaled to unit length and named after it.
    ## @item A contrast matrix with one row per response: each column as
    ## given, named @qcode{Contrast1}, @qcode{Contrast2}, @dots{}
    ## @end itemize
    ##
    ## MATLAB R2024a and R2026a fail on @qcode{'separatemeans'} given
    ## explicitly, though it is their documented default; here it is the
    ## default.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.ranova,
    ## RepeatedMeasuresModel.manova}
    ## @end deftypefn
    function tbl = anova (this, varargin)

      [WM, args] = parsePairedArguments ({'WithinModel'}, {'separatemeans'}, ...
                                         varargin(:));
      if (! isempty (args))
        error ("RepeatedMeasuresModel.anova: invalid optional paired argument.");
      endif
      k = numel (this.ResponseNames);
      [W, errmsg] = withinTerms (WM, this.WithinDesign, ...
                                 this.WithinFactorNames, k);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel.anova: %s", errmsg);
      endif

      ## One response per column, and its name
      switch (W.kind)
        case 'separatemeans'
          Cs = {ones(k, 1) / sqrt(k)};
          wn = {'Constant'};
        case 'matrix'
          Cs = num2cell (W.C{1}, 1);
          wn = arrayfun (@(i) sprintf ('Contrast%d', i), 1:numel (Cs), ...
                         'UniformOutput', false);
        case 'orthogonalcontrasts'
          Cs = [{ones(k, 1) / sqrt(k)}, W.C(2:end)];
          wn = W.cols;
        case 'formula'
          Cs = {};
          wn = {};
          for w = 1:numel (W.C)
            for j = 1:columns (W.C{w})
              c = W.C{w}(:,j);
              Cs{end+1} = c / norm (c);
              wn{end+1} = W.cols{w}{j};
            endfor
          endfor
      endswitch

      X = this.X_;
      E = this.R_' * this.R_;
      XtXi = inv (X' * X);
      bn = this.TermNames_;
      bn(strcmp (bn, '(Intercept)')) = {'constant'};
      within = {};
      between = {};
      vals = zeros (0, 5);
      for w = 1:numel (Cs)
        c = Cs{w};
        MSE = c' * E * c / this.DFE;
        for i = 1:numel (this.TermNames_)
          A = eye (columns (X))(this.TermCols_{i},:);
          Abc = A * this.B_ * c;
          SS = Abc' * ((A * XtXi * A') \ Abc);
          df = rows (A);
          F = (SS / df) / MSE;
          vals(end+1,:) = [SS, df, SS / df, F, fcdf(F, df, this.DFE, 'upper')];
          within{end+1} = wn{w};
          between{end+1} = bn{i};
        endfor
        vals(end+1,:) = [MSE * this.DFE, this.DFE, MSE, NaN, NaN];
        within{end+1} = wn{w};
        between{end+1} = 'Error';
      endfor

      tbl = table (categorical (within(:)), categorical (between(:)), ...
                   vals(:,1), vals(:,2), vals(:,3), vals(:,4), vals(:,5), ...
                   'VariableNames', {'Within', 'Between', 'SumSq', 'DF', ...
                   'MeanSq', 'F', 'pValue'});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{tbl} =} manova (@var{rm})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{tbl} =} manova (@var{rm}, @var{name}, @var{value})
    ## @deftypefnx {RepeatedMeasuresModel} {[@var{tbl}, @var{A}, @var{C}, @var{D}] =} manova (@dots{})
    ##
    ## Multivariate analysis of variance.
    ##
    ## @code{@var{tbl} = manova (@var{rm})} tests every term of the
    ## within-subject model of @var{rm}, as it was fitted, against every
    ## between-subject term, by the four multivariate statistics.  @var{tbl}
    ## has the variables @qcode{Within}, @qcode{Between}, @qcode{Statistic},
    ## one of @qcode{Pillai}, @qcode{Wilks}, @qcode{Hotelling} and
    ## @qcode{Roy}, its @qcode{Value}, the @qcode{F} statistic that
    ## approximates it, @qcode{RSquare}, the degrees of freedom @qcode{df1}
    ## and @qcode{df2}, and the @qcode{pValue}.  Under
    ## @qcode{'separatemeans'} the within-subject hypothesis is that of equal
    ## means, named @qcode{Constant}, and a contrast matrix is named
    ## @qcode{Specified contrast}.
    ##
    ## The name-value arguments are @qcode{'WithinModel'}, which takes the
    ## forms @code{fitrm} takes but @qcode{'orthogonalcontrasts'}, and
    ## @qcode{'By'}, the name of a between-subject factor, which tests the
    ## within-subject hypotheses at each of its levels in place of the
    ## between-subject terms.
    ##
    ## @code{[@var{tbl}, @var{A}, @var{C}, @var{D}] = manova (@dots{})} also
    ## returns the hypotheses as @math{A B C = D}: @var{A} a cell column of
    ## the between-subject hypothesis matrices, @var{C} the within-subject
    ## contrast, a cell row of one per term when there are several, and
    ## @var{D} zero.
    ##
    ## Pillai's trace, Wilks' lambda and Roy's root are approximated by F as
    ## in MATLAB.  The Hotelling-Lawley trace uses McKeon's F approximation
    ## with its own second degrees of freedom, as SAS does, and the
    ## Pillai-Samson one where McKeon's is undefined.  MATLAB R2024a computes
    ## McKeon's F but refers it to the Pillai-Samson degrees of freedom,
    ## which makes its p-values too small.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.coeftest,
    ## RepeatedMeasuresModel.ranova}
    ## @end deftypefn
    function [tbl, A, C, D] = manova (this, varargin)

      [WM, By, args] = parsePairedArguments ({'WithinModel', 'By'}, ...
                                             {this.WithinModel, ''}, ...
                                             varargin(:));
      if (! isempty (args))
        error ("RepeatedMeasuresModel.manova: invalid optional paired argument.");
      endif
      if ((ischar (WM) || isa (WM, 'string')) ...
          && strcmpi (char (WM), 'orthogonalcontrasts'))
        error (strcat ("RepeatedMeasuresModel.manova: the", ...
                       " 'orthogonalcontrasts' model cannot be used with", ...
                       " manova."));
      endif
      k = numel (this.ResponseNames);
      [W, errmsg] = withinTerms (WM, this.WithinDesign, ...
                                 this.WithinFactorNames, k);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel.manova: %s", errmsg);
      endif
      switch (W.kind)
        case 'separatemeans'
          wn = {'Constant'};
        case 'matrix'
          wn = {'Specified contrast'};
        otherwise
          wn = W.names;
          wn(cellfun (@isempty, wn)) = {'(Intercept)'};
      endswitch

      ## The between-subject hypotheses: each term, or each level of BY
      X = this.X_;
      nc = columns (X);
      if (isempty (By))
        A = cell (numel (this.TermNames_), 1);
        for i = 1:numel (A)
          A{i} = eye (nc)(this.TermCols_{i},:);
        endfor
        bn = this.TermNames_;
      else
        [A, bn, errmsg] = byHypotheses (this, By);
        if (! isempty (errmsg))
          error ("RepeatedMeasuresModel.manova: %s", errmsg);
        endif
      endif

      stat = {'Pillai'; 'Wilks'; 'Hotelling'; 'Roy'};
      within = {};
      between = {};
      vals = zeros (0, 6);
      for w = 1:numel (W.C)
        for i = 1:numel (A)
          vals = [vals; mvtests(this, A{i}, W.C{w}, 0)];
          within = [within; repmat(wn(w), 4, 1)];
          between = [between; repmat(bn(i), 4, 1)];
        endfor
      endfor
      tbl = table (categorical (within), categorical (between), ...
                   categorical (repmat (stat, rows (vals) / 4, 1)), ...
                   vals(:,1), vals(:,2), vals(:,3), vals(:,4), vals(:,5), ...
                   vals(:,6), 'VariableNames', {'Within', 'Between', ...
                   'Statistic', 'Value', 'F', 'RSquare', 'df1', 'df2', ...
                   'pValue'});
      if (numel (W.C) == 1)
        C = W.C{1};
      else
        C = W.C;
      endif
      D = 0;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{tbl} =} coeftest (@var{rm}, @var{A}, @var{C})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{tbl} =} coeftest (@var{rm}, @var{A}, @var{C}, @var{D})
    ##
    ## Linear hypothesis test on the coefficients of a repeated measures
    ## model.
    ##
    ## @code{@var{tbl} = coeftest (@var{rm}, @var{A}, @var{C})} tests the
    ## hypothesis @math{A B C = 0} on the coefficient matrix @math{B} of
    ## @var{rm}, @var{A} having one column per between-subject coefficient
    ## and @var{C} one row per response.  @var{tbl} holds the four
    ## multivariate statistics as @code{manova} reports them: the variables
    ## @qcode{Statistic}, @qcode{Value}, @qcode{F}, @qcode{RSquare},
    ## @qcode{df1}, @qcode{df2} and @qcode{pValue}.
    ##
    ## @code{@var{tbl} = coeftest (@var{rm}, @var{A}, @var{C}, @var{D})} tests
    ## @math{A B C = D} instead, @var{D} a scalar or a matrix with as many
    ## rows as @var{A} and as many columns as @var{C}.  The default is 0.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.manova}
    ## @end deftypefn
    function tbl = coeftest (this, A, C, D)

      if (nargin < 3)
        error ("RepeatedMeasuresModel.coeftest: too few input arguments.");
      endif
      nc = columns (this.X_);
      k = numel (this.ResponseNames);
      if (! (isnumeric (A) && isreal (A) && ismatrix (A) && columns (A) == nc ...
             && ! isempty (A)))
        error ("RepeatedMeasuresModel.coeftest: A must be a matrix with %d columns.", ...
               nc);
      endif
      if (! (isnumeric (C) && isreal (C) && ismatrix (C) && rows (C) == k ...
             && ! isempty (C)))
        error ("RepeatedMeasuresModel.coeftest: C must be a matrix with %d rows.", ...
               k);
      endif
      if (nargin < 4)
        D = 0;
      endif
      if (! (isnumeric (D) && isreal (D) && (isscalar (D) ...
             || isequal (size (D), [rows(A), columns(C)]))))
        error (strcat ("RepeatedMeasuresModel.coeftest: D must be a scalar", ...
                       " or a %d-by-%d matrix."), rows (A), columns (C));
      endif
      vals = mvtests (this, double (A), double (C), double (D));
      tbl = table (categorical ({'Pillai'; 'Wilks'; 'Hotelling'; 'Roy'}), ...
                   vals(:,1), vals(:,2), vals(:,3), vals(:,4), vals(:,5), ...
                   vals(:,6), 'VariableNames', {'Statistic', 'Value', 'F', ...
                   'RSquare', 'df1', 'df2', 'pValue'});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{tbl} =} margmean (@var{rm}, @var{vars})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{tbl} =} margmean (@var{rm}, @var{vars}, @qcode{'Alpha'}, @var{alpha})
    ##
    ## Estimated marginal means.
    ##
    ## @code{@var{tbl} = margmean (@var{rm}, @var{vars})} estimates the mean
    ## of the repeated measures at each combination of the levels of the
    ## factors named by @var{vars}, a character vector or a cell array of
    ## them, each a categorical between-subject factor or a within-subject
    ## factor.  The other between-subject factors are averaged over their
    ## levels with equal weights, a continuous predictor is held at its mean
    ## over the subjects the fit used, and the responses are averaged over
    ## the other within-subject factors.  @var{tbl} holds one variable per
    ## factor, the first varying slowest, and @qcode{Mean}, @qcode{StdErr},
    ## @qcode{Lower} and @qcode{Upper}, the limits of a @math{100 (1 -
    ## alpha)} per cent confidence interval on @code{DFE} degrees of freedom.
    ## The default @var{alpha} is 0.05.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.multcompare,
    ## RepeatedMeasuresModel.grpstats}
    ## @end deftypefn
    function tbl = margmean (this, vars, varargin)

      if (nargin < 2)
        error ("RepeatedMeasuresModel.margmean: too few input arguments.");
      endif
      [alpha, args] = parsePairedArguments ({'Alpha'}, {0.05}, varargin(:));
      if (! isempty (args))
        error ("RepeatedMeasuresModel.margmean: invalid optional paired argument.");
      endif
      checkAlpha (alpha, 'margmean');
      [F, errmsg] = factorInfo (this, vars);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel.margmean: %s", errmsg);
      endif
      L = combos (F);
      S = this.R_' * this.R_ / this.DFE;
      V = inv (this.X_' * this.X_);
      est = zeros (rows (L), 4);
      tq = tinv (1 - alpha / 2, this.DFE);
      for r = 1:rows (L)
        [a, w] = weights (this, F, L(r,:));
        m = a * this.B_ * w;
        se = sqrt ((a * V * a') * (w' * S * w));
        est(r,:) = [m, se, m - tq * se, m + tq * se];
      endfor
      tbl = levelTable (F, L, {});
      tbl = [tbl, array2table(est, 'VariableNames', ...
                              {'Mean', 'StdErr', 'Lower', 'Upper'})];

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{tbl} =} grpstats (@var{rm}, @var{g})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{tbl} =} grpstats (@var{rm}, @var{g}, @var{stats})
    ##
    ## Descriptive statistics of the repeated measures by group.
    ##
    ## @code{@var{tbl} = grpstats (@var{rm}, @var{g})} pools the repeated
    ## measures of every subject and summarizes them at each combination of
    ## the levels of the factors named by @var{g}, a character vector or a
    ## cell array of them, each a categorical between-subject factor or a
    ## within-subject factor.  @var{tbl} holds one variable per factor,
    ## @qcode{GroupCount}, the number of values in the group, and their
    ## @qcode{mean} and @qcode{std}.
    ##
    ## @code{@var{tbl} = grpstats (@var{rm}, @var{g}, @var{stats})} computes
    ## the statistics named by @var{stats} instead, a character vector or a
    ## cell array of them, each one of @qcode{'mean'}, @qcode{'median'},
    ## @qcode{'std'}, @qcode{'var'}, @qcode{'sem'}, @qcode{'min'},
    ## @qcode{'max'}, @qcode{'range'} and @qcode{'numel'}.  Missing values are
    ## left out.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.margmean}
    ## @end deftypefn
    function tbl = grpstats (this, g, stats)

      if (nargin < 2)
        error ("RepeatedMeasuresModel.grpstats: too few input arguments.");
      endif
      if (nargin < 3)
        stats = {'mean', 'std'};
      endif
      if (ischar (stats) || (isa (stats, 'string') && isscalar (stats)))
        stats = cellstr (stats);
      endif
      known = {'mean', 'median', 'std', 'var', 'sem', 'min', 'max', ...
               'range', 'numel'};
      if (! iscellstr (stats) || ! all (ismember (stats, known)))
        error (strcat ("RepeatedMeasuresModel.grpstats: STATS must name", ...
                       " statistics among 'mean', 'median', 'std',", ...
                       " 'var', 'sem', 'min', 'max', 'range' and", ...
                       " 'numel'."));
      endif
      [F, errmsg] = factorInfo (this, g);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel.grpstats: %s", errmsg);
      endif
      L = combos (F);
      Y = zeros (rows (this.BetweenDesign), numel (this.ResponseNames));
      for j = 1:columns (Y)
        Y(:,j) = double (this.BetweenDesign.(this.ResponseNames{j}));
      endfor
      keep = true (rows (L), 1);
      vals = zeros (rows (L), numel (stats) + 1);
      for r = 1:rows (L)
        [rows_, cols_] = groupMembers (this, F, L(r,:));
        y = Y(rows_, cols_)(:);
        y = y(! isnan (y));
        keep(r) = ! isempty (y);
        vals(r,1) = numel (y);
        for i = 1:numel (stats)
          vals(r,i+1) = groupStatistic (y, stats{i});
        endfor
      endfor
      tbl = levelTable (F, L(keep,:), {});
      tbl = [tbl, array2table(vals(keep,:), 'VariableNames', ...
                              [{'GroupCount'}, stats(:)'])];

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{tbl} =} multcompare (@var{rm}, @var{var})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{tbl} =} multcompare (@var{rm}, @var{var}, @var{name}, @var{value})
    ##
    ## Multiple comparison of estimated marginal means.
    ##
    ## @code{@var{tbl} = multcompare (@var{rm}, @var{var})} compares the
    ## estimated marginal means of every pair of levels of the factor
    ## @var{var}, a categorical between-subject factor or a within-subject
    ## factor, as @code{margmean} estimates them.  @var{tbl} holds one row per
    ## ordered pair, with the variables @qcode{@var{var}_1} and
    ## @qcode{@var{var}_2}, @qcode{Difference}, @qcode{StdErr}, the adjusted
    ## @qcode{pValue}, and @qcode{Lower} and @qcode{Upper}, the limits of the
    ## simultaneous confidence interval.
    ##
    ## The name-value arguments are:
    ##
    ## @multitable @columnfractions 0.25 0.75
    ## @headitem Name @tab Value
    ## @item @qcode{'By'} @tab The name of another factor; the comparisons are
    ## made at each of its levels, which become the first variable of
    ## @var{tbl}.
    ##
    ## @item @qcode{'ComparisonType'} @tab @qcode{'tukey-kramer'}, the
    ## default, @qcode{'bonferroni'}, @qcode{'dunn-sidak'}, @qcode{'lsd'} or
    ## @qcode{'scheffe'}.
    ##
    ## @item @qcode{'Alpha'} @tab The significance level of the intervals.
    ## The default is 0.05.
    ## @end multitable
    ##
    ## The Tukey-Kramer p-values and limits come from the studentized range
    ## distribution, @code{stdrcdf} and @code{stdrinv}.  MATLAB R2024a and
    ## R2026a floor their Tukey-Kramer p-values, 9.56e-10 over three groups
    ## whatever the difference, and R2024a's Dunn-Sidak p-values read 0 where
    ## the Bonferroni ones are near 1e-36; both are computed here.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.margmean, stdrcdf}
    ## @end deftypefn
    function tbl = multcompare (this, var, varargin)

      if (nargin < 2)
        error ("RepeatedMeasuresModel.multcompare: too few input arguments.");
      endif
      [by, ctype, alpha, args] = parsePairedArguments ( ...
                {'By', 'ComparisonType', 'Alpha'}, ...
                {'', 'tukey-kramer', 0.05}, varargin(:));
      if (! isempty (args))
        error ("RepeatedMeasuresModel.multcompare: invalid optional paired argument.");
      endif
      checkAlpha (alpha, 'multcompare');
      types = {'tukey-kramer', 'bonferroni', 'dunn-sidak', 'lsd', 'scheffe'};
      if (! ((ischar (ctype) || isa (ctype, 'string')) ...
             && any (strcmpi (char (ctype), types))))
        error (strcat ("RepeatedMeasuresModel.multcompare: 'ComparisonType'", ...
                       " must be 'tukey-kramer', 'bonferroni',", ...
                       " 'dunn-sidak', 'lsd' or 'scheffe'."));
      endif
      ctype = lower (char (ctype));
      if (isempty (by))
        names = {var};
      else
        names = {by, var};
        if (strcmp (char (by), char (var)))
          error (strcat ("RepeatedMeasuresModel.multcompare: 'By' must", ...
                         " differ from the factor compared."));
        endif
      endif
      [F, errmsg] = factorInfo (this, names);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel.multcompare: %s", errmsg);
      endif
      nb = 1;
      if (numel (F) == 2)
        nb = numel (F(1).levels);
      endif
      m = numel (F(end).levels);
      npairs = m * (m - 1) / 2;
      S = this.R_' * this.R_ / this.DFE;
      V = inv (this.X_' * this.X_);
      crit = critical (m, npairs, this.DFE, ctype, alpha);
      L = zeros (0, numel (F) + 1);
      vals = zeros (0, 5);
      for b = 1:nb
        ## A pair and its reverse share their p-value
        pp = NaN (m);
        for i = 1:m
          for j = [1:i-1, i+1:m]
            if (numel (F) == 2)
              [ai, wi] = weights (this, F, [b, i]);
              [aj, wj] = weights (this, F, [b, j]);
            else
              [ai, wi] = weights (this, F, i);
              [aj, wj] = weights (this, F, j);
            endif
            d = ai * this.B_ * wi - aj * this.B_ * wj;
            if (F(end).between)
              se = sqrt (((ai - aj) * V * (ai - aj)') * (wi' * S * wi));
            else
              se = sqrt ((ai * V * ai') * ((wi - wj)' * S * (wi - wj)));
            endif
            if (isnan (pp(j,i)))
              pp(i,j) = pairwise (abs (d) / se, m, npairs, this.DFE, ctype);
            else
              pp(i,j) = pp(j,i);
            endif
            pv = pp(i,j);
            vals(end+1,:) = [d, se, pv, d - crit * se, d + crit * se];
            if (numel (F) == 2)
              L(end+1,:) = [b, i, j];
            else
              L(end+1,:) = [i, j];
            endif
          endfor
        endfor
      endfor
      if (numel (F) == 1)
        G = [F, F];
      else
        G = [F, F(2)];
      endif
      vn = {sprintf('%s_1', F(end).name), sprintf('%s_2', F(end).name)};
      if (numel (F) == 2)
        vn = [{F(1).name}, vn];
      endif
      tbl = levelTable (G, L, vn);
      tbl = [tbl, array2table(vals, 'VariableNames', ...
             {'Difference', 'StdErr', 'pValue', 'Lower', 'Upper'})];

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{ypred} =} predict (@var{rm})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{ypred} =} predict (@var{rm}, @var{tnew})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{ypred} =} predict (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {RepeatedMeasuresModel} {[@var{ypred}, @var{yci}] =} predict (@dots{})
    ##
    ## Predict the repeated measures.
    ##
    ## @code{@var{ypred} = predict (@var{rm}, @var{tnew})} returns the
    ## predicted repeated measures of the subjects in the table @var{tnew},
    ## one row each, from their between-subject predictors.  @var{tnew}
    ## defaults to @code{BetweenDesign}.  A subject missing a predictor, or
    ## holding a level the model was not fitted on, is predicted as
    ## @code{NaN}.
    ##
    ## Under @qcode{'separatemeans'} the prediction is one column per
    ## response.  Under any other within-subject model the means are smoothed
    ## through that model: projected onto the terms of a formula, or onto the
    ## polynomial of @qcode{'orthogonalcontrasts'}, and evaluated at the
    ## within-subject design, which may then be a new one.
    ##
    ## The name-value arguments are @qcode{'WithinDesign'}, the within-subject
    ## design to predict at, as @code{fitrm} takes it; @qcode{'WithinModel'},
    ## the within-subject model, by default the one fitted; and
    ## @qcode{'Alpha'}, the significance level of the intervals, 0.05 by
    ## default.  Under @qcode{'separatemeans'} a new within-subject design is
    ## ignored with a warning, as in MATLAB.
    ##
    ## @code{[@var{ypred}, @var{yci}] = predict (@dots{})} also returns the
    ## confidence limits of the predicted means as an array of the size of
    ## @var{ypred} by 2, the lower limits first.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.random}
    ## @end deftypefn
    function [ypred, yci] = predict (this, varargin)

      tnew = this.BetweenDesign;
      if (numel (varargin) > 0 && istable (varargin{1}))
        tnew = varargin{1};
        varargin(1) = [];
      endif
      [WDn, WM, alpha, args] = parsePairedArguments ( ...
                {'WithinDesign', 'WithinModel', 'Alpha'}, ...
                {[], this.WithinModel, 0.05}, varargin(:));
      if (! isempty (args))
        error ("RepeatedMeasuresModel.predict: invalid optional paired argument.");
      endif
      checkAlpha (alpha, 'predict');
      [Xn, errmsg] = newDesign (this, tnew);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel.predict: %s", errmsg);
      endif
      [P, errmsg] = smoother (this, WM, WDn);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel.predict: %s", errmsg);
      endif
      ypred = Xn * this.B_ * P;
      if (nargout > 1)
        S = this.R_' * this.R_ / this.DFE;
        vx = sum ((Xn / (this.X_' * this.X_)) .* Xn, 2);
        se = sqrt (vx .* sum (P .* (S * P), 1));
        tq = tinv (1 - alpha / 2, this.DFE);
        yci = cat (3, ypred - tq * se, ypred + tq * se);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RepeatedMeasuresModel} {@var{ysim} =} random (@var{rm})
    ## @deftypefnx {RepeatedMeasuresModel} {@var{ysim} =} random (@var{rm}, @var{tnew})
    ##
    ## Generate new random repeated measures.
    ##
    ## @code{@var{ysim} = random (@var{rm}, @var{tnew})} draws one set of
    ## repeated measures for each subject of the table @var{tnew}, from a
    ## normal distribution with the subject's predicted means and the
    ## estimated covariance @code{Covariance}.  @var{tnew} defaults to
    ## @code{BetweenDesign}.  A subject missing a predictor gets a row of
    ## @code{NaN}.
    ##
    ## @seealso{fitrm, RepeatedMeasuresModel.predict, mvnrnd}
    ## @end deftypefn
    function ysim = random (this, tnew)

      if (nargin < 2)
        tnew = this.BetweenDesign;
      elseif (! istable (tnew))
        error ("RepeatedMeasuresModel.random: TNEW must be a table.");
      endif
      [Xn, errmsg] = newDesign (this, tnew);
      if (! isempty (errmsg))
        error ("RepeatedMeasuresModel.random: %s", errmsg);
      endif
      S = this.R_' * this.R_ / this.DFE;
      mu = Xn * this.B_;
      ysim = NaN (size (mu));
      ok = all (isfinite (mu), 2);
      if (any (ok))
        ysim(ok,:) = mvnrnd (mu(ok,:), S);
      endif

    endfunction

  endmethods

  methods (Access = private)

    ## An orthonormal basis of the contrast C, or of the contrasts between
    ## successive responses when C is empty.
    function Q = contrastBasis (this, C, meth)
      k = numel (this.ResponseNames);
      if (isempty (C))
        C = successive (k);
      elseif (! (isnumeric (C) && isreal (C) && ismatrix (C) && rows (C) == k))
        error ("RepeatedMeasuresModel.%s: C must be a matrix with %d rows.", ...
               meth, k);
      endif
      Q = orth (double (C));
    endfunction

    ## The factors named by NAMES, each categorical and either between or
    ## within: its name, side, level keys and native values.
    function [F, errmsg] = factorInfo (this, names)
      F = struct ('name', {}, 'between', {}, 'levels', {}, 'native', {});
      errmsg = '';
      if (ischar (names) || (isa (names, 'string') && isscalar (names)))
        names = cellstr (names);
      endif
      if (! iscellstr (names) || isempty (names))
        errmsg = "factor names must be text.";
        return;
      endif
      for i = 1:numel (names)
        nm = names{i};
        if (any (strcmp (nm, this.BetweenFactorNames)) ...
            && this.Coding_.(nm).iscat)
          c = this.Coding_.(nm);
          F(end+1) = struct ('name', nm, 'between', true, ...
                             'levels', {c.levels}, 'native', {c.native});
        elseif (any (strcmp (nm, this.WithinFactorNames)))
          [lev, ~, native, iscat] = levelsOf (this.WithinDesign.(nm));
          if (! iscat)
            v = this.WithinDesign.(nm);
            u = unique (v(:));
            lev = arrayfun (@num2str, u, 'UniformOutput', false);
            native = num2cell (u);
          endif
          F(end+1) = struct ('name', nm, 'between', false, ...
                             'levels', {lev}, 'native', {native});
        else
          errmsg = sprintf ("'%s' is not a categorical factor of the model.", nm);
          return;
        endif
      endfor
    endfunction

    ## The between-subject row A and the within-subject weights W at the
    ## levels IDX of the factors F, the other factors averaged.
    function [a, w] = weights (this, F, idx)
      fixed = struct ();
      k = numel (this.ResponseNames);
      w = ones (k, 1);
      for i = 1:numel (F)
        if (F(i).between)
          fixed.(F(i).name) = idx(i);
        else
          v = this.WithinDesign.(F(i).name);
          if (isnumeric (v))
            w &= (double (v(:)) == F(i).native{idx(i)});
          else
            w &= strcmp (cellstr (v(:)), F(i).levels{idx(i)});
          endif
        endif
      endfor
      w = double (w) / sum (w);
      a = zeros (1, columns (this.X_));
      off = numel (this.TermNames_) - numel (this.TermVars_);
      if (off)
        a(1) = 1;
      endif
      for t = 1:numel (this.TermVars_)
        M = 1;
        for v = this.TermVars_{t}
          c = this.Coding_.(v{1});
          if (! c.iscat)
            e = c.mean;
          elseif (isfield (fixed, v{1}))
            L = numel (c.levels);
            e = zeros (1, L - 1);
            l = fixed.(v{1});
            if (l < L)
              e(l) = 1;
            else
              e(:) = -1;
            endif
          else
            e = zeros (1, numel (c.levels) - 1);
          endif
          M = kron (M, e);
        endfor
        a(this.TermCols_{t+off}) = M;
      endfor
    endfunction

    ## The subjects and responses at the levels IDX of the factors F.
    function [r, c] = groupMembers (this, F, idx)
      r = true (rows (this.BetweenDesign), 1);
      c = true (numel (this.ResponseNames), 1);
      for i = 1:numel (F)
        if (F(i).between)
          [~, ic] = levelsOf (this.BetweenDesign.(F(i).name), F(i).levels);
          r &= (ic(:) == idx(i));
        else
          v = this.WithinDesign.(F(i).name);
          if (isnumeric (v))
            c &= (double (v(:)) == F(i).native{idx(i)});
          else
            c &= strcmp (cellstr (v(:)), F(i).levels{idx(i)});
          endif
        endif
      endfor
    endfunction

    ## The effects-coded design of the subjects of TNEW, NaN on a row missing
    ## a predictor or holding a level the model does not know.
    function [X, errmsg] = newDesign (this, tnew)
      X = [];
      errmsg = '';
      if (! istable (tnew))
        errmsg = "TNEW must be a table.";
        return;
      endif
      miss = setdiff (this.BetweenFactorNames, tnew.Properties.VariableNames);
      if (! isempty (miss))
        errmsg = sprintf ("TNEW has no variable '%s'.", miss{1});
        return;
      endif
      intercept = numel (this.TermNames_) > numel (this.TermVars_);
      [X, ~, ~, ~, bad] = effectsDesign (tnew, this.TermVars_, intercept, ...
                                         this.BetweenFactorNames, this.Coding_);
      X(bad,:) = NaN;
    endfunction

    ## The matrix that maps fitted means over the responses onto the
    ## predictions at the within-subject design WDN under the model WM.
    function [P, errmsg] = smoother (this, WM, WDn)
      errmsg = '';
      k = numel (this.ResponseNames);
      P = [];
      if (isa (WM, 'string') && isscalar (WM))
        WM = char (WM);
      endif
      if (isnumeric (WM) || strcmpi (WM, 'separatemeans'))
        if (! isempty (WDn))
          warning (strcat ("RepeatedMeasuresModel.predict: the", ...
                           " 'separatemeans' model does not use the", ...
                           " 'WithinDesign' given."));
        endif
        P = eye (k);
        return;
      endif
      WD = this.WithinDesign;
      wn = this.WithinFactorNames;
      if (isempty (WDn))
        WDn = WD;
      elseif (isnumeric (WDn) && isvector (WDn) && numel (wn) == 1)
        WDn = table (double (WDn(:)), 'VariableNames', wn);
      elseif (! istable (WDn) || ! all (ismember (wn, WDn.Properties.VariableNames)))
        errmsg = "'WithinDesign' must hold the within-subject factors.";
        return;
      endif
      [~, errmsg] = withinTerms (WM, WD, wn, k);
      if (! isempty (errmsg))
        return;
      endif
      if (strcmpi (WM, 'orthogonalcontrasts'))
        v = double (WD.(wn{1}));
        mu = mean (v);
        sc = max (1, std (v));
        Wd = ((v(:) - mu) / sc) .^ (0:k-1);
        Wn = ((double (WDn.(wn{1}))(:) - mu) / sc) .^ (0:k-1);
      else
        s = parseWilkinsonFormula (['~', WM]);
        model = s.model;
        empty = cellfun (@isempty, model);
        terms = cellfun (@cellstr, model(! empty), 'UniformOutput', false);
        vars = unique ([terms{:}], 'stable');
        intercept = any (empty) || isempty (model);
        [Wd, ~, ~, ~, ~, coding] = effectsDesign (WD, terms, intercept, vars);
        [Wn, ~, ~, ~, bad] = effectsDesign (WDn, terms, intercept, vars, coding);
        Wn(bad,:) = NaN;
      endif
      P = (Wn * pinv (Wd))';
    endfunction

    ## The four multivariate statistics of A B C = D, one row each of Value,
    ## F, RSquare, df1, df2 and pValue.
    function vals = mvtests (this, A, C, D)
      X = this.X_;
      Q = orth (C);
      T = C \ Q;
      E = Q' * (this.R_' * this.R_) * Q;
      M = (A * this.B_ * C - D) * T;
      H = M' * ((A / (X' * X) * A') \ M);
      lam = sort (max (real (eig (E \ H)), 0), 'descend');
      vals = multivariateF (lam, columns (Q), rank (A), this.DFE);
    endfunction

    ## The hypotheses at each level of the between-subject factor BY: the
    ## intercept and BY's own columns at that level's effects code.
    function [A, names, errmsg] = byHypotheses (this, by)
      A = {};
      names = {};
      errmsg = '';
      by = char (by);
      if (! any (strcmp (by, this.BetweenFactorNames)))
        errmsg = sprintf ("'By' must name a between-subject factor, not '%s'.", ...
                          by);
        return;
      endif
      v = this.BetweenDesign.(by);
      if (iscategorical (v))
        lev = categories (v);
        lev = lev(ismember (lev, cellstr (v)));
      elseif (iscellstr (v) || isa (v, 'string') || ischar (v))
        lev = unique (cellstr (v));
      else
        errmsg = sprintf ("'By' must name a categorical factor, not '%s'.", by);
        return;
      endif
      col = find (strcmp (this.TermNames_, by));
      L = numel (lev);
      for l = 1:L
        a = zeros (1, columns (this.X_));
        a(1) = 1;
        e = zeros (1, L - 1);
        if (l < L)
          e(l) = 1;
        else
          e(:) = -1;
        endif
        a(this.TermCols_{col}) = e;
        A{end+1,1} = a;
        names{end+1} = sprintf ('%s=%s', by, lev{l});
      endfor
    endfunction

  endmethods

endclassdef

## The contrasts between successive responses, one column per pair.
function C = successive (k)
  C = [eye(k - 1); zeros(1, k - 1)] - [zeros(1, k - 1); eye(k - 1)];
endfunction

## Greenhouse-Geisser, Huynh-Feldt and lower bound corrections from the
## covariance S of P orthonormal contrasts on DFE error degrees of freedom.
## The Huynh-Feldt value is Lecoutre's corrected form, capped at 1.
function [gg, hf, lb] = sphericity (S, p, dfe)
  gg = trace (S) ^ 2 / (p * trace (S ^ 2));
  hf = min (1, ((dfe + 1) * p * gg - 2) / (p * (dfe - p * gg)));
  lb = 1 / p;
endfunction

## Every combination of the levels of the factors F, the first varying
## slowest, one row of level indices each.
function L = combos (F)
  n = cellfun (@numel, {F.levels});
  L = zeros (prod (n), numel (n));
  for i = 1:numel (n)
    reps = prod (n(i+1:end));
    L(:,i) = repmat (kron ((1:n(i))', ones (reps, 1)), prod (n(1:i-1)), 1);
  endfor
endfunction

## A table of the levels L of the factors F, one variable per column of L,
## named NAMES or after the factors, in the type each factor was given in.
function tbl = levelTable (F, L, names)
  if (isempty (names))
    names = {F.name};
  endif
  cols = cell (1, columns (L));
  for i = 1:columns (L)
    nat = F(i).native(L(:,i));
    if (iscell (nat{1}))
      cols{i} = [nat{:}]';
    else
      cols{i} = vertcat (nat{:});
    endif
  endfor
  tbl = table (cols{:}, 'VariableNames', names);
endfunction

## One statistic of the values Y.
function v = groupStatistic (y, name)
  switch (name)
    case 'mean'
      v = mean (y);
    case 'median'
      v = median (y);
    case 'std'
      v = std (y);
    case 'var'
      v = var (y);
    case 'sem'
      v = std (y) / sqrt (numel (y));
    case 'min'
      v = min (y);
    case 'max'
      v = max (y);
    case 'range'
      v = max (y) - min (y);
    case 'numel'
      v = numel (y);
  endswitch
endfunction

## The adjusted p-value of a pairwise comparison with |t| = T among M means,
## NPAIRS pairs and DFE error degrees of freedom.  The two-sided t tail is the
## incomplete beta function, which keeps its digits far out where tcdf
## rounds to zero.
function p = pairwise (T, m, npairs, dfe, ctype)
  p2 = betainc (T ^ 2 / (dfe + T ^ 2), 1/2, dfe / 2, 'upper');
  switch (ctype)
    case 'tukey-kramer'
      p = stdrcdf (sqrt (2) * T, m, dfe, 'upper');
    case 'bonferroni'
      p = min (1, npairs * p2);
    case 'dunn-sidak'
      p = - expm1 (npairs * log1p (- p2));
    case 'lsd'
      p = p2;
    case 'scheffe'
      p = fcdf (T ^ 2 / (m - 1), m - 1, dfe, 'upper');
  endswitch
endfunction

## The multiple of the standard error that gives simultaneous limits at
## level ALPHA for the comparison type CTYPE.
function crit = critical (m, npairs, dfe, ctype, alpha)
  switch (ctype)
    case 'tukey-kramer'
      crit = stdrinv (1 - alpha, m, dfe) / sqrt (2);
    case 'bonferroni'
      crit = tinv (1 - alpha / (2 * npairs), dfe);
    case 'dunn-sidak'
      crit = tinv (1 - (1 - (1 - alpha) ^ (1 / npairs)) / 2, dfe);
    case 'lsd'
      crit = tinv (1 - alpha / 2, dfe);
    case 'scheffe'
      crit = sqrt ((m - 1) * finv (1 - alpha, m - 1, dfe));
  endswitch
endfunction

## Refuse an ALPHA outside (0, 1).
function checkAlpha (alpha, meth)
  if (! (isnumeric (alpha) && isscalar (alpha) && alpha > 0 && alpha < 1))
    error ("RepeatedMeasuresModel.%s: 'Alpha' must be a scalar in (0, 1).", ...
           meth);
  endif
endfunction

## Pillai's trace, Wilks' lambda, the Hotelling-Lawley trace and Roy's root
## from the eigenvalues LAM of inv (E) * H, for P contrasts, a hypothesis of
## rank Q and DFE error degrees of freedom: one row each of Value, F,
## RSquare, df1, df2 and pValue.  The Hotelling-Lawley trace uses McKeon's F
## approximation where it is defined and the Pillai-Samson one elsewhere,
## which is exact for a single nonzero eigenvalue.
function vals = multivariateF (lam, p, q, dfe)
  s = min (p, q);
  m = (abs (p - q) - 1) / 2;
  n = (dfe - p - 1) / 2;
  vals = zeros (4, 6);

  V = sum (lam ./ (1 + lam));
  r2 = V / s;
  d1 = s * (2 * m + s + 1);
  d2 = s * (2 * n + s + 1);
  vals(1,:) = [V, d2 / d1 * r2 / (1 - r2), r2, d1, d2, 0];

  L = prod (1 ./ (1 + lam));
  if (p ^ 2 + q ^ 2 - 5 > 0)
    t = sqrt ((p ^ 2 * q ^ 2 - 4) / (p ^ 2 + q ^ 2 - 5));
  else
    t = 1;
  endif
  r2 = 1 - L ^ (1 / t);
  d1 = p * q;
  d2 = (dfe - (p - q + 1) / 2) * t - (p * q - 2) / 2;
  vals(2,:) = [L, d2 / d1 * r2 / (1 - r2), r2, d1, d2, 0];

  T = sum (lam);
  r2 = (T / s) / (1 + T / s);
  if (s > 1 && n > 1)
    b = (p + 2 * n) * (q + 2 * n) / (2 * (2 * n + 1) * (n - 1));
    d2 = 4 + (p * q + 2) / (b - 1);
    c = p * q * (d2 - 2) / (d2 * (dfe - p - 1));
    d1 = p * q;
    F = T / c;
  else
    d1 = s * (2 * m + s + 1);
    d2 = 2 * (s * n + 1);
    F = d2 * T / (s * d1);
  endif
  vals(3,:) = [T, F, r2, d1, d2, 0];

  th = max ([lam; 0]);
  rr = max (p, q);
  d2 = dfe - rr + q;
  vals(4,:) = [th, th * d2 / rr, th / (1 + th), rr, d2, 0];

  vals(:,6) = fcdf (vals(:,2), vals(:,4), vals(:,5), 'upper');
endfunction

## Response names from the left side of the formula: ranges of the table's
## variables in the table's order, a comma list, or both.
function [resp, errmsg] = responseNames (lhs, vnames)
  resp = {};
  errmsg = '';
  for part = strtrim (strsplit (lhs, ','))
    ends = strtrim (strsplit (part{1}, '-'));
    idx = cellfun (@(n) find (strcmp (n, vnames), 1), ends, ...
                   'UniformOutput', false);
    miss = cellfun (@isempty, idx);
    if (any (miss) || numel (ends) > 2)
      bad = ends(miss);
      if (isempty (bad))
        bad = part;
      endif
      errmsg = sprintf (strcat ("the model formula names '%s', which is", ...
                                " not a variable of T."), bad{1});
      return;
    endif
    if (numel (ends) == 2)
      if (idx{1} > idx{2})
        errmsg = sprintf ("the response range '%s' runs backwards.", part{1});
        return;
      endif
      resp = [resp, vnames(idx{1}:idx{2})];
    else
      resp = [resp, vnames(idx{1})];
    endif
  endfor
  if (numel (unique (resp)) < numel (resp))
    errmsg = "a response is named more than once.";
  endif
endfunction

## Between-subject terms of the right side, each a cell of variable names,
## and whether the model carries an intercept.
function [terms, intercept, errmsg] = betweenTerms (rhs, vnames, resp)
  terms = {};
  errmsg = '';
  s = parseWilkinsonFormula (['~', rhs]);
  model = s.model;
  empty = cellfun (@isempty, model);
  intercept = any (empty) || isempty (model);
  terms = cellfun (@cellstr, model(! empty), 'UniformOutput', false);
  for i = 1:numel (terms)
    for v = terms{i}
      if (! any (strcmp (v{1}, vnames)))
        errmsg = sprintf (strcat ("the model formula names '%s', which is", ...
                                  " not a variable of T."), v{1});
        return;
      elseif (any (strcmp (v{1}, resp)))
        errmsg = sprintf ("'%s' is a response and a predictor.", v{1});
        return;
      endif
    endfor
  endfor
endfunction

## The between-subject model as MATLAB renders it, '1 + g + x'.
function txt = modelText (rhs, intercept)
  txt = regexprep (strtrim (rhs), '\s*([+-])\s*', ' $1 ');
  if (! strcmp (txt, '1') && intercept && isempty (regexp (txt, '^1\b')))
    txt = ['1 + ', txt];
  endif
endfunction

## Effects-coded design matrix of the terms, the coefficient names, and for
## each term (the intercept first) its name and its columns.  BAD marks the
## subjects missing a predictor or holding a level CODING does not know.
## CODING, one field per variable, is built from T unless given; for a
## categorical variable it holds the level keys and their native values, and
## for a continuous one its mean, set by the caller.
function [X, cnames, tcols, tnames, bad, coding] = effectsDesign (t, terms, intercept, vars, coding)
  n = rows (t);
  bad = false (n, 1);
  given = (nargin > 4);
  if (! given)
    coding = struct ();
  endif
  code = struct ();
  for i = 1:numel (vars)
    v = t.(vars{i});
    if (given)
      cd = coding.(vars{i});
      if (cd.iscat)
        [~, ic, ~, ~, miss] = levelsOf (v, cd.levels);
      else
        [~, ~, ~, ~, miss] = levelsOf (v);
      endif
    else
      [lev, ic, native, iscat, miss] = levelsOf (v);
      cd = struct ('iscat', iscat, 'levels', {lev}, 'native', {native}, ...
                   'mean', NaN);
      coding.(vars{i}) = cd;
    endif
    if (cd.iscat)
      miss |= (ic(:) == 0);
    endif
    bad |= miss(:);
    if (cd.iscat)
      L = numel (cd.levels);
      M = zeros (n, L - 1);
      for l = 1:L-1
        M(:,l) = (ic(:) == l) - (ic(:) == L);
      endfor
      code.(vars{i}) = struct ('M', M, 'names', ...
                               {strcat(vars{i}, '_', cd.levels(1:L-1)(:)')});
    else
      code.(vars{i}) = struct ('M', double (v(:)), 'names', {vars(i)});
    endif
  endfor

  X = zeros (n, 0);
  cnames = {};
  tcols = {};
  tnames = {};
  if (intercept)
    X = ones (n, 1);
    cnames = {'(Intercept)'};
    tcols = {1};
    tnames = {'(Intercept)'};
  endif
  for i = 1:numel (terms)
    M = ones (n, 1);
    nm = {''};
    for v = terms{i}
      c = code.(v{1});
      Mn = zeros (n, 0);
      nn = {};
      for a = 1:columns (M)
        for b = 1:columns (c.M)
          Mn(:,end+1) = M(:,a) .* c.M(:,b);
          if (isempty (nm{a}))
            nn{end+1} = c.names{b};
          else
            nn{end+1} = [nm{a}, ':', c.names{b}];
          endif
        endfor
      endfor
      M = Mn;
      nm = nn;
    endfor
    tcols{end+1} = columns (X) + (1:columns (M));
    tnames{end+1} = strjoin (terms{i}, ':');
    X = [X, M];
    cnames = [cnames, nm];
  endfor
endfunction

## Levels of the variable V: their keys as text, the code of each element
## (0 where it is missing or not among LEV when LEV is given), one native
## value per level, whether V is categorical, and which elements are
## missing.  A numeric V is continuous and has no levels.
function [lev, ic, native, iscat, miss] = levelsOf (v, lev)
  native = {};
  ic = [];
  if (iscategorical (v))
    miss = isundefined (v(:));
    keys = cellstr (v(:));
    iscat = true;
  elseif (iscellstr (v) || isa (v, 'string') || ischar (v) || islogical (v))
    if (islogical (v))
      keys = cellstr (num2str (double (v(:))));
    else
      keys = cellstr (v);
      keys = keys(:);
    endif
    miss = cellfun (@isempty, keys);
    iscat = true;
  else
    v = double (v(:));
    miss = isnan (v);
    iscat = false;
    if (nargin < 2)
      lev = {};
    endif
    return;
  endif
  if (nargin < 2)
    if (iscategorical (v))
      lev = categories (v);
      lev = lev(ismember (lev, keys(! miss)));
    else
      lev = unique (keys(! miss));
    endif
  endif
  [~, ic] = ismember (keys, lev);
  ic(miss) = 0;
  native = cell (numel (lev), 1);
  for l = 1:numel (lev)
    j = find (ic == l, 1);
    if (isempty (j))
      native{l} = lev{l};
    elseif (iscategorical (v) || islogical (v))
      native{l} = v(j);
    else
      native{l} = keys(j);
    endif
  endfor
endfunction

## The within-subject terms of the model WM over the design WD: their
## contrast matrices, W.C, and names, W.names, empty for the constant term,
## the names of the columns of each, W.cols, and the kind of model, W.kind.
function [W, errmsg] = withinTerms (WM, WD, wnames, k)
  W = struct ('C', {{}}, 'names', {{}}, 'cols', {{}}, 'kind', '');
  errmsg = '';
  if (numel (wnames) == 1)
    label = wnames{1};
  else
    label = 'Time';
  endif
  if (isa (WM, 'string') && isscalar (WM))
    WM = char (WM);
  endif
  if (isnumeric (WM))
    if (! (isreal (WM) && ismatrix (WM) && rows (WM) == k && ! isempty (WM)))
      errmsg = sprintf ("'WithinModel' must be a matrix with %d rows.", k);
      return;
    endif
    W.C = {double(WM)};
    W.names = {label};
    W.kind = 'matrix';
  elseif (! (ischar (WM) && isrow (WM)))
    errmsg = "invalid 'WithinModel'.";
  elseif (strcmpi (WM, 'separatemeans'))
    W.C = {successive(k)};
    W.names = {label};
    W.kind = 'separatemeans';
  elseif (strcmpi (WM, 'orthogonalcontrasts'))
    if (numel (wnames) != 1 || ! isnumeric (WD.(wnames{1})))
      errmsg = strcat ("the 'orthogonalcontrasts' model requires a", ...
                       " within-subject design with a single numeric factor.");
      return;
    endif
    v = double (WD.(wnames{1}));
    v = (v(:) - mean (v)) / max (1, std (v));
    [Qp, ~] = qr (v .^ (0:k-1), 0);
    W.C = [{ones(k, 1)}, num2cell(Qp(:,2:k), 1)];
    W.names = [{''}, {wnames{1}}, ...
               arrayfun(@(d) sprintf ('%s^%d', wnames{1}, d), 2:k-1, ...
                        'UniformOutput', false)];
    W.cols = W.names;
    W.cols{1} = 'Constant';
    W.kind = 'orthogonalcontrasts';
  else
    try
      s = parseWilkinsonFormula (['~', WM]);
    catch
      errmsg = "invalid 'WithinModel'.";
      return;
    end_try_catch
    model = s.model;
    empty = cellfun (@isempty, model);
    terms = cellfun (@cellstr, model(! empty), 'UniformOutput', false);
    if (any (! ismember ([terms{:}], wnames)))
      errmsg = "invalid 'WithinModel'.";
      return;
    endif
    [Xw, cn, tc, tn] = effectsDesign (WD, terms, true, unique ([terms{:}], ...
                                      'stable'));
    if (any (empty) || isempty (model))
      W.C = {ones(k, 1)};
      W.names = {''};
      W.cols = {{'(Intercept)'}};
    endif
    for i = 1:numel (terms)
      W.C{end+1} = Xw(:,tc{i+1});
      W.names{end+1} = tn{i+1};
      W.cols{end+1} = cn(tc{i+1});
    endfor
    W.kind = 'formula';
  endif
endfunction

## A cell row as MATLAB displays one, {'a'  'b'}.
function s = cellrow (c)
  if (isempty (c))
    s = sprintf ('{%dx%d cell}', size (c));
  else
    s = ['{', strjoin(strcat ("'", c, "'"), '  '), '}'];
  endif
endfunction

## Expected values are MATLAB R2024a's unless stated otherwise.
%!shared rm, t2, W
%! load fisheriris
%! t = table (species, meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!            'VariableNames', {'species', 'meas1', 'meas2', 'meas3', 'meas4'});
%! Meas = table ([1, 2, 3, 4]', 'VariableNames', {'Measurements'});
%! rm = fitrm (t, 'meas1-meas4 ~ species', 'WithinDesign', Meas);
%! g = categorical ([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]');
%! x = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8]';
%! [j, i] = meshgrid (1:6, 1:12);
%! Y = mod (7 * i .* j + 3 * i + j .^ 2, 13) + 2 * j + 3 * (double (g) == 2) ...
%!     + 0.5 * x;
%! t2 = array2table (Y, 'VariableNames', {'y1', 'y2', 'y3', 'y4', 'y5', 'y6'});
%! t2.g = g;
%! t2.x = x;
%! W = table (categorical ([1, 1, 1, 2, 2, 2]'), ...
%!            categorical ([1, 2, 3, 1, 2, 3]'), 'VariableNames', {'A', 'B'});

## Properties
%!assert_equal (rm.ResponseNames, {'meas1', 'meas2', 'meas3', 'meas4'})
%!assert_equal (rm.BetweenFactorNames, {'species'})
%!assert_equal (rm.BetweenModel, '1 + species')
%!assert_equal (rm.WithinFactorNames, {'Measurements'})
%!assert_equal (rm.WithinModel, 'separatemeans')
%!assert_equal (rm.DFE, 147)
%!assert_equal (size (rm.BetweenDesign), [150, 5])
%!assert_equal (rm.WithinDesign.Properties.RowNames', rm.ResponseNames)
%!assert_equal (rm.Coefficients.Properties.RowNames', ...
%!              {'(Intercept)', 'species_setosa', 'species_versicolor'})
%!assert_equal (table2array (rm.Coefficients), ...
%!              [5.84333333333334, 3.05733333333334, 3.758, 1.19933333333333; ...
%!               -0.837333333333333, 0.370666666666667, -2.296, ...
%!               -0.953333333333333; 0.0926666666666665, ...
%!               -0.287333333333334, 0.502, 0.126666666666667], -1e-12)
%!assert_equal (table2array (rm.Covariance)(:,1), ...
%!              [0.265008163265306; 0.0927210884353741; ...
%!               0.167514285714286; 0.0384013605442177], -1e-12)

## ranova under the default within model
%!test
%! tbl = ranova (rm);
%! assert_equal (tbl.Properties.RowNames', {'(Intercept):Measurements', ...
%!               'species:Measurements', 'Error(Measurements)'});
%!test
%! tbl = ranova (rm);
%! assert_equal (tbl.SumSq, [1656.26325; 282.4665; 35.42275], -1e-12);
%! assert_equal (tbl.DF, [3; 6; 441]);
%!test
%! tbl = ranova (rm);
%! assert_equal (tbl.F(1:2), [6873.28617202222; 586.100394520471], -1e-12);
%! assert_equal (tbl.pValue(2), 1.42713829312908e-206, -1e-10);
%!test
%! tbl = ranova (rm);
%! assert_equal (tbl.pValueGG(1:2), [9.44912100730417e-279; ...
%!                                   4.93131396537195e-156], -1e-10);
%! assert_equal (tbl.pValueHF(1:2), [2.92129935033006e-283; ...
%!                                   1.54056390379818e-158], -1e-10);
%! assert_equal (tbl.pValueLB(1:2), [2.58714529113126e-125; ...
%!                                   9.01514627639756e-71], -1e-10);
%!test
%! [~, A, C, D] = ranova (rm);
%! assert_equal (A, {[1, 0, 0]; [0, 1, 0; 0, 0, 1]});
%! assert_equal (C, [1, 0, 0; -1, 1, 0; 0, -1, 1; 0, 0, -1]);
%! assert_equal (D, 0);
%!test
%! tbl = ranova (rm, 'WithinModel', [1, -1, 0, 0; 0, 1, -1, 0]');
%! assert_equal (tbl.SumSq, [630.067244444444; 282.402488888889; ...
%!                           24.5102666666667], -1e-12);
%! assert_equal (tbl.pValueGG(1:2), [2.73416691642378e-192; ...
%!                                   3.34270901331439e-146], -1e-10);

## orthogonalcontrasts, which R2024a and R2026a fail on inside ranova; the
## expected values are R2026a's anova for the same contrasts
%!test
%! tbl = ranova (rm, 'WithinModel', 'orthogonalcontrasts');
%! assert_equal (tbl.SumSq, [7201.65615000052; 309.6067; 53.87465; ...
%!               1313.01136333336; 35.9780866666667; 13.75805; ...
%!               1.93801666666664; 0.469633333333338; 4.53985; ...
%!               341.313870000001; 246.01878; 17.12485], -1e-11);
%!test
%! tbl = ranova (rm, 'WithinModel', 'orthogonalcontrasts');
%! assert_equal (tbl.F([4, 5, 7, 8, 10, 11]), [14029.0717369107; ...
%!               192.206698623715; 62.7528332433881; 7.60334592552624; ...
%!               2929.84399221016; 1055.91466961754], -1e-11);
%!test
%! tbl = ranova (rm, 'WithinModel', 'orthogonalcontrasts');
%! assert_equal (tbl.Properties.RowNames([4, 7, 10])', ...
%!               {'(Intercept):Measurements', '(Intercept):Measurements^2', ...
%!                '(Intercept):Measurements^3'});

## mauchly and epsilon
%!test
%! tbl = mauchly (rm);
%! assert_equal (table2array (tbl), [0.558144130896497, 84.9761726331633, ...
%!                                   5, 7.61488382466235e-17], -1e-12);
%!test
%! tbl = epsilon (rm);
%! assert_equal (table2array (tbl), [1, 0.751790042430452, ...
%!                                   0.764092317002038, 1/3], -1e-12);
%!test
%! C = [1, -1, 0, 0; 0, 1, -1, 0; 0, 0, 1, -1]';
%! assert_equal (table2array (epsilon (rm, C)), ...
%!               [1, 0.751790042430453, 0.764092317002038, 1/3], -1e-12);
%!test
%! C = [1, -1, 0, 0; 0, 1, -1, 0; 0, 0, 1, -1]';
%! assert_equal (mauchly (rm, C).W, 0.558144130896497, -1e-12);

## A continuous covariate beside a between factor, and a factorial within
## design
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W, 'WithinModel', 'A*B');
%! assert_equal (rm2.BetweenModel, '1 + g + x');
%! assert_equal (rm2.BetweenFactorNames, {'g', 'x'});
%! assert_equal (rm2.DFE, 9);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! assert_equal (table2array (rm2.Coefficients)(:,1), ...
%!               [9.54310344827587; -0.0402298850574699; 0.586206896551724], ...
%!               -1e-12);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W, 'WithinModel', 'A*B');
%! tbl = ranova (rm2);
%! assert_equal (tbl.Properties.RowNames', {'(Intercept):Time', 'g:Time', ...
%!               'x:Time', 'Error(Time)'});
%! assert_equal (tbl.SumSq, [255.837635271832; 105.120215891605; ...
%!                           106.01539408867; 651.790161466886], -1e-12);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = ranova (rm2);
%! assert_equal (tbl.pValueGG(1:3), [0.028891532860066; 0.25048764667133; ...
%!                                   0.247169589665584], -1e-12);
%! assert_equal (tbl.pValueHF(1:3), [0.0115155408306671; 0.23038031489716; ...
%!                                   0.226411224877626], -1e-12);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = ranova (rm2, 'WithinModel', 'A*B');
%! assert_equal (tbl.SumSq, [3499.10306228529; 40.8161242267455; ...
%!               97.0498768472906; 63.3112342638205; 173.443725936885; ...
%!               6.54103670312192; 3.16262999452654; 229.976258894362; ...
%!               60.1151512858879; 32.4673167683903; 78.5480295566504; ...
%!               145.674192665572; 22.2787580490593; 66.1118624200928; ...
%!               24.3047345374931; 276.139709906951], -1e-11);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = ranova (rm2, 'WithinModel', 'A*B');
%! assert_equal (tbl.pValueGG(9:11), [0.0492037168776316; ...
%!               0.167760677557413; 0.0238347461584444], -1e-11);
%! assert_equal (tbl.pValueHF(13:15), [0.488916988695095; ...
%!               0.149284419253376; 0.460679398895663], -1e-11);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! assert_equal (table2array (mauchly (rm2)), [0.109478004315118, ...
%!               15.7054245285396, 14, 0.331692719814282], -1e-11);
%! assert_equal (table2array (epsilon (rm2)), [1, 0.589453986170189, ...
%!               0.907778770724919, 0.2], -1e-12);
%!test
%! rm9 = fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', W, 'WithinModel', 'A+B');
%! tbl = ranova (rm9, 'WithinModel', rm9.WithinModel);
%! assert_equal (tbl.SumSq, [19900.125; 74.0138888888887; 160.361111111111; ...
%!               583.680555555555; 5.01388888888888; 233.138888888889; ...
%!               137.25; 58.5277777777778; 224.222222222222], -1e-12);
%! assert_equal (tbl.pValueGG(7:8), [0.0102255822891222; ...
%!                                   0.103329724488498], -1e-11);

## A numeric within design, Huynh-Feldt capped at 1
%!test
%! rm3 = fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', [1, 2, 4, 8, 16, 32]');
%! assert_equal (rm3.WithinDesign.Time, [1; 2; 4; 8; 16; 32]);
%! assert_equal (rm3.WithinFactorNames, {'Time'});
%!test
%! rm3 = fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', [1, 2, 4, 8, 16, 32]');
%! tbl = ranova (rm3);
%! assert_equal (tbl.SumSq, [722.125; 124.569444444444; 757.805555555556], ...
%!               -1e-12);
%! assert_equal (tbl.pValueGG(1:2), [8.15757856520654e-05; ...
%!                                   0.195494399896943], -1e-11);
%! assert_equal (tbl.pValueHF(1:2), tbl.pValue(1:2));
%!test
%! rm3 = fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', [1, 2, 4, 8, 16, 32]');
%! assert_equal (table2array (epsilon (rm3)), [1, 0.650731538906031, 1, ...
%!                                             0.2], -1e-12);
%! assert_equal (table2array (mauchly (rm3)), [0.172269738544646, ...
%!               14.2454196441778, 14, 0.431591774636887], -1e-11);

## An intercept alone
%!test
%! rm4 = fitrm (t2, 'y1-y6 ~ 1');
%! assert_equal (rm4.BetweenFactorNames, cell (1, 0));
%! assert_equal (rm4.BetweenModel, '1');
%!test
%! rm4 = fitrm (t2, 'y1-y6 ~ 1');
%! tbl = ranova (rm4);
%! assert_equal (tbl.SumSq, [722.125; 882.375], -1e-12);
%! assert_equal (tbl.pValueGG(1), 6.37221717220761e-05, -1e-11);
%! assert_equal (epsilon (rm4).GreenhouseGeisser, 0.693384390469508, -1e-12);

## An interaction between subjects, and a comma list of responses
%!test
%! rm7 = fitrm (t2, 'y1-y6 ~ g*x', 'WithinDesign', W);
%! assert_equal (rm7.BetweenModel, '1 + g*x');
%! assert_equal (rm7.Coefficients.Properties.RowNames', ...
%!               {'(Intercept)', 'g_1', 'x', 'g_1:x'});
%! assert_equal (table2array (rm7.Coefficients)(:,1), [10.8028790057797; ...
%!               -3.00733997232247; 0.37163867256397; 0.659959840447182], ...
%!               -1e-12);
%!test
%! rm8 = fitrm (t2, 'y1,y2,y3 ~ g');
%! assert_equal (rm8.ResponseNames, {'y1', 'y2', 'y3'});
%! assert_equal (rm8.WithinDesign.Time, [1; 2; 3]);

## A subject missing a response is left out of the fit and kept in the design
%!test
%! t3 = t2;
%! t3.y2(3) = NaN;
%! rm5 = fitrm (t3, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! assert_equal (rm5.DFE, 8);
%! assert_equal (rows (rm5.BetweenDesign), 12);
%! assert_equal (table2array (rm5.Coefficients)(:,1), [9.79152291769344; ...
%!               0.248891079349433; 0.594627895515033], -1e-12);
%! assert_equal (ranova (rm5).SumSq, [242.715922479748; 90.7245227703455; ...
%!               106.321888176989; 631.122556267455], -1e-12);
%!test
%! t4 = t2;
%! t4.y2(12) = NaN;
%! t4.y3(12) = NaN;
%! assert_equal (fitrm (t4, 'y1-y6 ~ g', 'WithinDesign', W).DFE, 9);

## A singular covariance fails Mauchly's test outright
%!test
%! t5 = t2;
%! t5.y2 = t5.y1 + 1;
%! tbl = mauchly (fitrm (t5, 'y1-y6 ~ g'));
%! assert_equal (table2array (tbl), [0, Inf, 14, 0]);

## anova
%!test
%! tbl = anova (rm);
%! assert_equal (cellstr (tbl.Between)', {'constant', 'species', 'Error'});
%! assert_equal (cellstr (tbl.Within)', {'Constant', 'Constant', 'Constant'});
%!test
%! tbl = anova (rm);
%! assert_equal (tbl.SumSq, [7201.65615000065; 309.606699999999; ...
%!                           53.87465], -1e-12);
%! assert_equal (tbl.DF, [1; 2; 147]);
%! assert_equal (tbl.F(1:2), [19650.1221641365; 422.389610883782], -1e-12);
%!test
%! assert_equal (anova (rm, 'WithinModel', 'separatemeans').SumSq, ...
%!               anova (rm).SumSq);
%!test
%! tbl = anova (rm, 'WithinModel', [1, -1, 0, 0; 0, 1, -1, 0]');
%! assert_equal (cellstr (tbl.Within([1, 4]))', {'Contrast1', 'Contrast2'});
%! assert_equal (tbl.SumSq, [1164.2694; 114.4624; 28.6582; ...
%!               73.6400666666671; 562.926933333334; 27.943], -1e-12);
%! assert_equal (tbl.F([4, 5]), [387.39898364528; 1480.69747700676], -1e-12);
%!test
%! tbl = anova (rm, 'WithinModel', [1, 1, 1, 1]');
%! assert_equal (tbl.SumSq(1), 28806.6246000026, -1e-12);
## R2026a's values, as for ranova
%!test
%! tbl = anova (rm, 'WithinModel', 'orthogonalcontrasts');
%! assert_equal (cellstr (tbl.Within([1, 4, 7, 10]))', {'Constant', ...
%!               'Measurements', 'Measurements^2', 'Measurements^3'});
%! assert_equal (tbl.SumSq([4, 7, 10]), [1313.01136333336; ...
%!               1.93801666666664; 341.313870000001], -1e-11);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = anova (rm2, 'WithinModel', 'A*B');
%! assert_equal (cellstr (tbl.Within(1:4:end))', {'(Intercept)', 'A_1', ...
%!               'B_1', 'B_2', 'A_1:B_1', 'A_1:B_2'});
%! assert_equal (tbl.SumSq(9:12), [10.539317448175; 0.186084518387707; ...
%!               7.3161945812808; 73.7254720853859], -1e-12);
%! assert_equal (tbl.pValue(22:23), [0.153846787233791; ...
%!               0.29518761419548], -1e-11);

## manova
%!test
%! tbl = manova (rm);
%! assert_equal (cellstr (tbl.Statistic(1:4))', {'Pillai', 'Wilks', ...
%!               'Hotelling', 'Roy'});
%! assert_equal (cellstr (tbl.Between([1, 5]))', {'(Intercept)', 'species'});
%! assert_equal (cellstr (tbl.Within(1)), {'Constant'});
%!test
%! tbl = manova (rm);
%! assert_equal (tbl.Value(5:8), [0.969092455350808; 0.0411531658067895; ...
%!               23.0505039998432; 23.0396981666769], -1e-12);
%! assert_equal (tbl.RSquare(5:8), [0.484546227675404; 0.797137569257416; ...
%!               0.920161286974006; 0.958402139949237], -1e-12);
%!test
%! tbl = manova (rm);
%! assert_equal (tbl.F([5, 6, 8]), [45.7485249172255; 189.92336681764; ...
%!               1121.26531077827], -1e-11);
%! assert_equal ([tbl.df1([5, 6, 8]), tbl.df2([5, 6, 8])], ...
%!               [6, 292; 6, 290; 3, 146]);
%! assert_equal (tbl.pValue([5, 6, 8]), [2.47288601081631e-39; ...
%!               2.39583203496623e-97; 1.47713689662351e-100], -1e-10);
## The Hotelling-Lawley F is R2024a's; its second degrees of freedom are
## McKeon's, and the p-value follows from them
%!test
%! tbl = manova (rm);
%! assert_equal (tbl.F(7), 555.166436056617, -1e-11);
%! n = (147 - 3 - 1) / 2;
%! b = (3 + 2 * n) * (2 + 2 * n) / (2 * (2 * n + 1) * (n - 1));
%! d2 = 4 + 8 / (b - 1);
%! assert_equal ([tbl.df1(7), tbl.df2(7)], [6, d2], -1e-14);
%! assert_equal (tbl.pValue(7), fcdf (tbl.F(7), 6, d2, 'upper'), -1e-12);
%!test
%! tbl = manova (rm, 'By', 'species');
%! assert_equal (cellstr (tbl.Between(1:4:end))', {'species=setosa', ...
%!               'species=versicolor', 'species=virginica'});
%! assert_equal (tbl.Value(1:4:end), [0.982302126154032; ...
%!               0.97000070863222; 0.972606348435736], -1e-12);
%! assert_equal (tbl.F(1), 2682.69152049929, -1e-11);
%!test
%! tbl = manova (rm, 'WithinModel', [1, -1, 0, 0; 0, 1, -1, 0]');
%! assert_equal (cellstr (tbl.Within(1)), {'Specified contrast'});
%! assert_equal (tbl.Value(5), 0.960761007306465, -1e-12);
%! assert_equal (tbl.F(1), 4126.8008429193, -1e-11);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W, 'WithinModel', 'A*B');
%! [tbl, A, C, D] = manova (rm2);
%! assert_equal (rows (tbl), 48);
%! assert_equal (unique (cellstr (tbl.Within), 'stable')', ...
%!               {'(Intercept)', 'A', 'B', 'A:B'});
%! assert_equal (A, {[1, 0, 0]; [0, 1, 0]; [0, 0, 1]});
%! assert_equal (cellfun (@columns, C), [1, 1, 2, 2]);
%! assert_equal (D, 0);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W, 'WithinModel', 'A*B');
%! tbl = manova (rm2);
%! assert_equal (tbl.Value(41:44), [0.486347963843071; 0.513652036156928; ...
%!               0.946843251088535; 0.946843251088535], -1e-12);
%! assert_equal (tbl.pValue(41), 0.069610708832983, -1e-11);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W, 'WithinModel', 'A*B');
%! tbl = manova (rm2, 'WithinModel', 'separatemeans');
%! assert_equal (tbl.Value(1), 0.786551501991551, -1e-12);
%! assert_equal ([tbl.F(3), tbl.df1(3), tbl.df2(3)], ...
%!               [3.68497089148136, 5, 5], -1e-12);
%! assert_equal (tbl.pValue(3), 0.0893171383745427, -1e-11);

## coeftest
%!test
%! tbl = coeftest (rm, [0, 1, 0; 0, 0, 1], [1, -1, 0, 0; 0, 1, -1, 0]');
%! assert_equal (tbl.Value, [0.960761007306465; 0.0420065576612523; ...
%!               22.739922777632; 22.7370251198758], -1e-12);
%! assert_equal (tbl.F([1, 2, 4]), [67.9496579068885; 283.175722000404; ...
%!               1671.17134631087], -1e-11);
%! assert_equal (tbl.F(3), 828.147118575719, -1e-11);
%!test
%! tbl = coeftest (rm, [0, 1, 0; 0, 0, 1], [1, -1, 0, 0; 0, 1, -1, 0]', ...
%!                 [0.5, 0.1; -0.2, 0.3]);
%! assert_equal (tbl.Value, [1.07175176658522; 0.0463873977279332; ...
%!               18.0107848010586; 17.8682529938836], -1e-12);
%! assert_equal (tbl.pValue(4), 1.71285549134517e-94, -1e-10);
%!test
%! tbl = coeftest (rm, [1, 0, 0], [1, -1, 0, 0]', 1);
%! assert_equal (tbl.F, 2454.27144063479 * ones (4, 1), -1e-12);
%! assert_equal ([tbl.df1, tbl.df2], repmat ([1, 147], 4, 1));
%!test
%! tbl = coeftest (rm, [0, 1, 0], [1, -1, 0, 0]');
%! assert_equal (tbl.Value(1), 0.792486767123089, -1e-12);
%! assert_equal (tbl.F(1), 561.38855894648, -1e-12);

## margmean
%!test
%! tbl = margmean (rm, 'species');
%! assert_equal (tbl.species, {'setosa'; 'versicolor'; 'virginica'});
%! assert_equal (tbl.Mean, [2.5355; 3.573; 4.285], -1e-13);
%! assert_equal (tbl.StdErr, 0.0428073718935813 * ones (3, 1), -1e-12);
%! assert_equal (tbl.Lower(1), 2.45090264579764, -1e-12);
%!test
%! tbl = margmean (rm, 'Measurements');
%! assert_equal (tbl.Measurements, [1; 2; 3; 4]);
%! assert_equal (tbl.StdErr, [0.0420323814271256; 0.0277353871557668; ...
%!               0.0351366622491893; 0.0167096045540803], -1e-12);
%!test
%! tbl = margmean (rm, {'species', 'Measurements'});
%! assert_equal (rows (tbl), 12);
%! assert_equal (tbl.Mean(1:4), [5.006; 3.428; 1.462; 0.246], -1e-12);
%! assert_equal (tbl.StdErr(5), 0.072802220194896, -1e-12);
%!test
%! tbl = margmean (rm, 'species', 'Alpha', 0.01);
%! assert_equal (tbl.Lower(1), 2.42378611949266, -1e-12);
## A continuous predictor enters at its mean
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = margmean (rm2, 'g');
%! assert_equal (tbl.Mean, [15.8555692391899; 17.3944307608101], -1e-12);
%! assert_equal (tbl.StdErr, 0.446919101103406 * ones (2, 1), -1e-12);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = margmean (rm2, {'g', 'A'});
%! assert_equal (tbl.Mean, [13.3163656267105; 18.3947728516694; ...
%!               14.2391899288451; 20.549671592775], -1e-12);
%! assert_equal (tbl.StdErr(1:2), [1.12245927784633; 0.768527281387124], ...
%!               -1e-12);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = margmean (rm2, {'B', 'g'});
%! assert_equal (tbl.Mean([1, 4, 6]), [13.8103448275862; 16.301724137931; ...
%!               19.6919129720854], -1e-12);
%! assert_equal (tbl.StdErr([1, 3, 5]), [0.863939313862951; ...
%!               0.870431298727607; 0.688406615362842], -1e-12);

## grpstats
%!test
%! tbl = grpstats (rm, 'species');
%! assert_equal (tbl.GroupCount, [200; 200; 200]);
%! assert_equal (tbl.mean, [2.5355; 3.573; 4.285], -1e-13);
%! assert_equal (tbl.std, [1.84834293572383; 1.76238503313695; ...
%!               1.91538993235446], -1e-12);
%!test
%! tbl = grpstats (rm, 'Measurements');
%! assert_equal (tbl.GroupCount, 150 * ones (4, 1));
%! assert_equal (tbl.std, [0.828066127977863; 0.435866284936698; ...
%!               1.76529823325947; 0.762237668960347], -1e-12);
%!test
%! tbl = grpstats (rm, 'species', {'min', 'max'});
%! assert_equal (tbl.Properties.VariableNames, {'species', 'GroupCount', ...
%!               'min', 'max'});
%! assert_equal ([tbl.min, tbl.max], [0.1, 5.8; 1, 7; 1.4, 7.9]);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = grpstats (rm2, {'g', 'B'});
%! assert_equal (tbl.GroupCount, 12 * ones (6, 1));
%! assert_equal (tbl.std([1, 5]), [4.56767297030297; 3.29944898981082], ...
%!               -1e-12);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = grpstats (rm2, 'g', {'mean', 'sem', 'var', 'range', 'numel', ...
%!                            'median'});
%! assert_equal ([tbl.sem, tbl.var], [0.797344578852233, ...
%!               22.8873015873016; 0.874599933736397, 27.5373015873016], ...
%!               -1e-12);
%! assert_equal ([tbl.range, tbl.numel, tbl.median], [19, 36, 16.25; ...
%!                                                    22, 36, 17.5]);

## multcompare
%!test
%! tbl = multcompare (rm, 'species', 'ComparisonType', 'bonferroni');
%! assert_equal (tbl.Difference(1:2), [-1.0375; -1.7495], -1e-12);
%! assert_equal (tbl.StdErr(1), 0.0605387659014515, -1e-12);
%! assert_equal (tbl.pValue(1:2), [2.15970302648855e-36; ...
%!               5.03902053290415e-62], -1e-9);
%! assert_equal (tbl.Lower(1), -1.18410589638625, -1e-12);
%!test
%! tbl = multcompare (rm, 'species', 'ComparisonType', 'lsd');
%! assert_equal (tbl.pValue(1), 7.19901008829516e-37, -1e-9);
%! assert_equal (tbl.Lower(1), -1.15713872565386, -1e-12);
%!test
%! tbl = multcompare (rm, 'species', 'ComparisonType', 'scheffe');
%! assert_equal (tbl.pValue(1), 8.97549175066402e-36, -1e-9);
%! assert_equal (tbl.Lower(1), -1.18720639857468, -1e-12);
## R2024a floors these Tukey-Kramer p-values, and takes its limits from the
## studentized range on infinite degrees of freedom; the quantile here is
## R 4.5.0's qtukey (0.95, 3, 147)
%!test
%! tbl = multcompare (rm, 'species');
%! assert_equal (tbl.Lower(1), -1.0375 - 3.34842406186643 / sqrt (2) ...
%!               * 0.0605387659014515, -1e-10);
%! assert_equal (all (tbl.pValue < 1e-20), true);
%!test
%! tbl = multcompare (rm, 'Measurements', 'By', 'species');
%! assert_equal (rows (tbl), 36);
%! assert_equal (tbl.Properties.VariableNames(1:3), {'species', ...
%!               'Measurements_1', 'Measurements_2'});
%! assert_equal (tbl.StdErr(1:3), [0.0624425722558894; 0.047993196796791; ...
%!               0.0678361370996215], -1e-12);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = multcompare (rm2, 'B');
%! assert_equal (tbl.StdErr([1, 2, 4]), [0.914225240599476; ...
%!               0.826222282469959; 0.710493930524169], -1e-12);
%! ## R2024a's studentized range is good to about 2e-7 here
%! assert_equal (tbl.pValue([1, 2, 4]), [0.278857228960382; ...
%!               0.00693824145665634; 0.0634575438808659], -1e-6);
%! assert_equal (tbl.Lower(1:2), [-4.05252199083511; -5.68181723897677], ...
%!               -1e-7);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = multcompare (rm2, 'B', 'By', 'g');
%! assert_equal (tbl.Difference([1, 8]), [-2.88793103448276; ...
%!               -3.50225779967159], -1e-12);
%! assert_equal (tbl.pValue([1, 10]), [0.122845632074643; ...
%!               0.0214216738123368], -1e-7);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = multcompare (rm2, 'g', 'By', 'A');
%! assert_equal (tbl.StdErr([1, 3]), [1.60451755027924; 1.09858373946539], ...
%!               -1e-12);
%! assert_equal (tbl.pValue([1, 3]), [0.579288143050048; ...
%!               0.0814446080647632], -1e-7);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! tbl = multcompare (rm2, 'B', 'ComparisonType', 'dunn-sidak');
%! assert_equal (tbl.pValue([1, 2, 4]), [0.353395443395346; ...
%!               0.00819161456139628; 0.0787120565287793], -1e-10);
%! assert_equal (tbl.Lower(1), -4.17216326867877, -1e-10);

## predict and random
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! assert_equal (predict (rm2, t2(1:2,:))(1,:), [11.2614942528736, ...
%!               13.2783251231527, 13.1005747126437, 14.7040229885058, ...
%!               20.3940886699507, 18.4835796387521], -1e-12);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! [~, yci] = predict (rm2, t2(1:2,:));
%! assert_equal (size (yci), [2, 6, 2]);
%! assert_equal ([yci(1,1,1), yci(1,1,2)], [7.64052800884495, ...
%!                                          14.8824604969022], -1e-12);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! t4 = t2(1:2,:);
%! t4.x(1) = NaN;
%! yp = predict (rm2, t4);
%! assert_equal (isnan (yp(:,1)), [true; false]);
%!test
%! rm2 = fitrm (t2, 'y1-y6 ~ g + x', 'WithinDesign', W);
%! assert_equal (size (predict (rm2)), [12, 6]);
%!test
%! rm3 = fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', [1, 2, 4, 8, 16, 32]', ...
%!              'WithinModel', 'orthogonalcontrasts');
%! assert_equal (predict (rm3, t2(1,:), 'WithinDesign', [3, 5]'), ...
%!               [13.8420137222029, 14.4827614967302], -1e-11);
%!test
%! rm9 = fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', W, 'WithinModel', 'A+B');
%! assert_equal (predict (rm9, t2(1,:)), [10.9166666666667, ...
%!               14.1666666666667, 14, 16.0833333333333, 19.3333333333333, ...
%!               19.1666666666667], -1e-12);
%!warning<RepeatedMeasuresModel.predict: the 'separatemeans' model does not use the 'WithinDesign' given.> ...
%! predict (rm, rm.BetweenDesign(1,:), 'WithinDesign', [1.5, 2.5]');
%!test
%! assert_equal (size (random (rm)), [150, 4]);
%! assert_equal (size (random (rm, rm.BetweenDesign([1, 51, 101],:))), [3, 4]);

## Test input validation
%!error<RepeatedMeasuresModel: too few input arguments.> RepeatedMeasuresModel (1)
%!error<RepeatedMeasuresModel: T must be a table.> fitrm ([1, 2, 3], 'y1-y6 ~ g')
%!error<RepeatedMeasuresModel: MODELSPEC must be a character vector or a string scalar.> ...
%! fitrm (t2, 5)
%!error<RepeatedMeasuresModel: MODELSPEC must be of the form 'responses ~ terms'.> ...
%! fitrm (t2, 'y1-y6')
%!error<RepeatedMeasuresModel: the model formula names 'q', which is not a variable of T.> ...
%! fitrm (t2, 'y1-y6 ~ q')
%!error<RepeatedMeasuresModel: the model formula names 'y9', which is not a variable of T.> ...
%! fitrm (t2, 'y1-y9 ~ g')
%!error<RepeatedMeasuresModel: the response range 'y6-y1' runs backwards.> ...
%! fitrm (t2, 'y6-y1 ~ g')
%!error<RepeatedMeasuresModel: a response is named more than once.> ...
%! fitrm (t2, 'y1,y1 ~ g')
%!error<RepeatedMeasuresModel: 'y1' is a response and a predictor.> ...
%! fitrm (t2, 'y1-y6 ~ y1')
%!error<RepeatedMeasuresModel: response 'g' must be numeric.> ...
%! fitrm (t2, 'g,y1 ~ x')
%!error<RepeatedMeasuresModel: invalid optional paired argument.> ...
%! fitrm (t2, 'y1-y6 ~ g', 'Nonsense', 1)
%!error<RepeatedMeasuresModel: 'WithinDesign' has 3 points, where 6 are required.> ...
%! fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', [1, 2, 3]')
%!error<RepeatedMeasuresModel: 'WithinDesign' must be a table or a numeric vector.> ...
%! fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', {1})
%!error<RepeatedMeasuresModel: invalid 'WithinModel'.> ...
%! fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', W, 'WithinModel', 'C*D')
%!error<RepeatedMeasuresModel: the 'orthogonalcontrasts' model requires a within-subject design with a single numeric factor.> ...
%! fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', W, 'WithinModel', 'orthogonalcontrasts')
%!error<RepeatedMeasuresModel: the between-subject design is rank deficient or leaves no error degrees of freedom.> ...
%! fitrm ([t2, table(2 * t2.x, 'VariableNames', {'x2'})], 'y1-y6 ~ x + x2')
%!error<RepeatedMeasuresModel.ranova: 'WithinModel' must be a matrix with 4 rows.> ...
%! ranova (rm, 'WithinModel', ones (5, 1))
%!error<RepeatedMeasuresModel.ranova: invalid 'WithinModel'.> ...
%! ranova (rm, 'WithinModel', 'C*D')
%!error<RepeatedMeasuresModel.ranova: invalid optional paired argument.> ...
%! ranova (rm, 'Nonsense', 1)
%!error<RepeatedMeasuresModel.mauchly: C must be a matrix with 4 rows.> ...
%! mauchly (rm, ones (5, 1))
%!error<RepeatedMeasuresModel.epsilon: C must be a matrix with 4 rows.> ...
%! epsilon (rm, ones (5, 1))
%!error<RepeatedMeasuresModel.anova: invalid 'WithinModel'.> ...
%! anova (rm, 'WithinModel', 'C*D')
%!error<RepeatedMeasuresModel.anova: invalid optional paired argument.> ...
%! anova (rm, 'Nonsense', 1)
%!error<RepeatedMeasuresModel.manova: the 'orthogonalcontrasts' model cannot be used with manova.> ...
%! manova (rm, 'WithinModel', 'orthogonalcontrasts')
%!error<RepeatedMeasuresModel.manova: 'By' must name a between-subject factor, not 'g'.> ...
%! manova (rm, 'By', 'g')
%!error<RepeatedMeasuresModel.manova: 'By' must name a categorical factor, not 'x'.> ...
%! manova (fitrm (t2, 'y1-y6 ~ x'), 'By', 'x')
%!error<RepeatedMeasuresModel.manova: invalid optional paired argument.> ...
%! manova (rm, 'Nonsense', 1)
%!error<RepeatedMeasuresModel.coeftest: too few input arguments.> ...
%! coeftest (rm, [0, 1, 0])
%!error<RepeatedMeasuresModel.coeftest: A must be a matrix with 3 columns.> ...
%! coeftest (rm, [1, 0], [1, -1, 0, 0]')
%!error<RepeatedMeasuresModel.coeftest: C must be a matrix with 4 rows.> ...
%! coeftest (rm, [0, 1, 0], [1, -1, 0]')
%!error<RepeatedMeasuresModel.coeftest: D must be a scalar or a 1-by-1 matrix.> ...
%! coeftest (rm, [0, 1, 0], [1, -1, 0, 0]', [1, 2])
%!error<RepeatedMeasuresModel.margmean: too few input arguments.> margmean (rm)
%!error<RepeatedMeasuresModel.margmean: 'nosuch' is not a categorical factor of the model.> ...
%! margmean (rm, 'nosuch')
%!error<RepeatedMeasuresModel.margmean: 'x' is not a categorical factor of the model.> ...
%! margmean (fitrm (t2, 'y1-y6 ~ g + x'), 'x')
%!error<RepeatedMeasuresModel.margmean: 'Alpha' must be a scalar in \(0, 1\).> ...
%! margmean (rm, 'species', 'Alpha', 1)
%!error<RepeatedMeasuresModel.margmean: invalid optional paired argument.> ...
%! margmean (rm, 'species', 'Nonsense', 1)
%!error<RepeatedMeasuresModel.grpstats: too few input arguments.> grpstats (rm)
%!error<RepeatedMeasuresModel.grpstats: 'nosuch' is not a categorical factor of the model.> ...
%! grpstats (rm, 'nosuch')
%!error<RepeatedMeasuresModel.grpstats: STATS must name statistics among 'mean', 'median', 'std', 'var', 'sem', 'min', 'max', 'range' and 'numel'.> ...
%! grpstats (rm, 'species', 'mode')
%!error<RepeatedMeasuresModel.multcompare: too few input arguments.> multcompare (rm)
%!error<RepeatedMeasuresModel.multcompare: 'nosuch' is not a categorical factor of the model.> ...
%! multcompare (rm, 'nosuch')
%!error<RepeatedMeasuresModel.multcompare: 'ComparisonType' must be 'tukey-kramer', 'bonferroni', 'dunn-sidak', 'lsd' or 'scheffe'.> ...
%! multcompare (rm, 'species', 'ComparisonType', 'hsd')
%!error<RepeatedMeasuresModel.multcompare: 'By' must differ from the factor compared.> ...
%! multcompare (rm, 'species', 'By', 'species')
%!error<RepeatedMeasuresModel.multcompare: 'Alpha' must be a scalar in \(0, 1\).> ...
%! multcompare (rm, 'species', 'Alpha', 0)
%!error<RepeatedMeasuresModel.multcompare: invalid optional paired argument.> ...
%! multcompare (rm, 'species', 'Nonsense', 1)
%!error<RepeatedMeasuresModel.predict: TNEW has no variable 'species'.> ...
%! predict (rm, table ([1; 2], 'VariableNames', {'z'}))
%!error<RepeatedMeasuresModel.predict: 'Alpha' must be a scalar in \(0, 1\).> ...
%! predict (rm, 'Alpha', 2)
%!error<RepeatedMeasuresModel.predict: invalid optional paired argument.> ...
%! predict (rm, 'Nonsense', 1)
%!error<RepeatedMeasuresModel.predict: 'WithinDesign' must hold the within-subject factors.> ...
%! predict (fitrm (t2, 'y1-y6 ~ g', 'WithinDesign', W, 'WithinModel', 'A+B'), ...
%!          'WithinDesign', {1})
%!error<RepeatedMeasuresModel.random: TNEW must be a table.> random (rm, 1)
%!error<RepeatedMeasuresModel.random: TNEW has no variable 'species'.> ...
%! random (rm, table ([1; 2], 'VariableNames', {'z'}))
