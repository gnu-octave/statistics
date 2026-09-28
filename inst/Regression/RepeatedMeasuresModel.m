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
    ## matrix, and for each between term its name and its columns.
    X_          = [];
    R_          = [];
    B_          = [];
    TermNames_  = {};
    TermCols_   = {};
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
      [X, cnames, tcols, tnames, bad] = effectsDesign (t, terms, intercept, vars);
      ok = ! bad & all (! isnan (Y), 2);
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
      if (! (W > 0))
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
## subjects missing a predictor.
function [X, cnames, tcols, tnames, bad] = effectsDesign (t, terms, intercept, vars)
  n = rows (t);
  bad = false (n, 1);
  code = struct ();
  for i = 1:numel (vars)
    v = t.(vars{i});
    if (iscategorical (v))
      miss = isundefined (v);
      lev = categories (v(! miss));
      lev = lev(ismember (lev, cellstr (v(! miss))));
      [~, ic] = ismember (cellstr (v), lev);
      iscat = true;
    elseif (iscellstr (v) || isa (v, 'string') || ischar (v) || islogical (v))
      if (ischar (v))
        v = cellstr (v);
      elseif (islogical (v))
        v = cellstr (num2str (double (v(:))));
      else
        v = cellstr (v);
      endif
      miss = cellfun (@isempty, v);
      lev = unique (v(! miss));
      [~, ic] = ismember (v, lev);
      iscat = true;
    else
      v = double (v(:));
      miss = isnan (v);
      iscat = false;
    endif
    bad |= miss(:);
    if (iscat)
      L = numel (lev);
      M = zeros (n, L - 1);
      for l = 1:L-1
        M(:,l) = (ic(:) == l) - (ic(:) == L);
      endfor
      code.(vars{i}) = struct ('M', M, 'names', ...
                               {strcat(vars{i}, '_', lev(1:L-1)(:)')});
    else
      code.(vars{i}) = struct ('M', v, 'names', {vars(i)});
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

## The within-subject terms of the model WM over the design WD: their
## contrast matrices, W.C, and names, W.names, empty for the constant term.
function [W, errmsg] = withinTerms (WM, WD, wnames, k)
  W = struct ('C', {{}}, 'names', {{}});
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
  elseif (! (ischar (WM) && isrow (WM)))
    errmsg = "invalid 'WithinModel'.";
  elseif (strcmpi (WM, 'separatemeans'))
    W.C = {successive(k)};
    W.names = {label};
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
    [Xw, ~, tc, tn] = effectsDesign (WD, terms, true, unique ([terms{:}], ...
                                     'stable'));
    if (any (empty) || isempty (model))
      W.C = {ones(k, 1)};
      W.names = {''};
    endif
    for i = 1:numel (terms)
      W.C{end+1} = Xw(:,tc{i+1});
      W.names{end+1} = tn{i+1};
    endfor
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
