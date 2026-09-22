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

## -*- texinfo -*-
## @deftp {statistics} lime
##
## Local interpretable model-agnostic explanations for a fitted model.
##
## A @code{lime} object explains one prediction by fitting a simple model,
## a linear one or a shallow decision tree, over observations drawn around
## the query point and weighted by how near they lie to it.  The simple
## model is readable where the fitted one is not, and it is accurate near
## the query point rather than everywhere.
##
## @code{@var{explainer} = lime (@var{Mdl})} builds an explainer for the
## fitted model @var{Mdl} over the observations it was fitted on.  A compact
## model keeps none, so it must be given them as @var{X}, and so must a
## function handle.  The observations to fit the simple model on are drawn
## at once, and nothing is explained until a query point is given, either to
## the constructor as @qcode{'QueryPoint'} or afterwards to @code{fit}.
##
## @code{@var{explainer} = lime (@var{Mdl}, @var{X})} takes the observations
## the draw is fitted to as @var{X}, a real numeric matrix of one column per
## predictor.
##
## @code{@var{explainer} = lime (@var{fun}, @var{X})} takes a function
## handle in place of a model.  @var{fun} is called with a matrix of
## observations and answers with one column holding one value for each, and
## @qcode{'Type'} must say whether those values are a response or a label.
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{'Type'} @tab @tab Whether the model answers with a response,
## @qcode{'regression'}, or with a label, @qcode{'classification'}.  It is
## taken from a fitted model and is required for a function handle.
##
## @item @qcode{'DataLocality'} @tab @tab Where the observations are drawn
## from: @qcode{'global'}, the default, fits the distribution to the whole
## of @var{X}; @qcode{'local'} fits it to the @qcode{'NumNeighbors'}
## observations nearest the query point, which must then be known.
##
## @item @qcode{'NumNeighbors'} @tab @tab How many neighbours
## @qcode{'local'} fits to, 1500 by default.
##
## @item @qcode{'NumSyntheticData'} @tab @tab How many observations to draw,
## 5000 by default.
##
## @item @qcode{'CustomSyntheticData'} @tab @tab Observations to use instead
## of drawing any, one row each.  Nothing is drawn where it is given.
##
## @item @qcode{'CategoricalPredictors'} @tab @tab The predictors whose
## values are levels, taken as by every learner of this package.  It applies
## only to a function handle, a model being asked for its own.
##
## @item @qcode{'QueryPoint'} @tab @tab The observation to explain.  Given
## with @qcode{'NumImportantPredictors'} it is explained at once, otherwise
## it is left to @code{fit}.
##
## @item @qcode{'NumImportantPredictors'} @tab @tab How many predictors the
## simple model is fitted on.
## @end multitable
##
## Every option @code{fit} takes may also be given here, in which case it
## stands as the default for every later @code{fit}.
##
## @code{'UseParallel'} is not implemented and is refused rather than
## ignored.
##
## @seealso{shapley, partialDependence, PredictiveModel}
## @end deftp

classdef lime

  properties (SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {lime} {property} BlackboxModel
    ##
    ## The model being explained
    ##
    ## The fitted model, or the function handle, the explainer was built on.
    ## This property is read-only.
    ##
    ## @end deftp
    BlackboxModel = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} DataLocality
    ##
    ## Where the observations were drawn from
    ##
    ## @qcode{'global'} where the distribution was fitted to the whole of
    ## @qcode{X}, @qcode{'local'} where it was fitted to the neighbours of
    ## the query point.  This property is read-only.
    ##
    ## @end deftp
    DataLocality = 'global';

    ## -*- texinfo -*-
    ## @deftp {lime} {property} CategoricalPredictors
    ##
    ## The categorical predictors
    ##
    ## The indices of the predictors whose values are levels, empty where
    ## there are none.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} Type
    ##
    ## What the model answers with
    ##
    ## @qcode{'regression'} or @qcode{'classification'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Type = '';

    ## -*- texinfo -*-
    ## @deftp {lime} {property} X
    ##
    ## The observations the draw was fitted to
    ##
    ## A real numeric matrix of one row per observation and one column per
    ## predictor.  This property is read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} QueryPoint
    ##
    ## The observation explained
    ##
    ## One row holding one value per predictor, empty until a query point is
    ## given.  This property is read-only.
    ##
    ## @end deftp
    QueryPoint = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} NumImportantPredictors
    ##
    ## How many predictors the simple model was asked for
    ##
    ## Empty until a query point is given.  Fewer may be used, where a
    ## predictor adds nothing.  This property is read-only.
    ##
    ## @end deftp
    NumImportantPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} NumSyntheticData
    ##
    ## How many observations were drawn
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumSyntheticData = 5000;

    ## -*- texinfo -*-
    ## @deftp {lime} {property} SyntheticData
    ##
    ## The observations the simple model is fitted over
    ##
    ## One row each, drawn around the query point or given outright as
    ## @qcode{'CustomSyntheticData'}.  This property is read-only.
    ##
    ## @end deftp
    SyntheticData = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} Fitted
    ##
    ## What the model answers over the drawn observations
    ##
    ## One response per observation for a regression model, or one label for
    ## a classifier, keeping the type of the response the model was fitted
    ## with.  This property is read-only.
    ##
    ## @end deftp
    Fitted = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} SimpleModel
    ##
    ## The simple model fitted around the query point
    ##
    ## A @code{RegressionLinear}, @code{ClassificationLinear},
    ## @code{RegressionTree} or @code{ClassificationTree}, fitted on the
    ## important predictors alone and weighted by nearness to the query
    ## point.  For a classifier it answers @math{1} for the class the
    ## explained model predicted and @math{-1} for any other, which is what
    ## makes it a single model however many classes there are.  It is empty
    ## until a query point is given.  This property is read-only.
    ##
    ## @end deftp
    SimpleModel = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} ImportantPredictors
    ##
    ## The predictors the simple model was fitted on
    ##
    ## Their indices, in increasing order, empty until a query point is
    ## given.  This property is read-only.
    ##
    ## @end deftp
    ImportantPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} BlackboxFitted
    ##
    ## What the explained model answers at the query point
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    BlackboxFitted = [];

    ## -*- texinfo -*-
    ## @deftp {lime} {property} SimpleModelFitted
    ##
    ## What the simple model answers at the query point
    ##
    ## Where the two agree the simple model is worth reading; where they do
    ## not, it explains nothing.  This property is read-only.
    ##
    ## @end deftp
    SimpleModelFitted = [];

  endproperties

  properties (Access = protected, Hidden)
    PredictorNames_ = {};   # one name per column of X
    ClassNames_ = [];       # the classes, in the model's own order
    Opts_ = [];             # the fit options given at construction
    NumNeighbors_ = 1500;   # how many neighbours 'local' fits to
    Custom_ = false;        # whether the observations were given outright
  endproperties

  methods (Hidden)

    function disp (this)

      printf ('\n  lime with properties:\n\n');
      printf ('%26s: %s\n', 'Type', this.Type);
      printf ('%26s: %s\n', 'DataLocality', this.DataLocality);
      printf ('%26s: %d\n', 'NumSyntheticData', this.NumSyntheticData);
      if (isempty (this.SimpleModel))
        printf ('%26s: %s\n', 'SimpleModel', '[]');
      else
        printf ('%26s: %s\n', 'SimpleModel', class (this.SimpleModel));
        printf ('%26s: %s\n', 'ImportantPredictors', ...
                mat2str (this.ImportantPredictors));
      endif
      printf ('\n');

    endfunction

    function display (this)

      inName = inputname (1);
      if (! isempty (inName))
        printf ('%s =\n', inName);
      endif
      disp (this);

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {lime} {@var{obj} =} lime (@var{Mdl})
    ## @deftypefnx {lime} {@var{obj} =} lime (@var{Mdl}, @var{X})
    ## @deftypefnx {lime} {@var{obj} =} lime (@var{fun}, @var{X})
    ## @deftypefnx {lime} {@var{obj} =} lime (@dots{}, @var{name}, @var{value})
    ##
    ## Build a local interpretable model-agnostic explainer.
    ##
    ## The arguments are those described for the class.  The observations
    ## the simple model is fitted over are drawn here; where
    ## @qcode{'QueryPoint'} and @qcode{'NumImportantPredictors'} are both
    ## given the explanation is computed at once, otherwise it is left to
    ## @code{fit}.
    ##
    ## @seealso{lime, lime.fit, lime.plot}
    ## @end deftypefn

    function this = lime (blackbox, varargin)

      ## Input validation
      if (nargin < 1)
        error ("lime: too few input arguments.");
      endif
      if (! (isa (blackbox, 'PredictiveModel')
             || is_function_handle (blackbox)))
        error (strcat ("lime: BLACKBOX must be a fitted model that", ...
                       " predicts, or a function handle."));
      endif

      ## The observations come before the options, as MATLAB orders them
      args = varargin;
      Data = [];
      if (! isempty (args) && ! (ischar (args{1}) || isa (args{1}, 'string')))
        Data = args{1};
        args(1) = [];
      endif

      optNames = {'Type', 'DataLocality', 'NumNeighbors', ...
                  'NumSyntheticData', 'CustomSyntheticData', ...
                  'CategoricalPredictors', 'QueryPoint', ...
                  'NumImportantPredictors', 'Distance', 'KernelWidth', ...
                  'SimpleModelType', 'BetaTolerance', 'Cov', 'P', ...
                  'Scale', 'UseParallel'};
      dfValues = repmat ({[]}, 1, numel (optNames));
      [Type, Loc, NumNb, NumSyn, Custom, CatPred, QP, NIP, Dist, KW, ...
       SMT, BT, Cov, PP, Scale, Par, rem] = ...
                          parsePairedArguments (optNames, dfValues, args);
      if (! isempty (rem))
        error ("lime: unknown optional argument or misplaced value.");
      endif
      if (! isempty (Par))
        error ("lime: 'UseParallel' is not implemented.");
      endif

      [F, errmsg] = limeFrame (blackbox, Data, CatPred, Type);
      if (! isempty (errmsg))
        error ("lime: %s", errmsg);
      endif

      this.BlackboxModel = blackbox;
      this.X = F.X;
      this.Type = F.Type;
      this.CategoricalPredictors = F.Cat;
      this.PredictorNames_ = F.PredictorNames;
      this.ClassNames_ = F.ClassNames;

      M = columns (F.X);
      this.DataLocality = limeCheckWord (Loc, {'global', 'local'}, ...
                                         'lime', 'DataLocality', 'global');
      this.NumNeighbors_ = limeCheckCount (NumNb, 'lime', 'NumNeighbors', ...
                                           1500);
      this.Opts_ = limeFitOpts (Dist, KW, SMT, BT, Cov, PP, Scale, ...
                                'lime', F);

      if (! isempty (QP))
        this.QueryPoint = limeCheckPoint (QP, M, 'lime');
      endif

      ## The observations are drawn here, not at the explanation, so that
      ## an explainer carries them before any query point is given
      if (! isempty (Custom))
        if (! (isnumeric (Custom) && isreal (Custom) && ismatrix (Custom)
               && ndims (Custom) == 2 && ! isempty (Custom)))
          error (strcat ("lime: 'CustomSyntheticData' must be a real", ...
                         " numeric matrix."));
        endif
        if (columns (Custom) != M)
          error (strcat ("lime: 'CustomSyntheticData' must have one", ...
                         " column per predictor of the model."));
        endif
        this.SyntheticData = double (Custom);
        this.NumSyntheticData = rows (Custom);
        this.Custom_ = true;
      else
        this.NumSyntheticData = limeCheckCount (NumSyn, 'lime', ...
                                                'NumSyntheticData', 5000);
        this.SyntheticData = limeSynth (this);
      endif
      this.Fitted = limePredict (blackbox, this.SyntheticData, F.Type);

      ## A query point with a count of predictors is explained at once
      if (! isempty (QP) && ! isempty (NIP))
        this = fit (this, this.QueryPoint, NIP);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {lime} {@var{obj} =} fit (@var{obj}, @var{queryPoint}, @var{numImportantPredictors})
    ## @deftypefnx {lime} {@var{obj} =} fit (@dots{}, @var{name}, @var{value})
    ##
    ## Fit the simple model around one query point.
    ##
    ## @var{queryPoint} is one row holding one value per predictor and
    ## @var{numImportantPredictors} how many predictors the simple model is
    ## fitted on.  Fewer are used where a predictor adds nothing, and
    ## @code{ImportantPredictors} says which were.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'Distance'} @tab @tab How nearness to the query point is
    ## measured.  Where the predictors hold numbers it is one of
    ## @qcode{'euclidean'}, the default, @qcode{'squaredeuclidean'},
    ## @qcode{'seuclidean'}, @qcode{'mahalanobis'}, @qcode{'cityblock'},
    ## @qcode{'minkowski'}, @qcode{'chebychev'}, @qcode{'cosine'},
    ## @qcode{'correlation'} and @qcode{'spearman'}; where they hold levels
    ## it is @qcode{'goodall3'}, the default, or @qcode{'ofd'}.  A function
    ## handle is also taken, which MATLAB does not; see below.
    ##
    ## @item @qcode{'KernelWidth'} @tab @tab How sharply the weight falls
    ## away with distance, from just above 0 to 1.  The default is 0.75.
    ##
    ## @item @qcode{'SimpleModelType'} @tab @tab @qcode{'linear'}, the
    ## default, or @qcode{'tree'}.
    ##
    ## @item @qcode{'BetaTolerance'} @tab @tab Relative tolerance the linear
    ## simple model is fitted to, @math{10^{-4}} by default.
    ##
    ## @item @qcode{'Cov'}, @qcode{'P'}, @qcode{'Scale'} @tab @tab Passed to
    ## the distance that takes them, as @code{pdist2} takes them.
    ## @end multitable
    ##
    ## The weight of a drawn observation is
    ## @math{exp (-0.5 (d / max (d) / w)^2)}, with @math{d} its distance
    ## from the query point and @math{w} the kernel width.
    ##
    ## @seealso{lime, lime.plot}
    ## @end deftypefn

    function this = fit (this, QueryPoint, NumImportant, varargin)

      ## Input validation
      if (nargin < 3)
        error ("lime.fit: too few input arguments.");
      endif
      M = numel (this.PredictorNames_);
      QueryPoint = limeCheckPoint (QueryPoint, M, 'lime.fit');
      if (! (isnumeric (NumImportant) && isscalar (NumImportant)
             && isreal (NumImportant) && isfinite (NumImportant)
             && NumImportant == fix (NumImportant) && NumImportant > 0))
        error (strcat ("lime.fit: NUMIMPORTANTPREDICTORS must be a", ...
                       " positive integer."));
      endif
      if (NumImportant > M)
        error (strcat ("lime.fit: NUMIMPORTANTPREDICTORS must not exceed", ...
                       " the %d predictors of the model."), M);
      endif

      optNames = {'Distance', 'KernelWidth', 'SimpleModelType', ...
                  'BetaTolerance', 'Cov', 'P', 'Scale'};
      dfValues = repmat ({[]}, 1, numel (optNames));
      [Dist, KW, SMT, BT, Cov, PP, Scale, rem] = ...
                    parsePairedArguments (optNames, dfValues, varargin);
      if (! isempty (rem))
        error ("lime.fit: unknown optional argument or misplaced value.");
      endif
      F = struct ('Cat', this.CategoricalPredictors, 'X', this.X);
      opts = limeFitOpts (Dist, KW, SMT, BT, Cov, PP, Scale, 'lime.fit', F);
      opts = limeMergeOpts (this.Opts_, opts);

      this.QueryPoint = QueryPoint;
      this.NumImportantPredictors = double (NumImportant);

      ## A local draw follows the query point, so it is taken again wherever
      ## the point has moved away from the one it was drawn around
      if (strcmp (this.DataLocality, 'local') && ! this.Custom_)
        this.SyntheticData = limeSynth (this);
        this.Fitted = limePredict (this.BlackboxModel, ...
                                   this.SyntheticData, this.Type);
      endif

      ## Nearness to the query point, and the weight it carries
      d = limeDistance (this.SyntheticData, QueryPoint, opts, ...
                        this.CategoricalPredictors, 'lime.fit');
      spread = limeSpread (this.SyntheticData, opts, ...
                           this.CategoricalPredictors, 'lime.fit');
      w = limeWeights (d, spread, opts.KernelWidth);

      ## A classifier is explained one class at a time: the simple model
      ## separates the class the explained model predicted from every other
      bbf = limePredict (this.BlackboxModel, QueryPoint, this.Type);
      this.BlackboxFitted = bbf;
      if (strcmp (this.Type, 'classification'))
        y = 2 * double (limeSameLabel (this.Fitted, bbf)) - 1;
      else
        y = double (this.Fitted(:));
      endif

      ## The predictors the simple model is fitted on, chosen a group at a
      ## time so that the levels of one predictor are taken or left together
      [A, grp] = limeExpand (this.SyntheticData, ...
                             this.CategoricalPredictors, this.SyntheticData);
      sel = limeOMP (A, y, w, grp, NumImportant);
      if (isempty (sel))
        sel = 1;
      endif
      this.ImportantPredictors = sel(:);

      this.SimpleModel = limeSimple (this, sel, y, w, opts);
      this.SimpleModelFitted = limeSimpleFitted (this, sel, QueryPoint, bbf);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {lime} {} plot (@var{obj})
    ## @deftypefnx {lime} {@var{f} =} plot (@var{obj})
    ##
    ## Plot what the simple model says about the query point.
    ##
    ## A horizontal bar per important predictor, holding the coefficient of
    ## a linear simple model or the predictor importance of a tree, in the
    ## order the predictors come in rather than by size.  @var{f} is the
    ## figure drawn.
    ##
    ## @seealso{lime, lime.fit}
    ## @end deftypefn

    function varargout = plot (this)

      if (isempty (this.SimpleModel))
        error (strcat ("lime.plot: the simple model is not fitted; use", ...
                       " fit to compute it."));
      endif

      [val, lab, ttl] = limePlotData (this);
      hf = figure ();
      ax = axes (hf);
      barh (ax, 1:numel (val), val);
      set (ax, 'ytick', 1:numel (val), 'yticklabel', lab);
      title (ax, ttl);
      if (strcmp (limeModelType (this), 'tree'))
        xlabel (ax, 'Predictor Importance');
      else
        xlabel (ax, 'Coefficient');
      endif
      ylabel (ax, 'Predictor');

      if (nargout > 0)
        varargout{1} = hf;
      endif

    endfunction

  endmethods

  methods (Access = private)

    ## Draw the observations the simple model is fitted over.  A predictor
    ## holding numbers is drawn from a normal fitted to the observations, and
    ## one holding levels from the levels as often as they appear among them.
    function S = limeSynth (this)

      X = this.X;
      n = this.NumSyntheticData;
      if (strcmp (this.DataLocality, 'local'))
        if (isempty (this.QueryPoint))
          error (strcat ("lime: 'DataLocality' of 'local' needs a query", ...
                         " point to draw around."));
        endif
        k = min (this.NumNeighbors_, rows (X));
        idx = knnsearch (X, this.QueryPoint, 'K', k);
        X = X(idx(:),:);
      endif

      M = columns (X);
      cat = this.CategoricalPredictors;
      num = setdiff (1:M, cat);
      S = zeros (n, M);
      if (! isempty (num))
        mu = mean (X(:,num), 1);
        if (rows (X) > 1)
          sigma = cov (X(:,num));
        else
          sigma = zeros (numel (num));
        endif
        S(:,num) = mvnrnd (mu, sigma, n);
      endif
      for j = cat(:)'
        lev = unique (X(:,j));
        cnt = arrayfun (@(v) sum (X(:,j) == v), lev);
        S(:,j) = lev(randsample (numel (lev), n, true, cnt));
      endfor

    endfunction

    ## Fit the simple model over the predictors taken, weighted by nearness.
    function Mdl = limeSimple (this, sel, y, w, opts)

      cat = this.CategoricalPredictors;
      Z = this.SyntheticData(:,sel);
      catsel = find (ismember (sel, cat));
      isclass = strcmp (this.Type, 'classification');
      args = {};
      if (! isempty (catsel))
        args = {'CategoricalPredictors', catsel};
      endif
      names = this.PredictorNames_(sel);

      if (strcmp (opts.SimpleModelType, 'tree'))
        if (isclass)
          Mdl = fitctree (Z, y, 'Weights', w, 'PredictorNames', names, args{:});
        else
          Mdl = fitrtree (Z, y, 'Weights', w, 'PredictorNames', names, args{:});
        endif
        return;
      endif

      ## A predictor holding levels is coded here rather than by the
      ## learner, which gives a column to every level where one is left out
      ## as the level the others are read against, as MATLAB leaves it
      [A, grp] = limeExpand (this.SyntheticData(:,sel), catsel, ...
                             this.SyntheticData(:,sel));
      enames = cell (1, columns (A));
      for k = 1:columns (A)
        enames{k} = names{grp(k)};
      endfor
      if (isclass)
        Mdl = fitclinear (A, y, 'Lambda', 0, 'Learner', 'logistic', ...
                          'Weights', w, 'BetaTolerance', opts.BetaTolerance, ...
                          'PredictorNames', enames);
      else
        Mdl = fitrlinear (A, y, 'Lambda', 0, 'Learner', 'leastsquares', ...
                          'Weights', w, 'BetaTolerance', opts.BetaTolerance, ...
                          'PredictorNames', enames);
      endif

    endfunction

    ## What a plot of the simple model draws, and what it is titled.
    function [val, lab, ttl] = limePlotData (this)

      sel = this.ImportantPredictors(:)';
      if (strcmp (limeModelType (this), 'tree'))
        val = predictorImportance (this.SimpleModel)(:);
        lab = this.PredictorNames_(sel);
        ttl = 'LIME with Decision Tree Model';
        return;
      endif
      val = this.SimpleModel.Beta(:);
      ttl = 'LIME with Linear Model';
      [~, grp] = limeExpand (this.SyntheticData(:,sel), ...
                             find (ismember (sel, ...
                                   this.CategoricalPredictors)), ...
                             this.SyntheticData(:,sel));
      lab = cell (1, numel (val));
      for k = 1:numel (val)
        lab{k} = this.PredictorNames_{sel(grp(k))};
      endfor

    endfunction

  endmethods

endclassdef

## The model, the observations and what the model answers with.
function [F, errmsg] = limeFrame (blackbox, Data, CatPred, Type)

  F = [];
  errmsg = '';
  isfh = is_function_handle (blackbox);
  if (isfh)
    props = {};
  else
    props = properties (blackbox);
  endif
  has = @(n) any (strcmp (props, n));

  ## The observations, from the model where it kept them
  if (isempty (Data))
    if (isfh)
      errmsg = "X is required when the model is a function handle.";
      return;
    endif
    if (! (has ('X') && ! isempty (blackbox.X)))
      errmsg = strcat ("X is required for a model that does not keep the", ...
                       " observations it was fitted on.");
      return;
    endif
    Data = blackbox.X;
  endif
  if (! (isnumeric (Data) && isreal (Data) && ismatrix (Data)
         && ndims (Data) == 2 && ! isempty (Data)))
    errmsg = "X must be a real numeric matrix.";
    return;
  endif
  p = columns (Data);

  ## What the model answers with
  if (isfh)
    if (isempty (Type))
      errmsg = strcat ("'Type' is required when the model is a function", ...
                       " handle.");
      return;
    endif
  endif
  if (isempty (Type))
    if (isa (blackbox, 'ClassificationPartitionedModel')
        || has ('ClassNames'))
      Type = 'classification';
    else
      Type = 'regression';
    endif
  endif
  if (! (ischar (Type) && isrow (Type)
         && any (strcmpi ({'regression', 'classification'}, Type))))
    errmsg = "'Type' must be 'regression' or 'classification'.";
    return;
  endif
  Type = lower (Type);

  ## The predictor names, from the model where it has them
  if (has ('PredictorNames') && ! isempty (blackbox.PredictorNames))
    pnames = blackbox.PredictorNames(:)';
    if (numel (pnames) != p)
      errmsg = "X must have one column per predictor of the model.";
      return;
    endif
  else
    pnames = arrayfun (@(k) sprintf ('x%d', k), 1:p, 'UniformOutput', false);
  endif

  ## The categorical predictors: the model's own, or those named for a handle
  if (isfh)
    spec = CatPred;
  else
    if (has ('CategoricalPredictors'))
      spec = blackbox.CategoricalPredictors;
    else
      spec = [];
    endif
    if (! isempty (CatPred) && ! isequal (CatPred, spec))
      errmsg = strcat ("'CategoricalPredictors' must match the", ...
                       " CategoricalPredictors property of the model.");
      return;
    endif
  endif
  cat = [];
  if (! isempty (spec))
    [C, msg] = dummyCoding (Data, spec, pnames);
    if (! isempty (msg))
      errmsg = msg;
      return;
    endif
    cat = C.Index(:)';
  endif

  cnames = [];
  if (strcmp (Type, 'classification') && has ('ClassNames'))
    cnames = blackbox.ClassNames;
  endif

  F = struct ('X', double (Data), 'Type', Type, 'Cat', cat, ...
              'PredictorNames', {pnames}, 'ClassNames', {cnames});

endfunction


## What the model answers over a matrix of observations.
function v = limePredict (blackbox, Z, Type)

  if (is_function_handle (blackbox))
    v = blackbox (Z);
    if (! (isvector (v) && numel (v) == rows (Z)))
      error (strcat ("lime: the function must answer with one value for", ...
                     " each observation."));
    endif
    if (strcmp (Type, 'regression'))
      v = double (v(:));
    else
      v = v(:);
    endif
    return;
  endif
  v = predict (blackbox, Z);
  if (strcmp (Type, 'regression'))
    v = double (v(:));
  endif

endfunction

## Whether each label is the one the explained model gave the query point.
function tf = limeSameLabel (labels, one)

  if (iscellstr (labels))
    tf = strcmp (labels(:), limeLabelText (one));
  elseif (ischar (labels))
    tf = strcmp (cellstr (labels), limeLabelText (one));
  elseif (isa (labels, 'string') || isa (labels, 'categorical'))
    tf = strcmp (cellstr (labels(:)), limeLabelText (one));
  else
    tf = (labels(:) == one(1));
  endif

endfunction

## One label as the text it compares as.
function s = limeLabelText (one)

  if (iscellstr (one))
    s = one{1};
  elseif (ischar (one))
    s = one;
  elseif (isa (one, 'string') || isa (one, 'categorical'))
    c = cellstr (one);
    s = c{1};
  else
    s = one(1);
  endif

endfunction

## The options a fit is made with, checked as they are given.
function opts = limeFitOpts (Dist, KW, SMT, BT, Cov, PP, Scale, caller, F)

  opts = struct ();
  opts.Distance = [];
  if (! isempty (Dist))
    opts.Distance = limeCheckDistance (Dist, F, caller);
  endif
  opts.KernelWidth = [];
  if (! isempty (KW))
    if (! (isnumeric (KW) && isscalar (KW) && isreal (KW) && KW > 0
           && KW <= 1))
      error (strcat ("%s: 'KernelWidth' must be a scalar greater than 0", ...
                     " and not greater than 1."), caller);
    endif
    opts.KernelWidth = double (KW);
  endif
  opts.SimpleModelType = [];
  if (! isempty (SMT))
    opts.SimpleModelType = limeCheckWord (SMT, {'linear', 'tree'}, ...
                                          caller, 'SimpleModelType', []);
  endif
  opts.BetaTolerance = [];
  if (! isempty (BT))
    if (! (isnumeric (BT) && isscalar (BT) && isreal (BT) && BT >= 0))
      error ("%s: 'BetaTolerance' must be a nonnegative scalar.", caller);
    endif
    opts.BetaTolerance = double (BT);
  endif
  opts.Cov = Cov;
  opts.P = PP;
  opts.Scale = Scale;

endfunction

## The options of a fit, over those the explainer was built with.
function opts = limeMergeOpts (base, over)

  opts = base;
  for f = fieldnames (over)'
    if (! isempty (over.(f{1})))
      opts.(f{1}) = over.(f{1});
    endif
  endfor
  if (isempty (opts.KernelWidth))
    opts.KernelWidth = 0.75;
  endif
  if (isempty (opts.SimpleModelType))
    opts.SimpleModelType = 'linear';
  endif
  if (isempty (opts.BetaTolerance))
    opts.BetaTolerance = 1e-4;
  endif

endfunction

## The named distances, which are two lists chosen by what the predictors
## hold, and a function handle, which is ours and not MATLAB's.
function d = limeCheckDistance (Dist, F, caller)

  numeric = {'euclidean', 'squaredeuclidean', 'seuclidean', ...
             'mahalanobis', 'cityblock', 'minkowski', 'chebychev', ...
             'cosine', 'correlation', 'spearman'};
  levels = {'goodall3', 'ofd'};
  iscat = ! isempty (F.Cat);

  if (is_function_handle (Dist))
    d = limeCheckHandle (Dist, F.X, caller);
    return;
  endif
  if (isa (Dist, 'string') && isscalar (Dist))
    Dist = char (Dist);
  endif
  if (iscat)
    allowed = levels;
  else
    allowed = numeric;
  endif
  if (! (ischar (Dist) && isrow (Dist) && any (strcmpi (allowed, Dist))))
    error ("%s: 'Distance' must be one of %s, or a function handle.", ...
           caller, strjoin (strcat ("'", allowed, "'"), ', '));
  endif
  j = find (strcmpi (allowed, Dist), 1);
  d = allowed{j};

endfunction

## A distance given as a function handle is tried before it is trusted, as
## pdist tries one, so that a handle that cannot answer is caught before any
## observation is drawn.
function h = limeCheckHandle (Dist, X, caller)

  if (rows (X) < 2)
    h = Dist;
    return;
  endif
  out = [];
  try
    out = Dist (X(1,:), X(2:end,:));
  catch
    error ("%s: 'Distance' is not a usable function handle.", caller);
  end_try_catch
  if (! isequal (size (out), [rows(X) - 1, 1]))
    error (strcat ("%s: 'Distance' must answer with one column holding", ...
                   " one distance per observation."), caller);
  endif
  h = Dist;

endfunction

## How far each drawn observation lies from the query point.
function d = limeDistance (S, q, opts, cat, caller)

  Dist = opts.Distance;
  if (isempty (Dist))
    if (isempty (cat))
      Dist = 'euclidean';
    else
      Dist = 'goodall3';
    endif
  endif
  if (is_function_handle (Dist))
    d = Dist (q, S);
    if (! isequal (size (d), [rows(S), 1]))
      error (strcat ("%s: 'Distance' must answer with one column holding", ...
                     " one distance per observation."), caller);
    endif
    d = double (d);
    return;
  endif
  switch (Dist)
    case 'goodall3'
      d = limeGoodall3 (S, q);
    case 'ofd'
      d = limeOFD (S, q);
    case 'seuclidean'
      d = limeNamed (S, q, 'seuclidean', opts.Scale);
    case 'mahalanobis'
      d = limeNamed (S, q, 'mahalanobis', opts.Cov);
    case 'minkowski'
      d = limeNamed (S, q, 'minkowski', opts.P);
    otherwise
      d = pdist2 (S, q, Dist);
  endswitch
  d = double (d(:));

endfunction

## A named distance that carries a parameter of its own.
function d = limeNamed (S, q, name, par)

  if (isempty (par))
    d = pdist2 (S, q, name);
  else
    d = pdist2 (S, q, name, par);
  endif

endfunction

## The Goodall 3 measure of Boriah, Chandola and Kumar (2008): a match on a
## rare level says more than a match on a common one, and a mismatch says
## nothing at all.  MATLAB's own values differ from the published measure.
function d = limeGoodall3 (S, q)

  n = rows (S);
  K = columns (S);
  s = zeros (n, 1);
  for k = 1:K
    m = (S(:,k) == q(k));
    f = sum (m);
    if (n > 1)
      p2 = f * (f - 1) / (n * (n - 1));
    else
      p2 = 0;
    endif
    s += m * (1 - p2);
  endfor
  d = 1 - s / K;

endfunction

## The occurrence frequency measure: a mismatch between two rare levels
## says less than one between two common ones.  The dissimilarity is the
## reciprocal of the similarity less one, not one less the similarity.
function d = limeOFD (S, q)

  n = rows (S);
  K = columns (S);
  s = zeros (n, 1);
  for k = 1:K
    v = S(:,k);
    m = (v == q(k));
    fq = sum (m);
    cnt = zeros (n, 1);
    lev = unique (v);
    for l = lev(:)'
      cnt(v == l) = sum (v == l);
    endfor
    sk = ones (n, 1);
    off = ! m;
    sk(off) = 1 ./ (1 + log (n / fq) * log (n ./ cnt(off)));
    s += sk;
  endfor
  d = 1 ./ (s / K) - 1;

endfunction

## The weight a drawn observation carries, which falls away with distance
## once the distances are put on a common scale.  The scale is how far apart
## the two furthest drawn observations lie, which belongs to the draw and
## not to the query point, so moving the query point does not rescale the
## weights of everything around it.
function w = limeWeights (d, scale, kw)

  if (! (scale > 0))
    w = ones (size (d));
    return;
  endif
  w = exp (-0.5 * (d ./ scale ./ kw) .^ 2);

endfunction

## How far apart the two furthest drawn observations lie.  A draw of levels
## repeats itself, so only its distinct rows are ever the furthest apart and
## the frequencies still come from the whole of it; a draw of numbers is
## taken a block at a time, so that a large one is never held as a full
## matrix of pairs.
function m = limeSpread (S, opts, cat, caller)

  n = rows (S);
  m = 0;
  if (n < 2)
    return;
  endif

  if (! isempty (cat) || is_function_handle (opts.Distance))
    if (! isempty (cat))
      U = unique (S, 'rows');
    else
      U = S;
    endif
    for k = 1:rows (U)
      d = limeDistance (S, U(k,:), opts, cat, caller);
      m = max (m, max (d));
    endfor
    return;
  endif

  Dist = opts.Distance;
  if (isempty (Dist))
    Dist = 'euclidean';
  endif
  step = max (1, min (n, floor (2e6 / n)));
  for i = 1:step:n
    j = min (i + step - 1, n);
    switch (Dist)
      case 'seuclidean'
        d = limeNamed (S, S(i:j,:), 'seuclidean', opts.Scale);
      case 'mahalanobis'
        d = limeNamed (S, S(i:j,:), 'mahalanobis', opts.Cov);
      case 'minkowski'
        d = limeNamed (S, S(i:j,:), 'minkowski', opts.P);
      otherwise
        d = pdist2 (S, S(i:j,:), Dist);
    endswitch
    m = max (m, max (d(:)));
  endfor

endfunction

## The observations as the columns a linear model is fitted on: a predictor
## holding levels becomes one column per level beyond the first, each
## saying whether the observation carries that level rather than the first.
function [A, grp] = limeExpand (Z, cat, ref)

  M = columns (Z);
  A = [];
  grp = [];
  for j = 1:M
    if (any (cat == j))
      lev = unique (ref(:,j));
      for l = lev(2:end)'
        A = [A, double(Z(:,j) == l)];
        grp = [grp, j];
      endfor
    else
      A = [A, double(Z(:,j))];
      grp = [grp, j];
    endif
  endfor

endfunction

## Group orthogonal matching pursuit: take the predictor whose columns
## explain most of what is left unexplained, fit again over everything taken
## so far, and stop where nothing is left to explain.
function sel = limeOMP (A, y, w, grp, k)

  sw = sqrt (w(:));
  Aw = A .* sw;
  yw = double (y(:)) .* sw;
  one = sw;
  den = one' * one;
  proj = @(v) v - one * ((one' * v) / den);
  r = proj (yw);
  cand = unique (grp);
  sel = [];
  for step = 1:k
    best = 0;
    bestval = 0;
    for g = cand(:)'
      if (any (sel == g))
        continue;
      endif
      B = proj (Aw(:, grp == g));
      nb = norm (B, 'fro');
      if (! (nb > 0))
        continue;
      endif
      val = norm (B * (pinv (B) * r));
      if (val > bestval)
        bestval = val;
        best = g;
      endif
    endfor
    if (best == 0 || bestval <= 1e-12 * max (1, norm (r)))
      break;
    endif
    sel(end+1) = best;
    S = [one, Aw(:, ismember (grp, sel))];
    r = yw - S * (S \ yw);
    if (norm (r) <= 1e-12 * max (1, norm (yw)))
      break;
    endif
  endfor
  sel = sort (sel);

endfunction


## What the simple model answers at the query point, put through the same
## coding the model was fitted on.
function v = limeSimpleFitted (this, sel, q, bbf)

  if (any (strcmp (class (this.SimpleModel), ...
                   {'RegressionTree', 'ClassificationTree'})))
    z = q(sel);
  else
    catsel = find (ismember (sel, this.CategoricalPredictors));
    z = limeExpand (q(sel), catsel, this.SyntheticData(:,sel));
  endif
  v = predict (this.SimpleModel, z);
  if (strcmp (this.Type, 'regression'))
    v = double (v(1));
  else
    v = v(1);
  endif

endfunction

## Whether the simple model is a tree or a weighted sum.
function t = limeModelType (this)

  if (any (strcmp (class (this.SimpleModel), ...
                   {'RegressionTree', 'ClassificationTree'})))
    t = 'tree';
  else
    t = 'linear';
  endif

endfunction


## One of a short list of words, or the default where none was given.
function v = limeCheckWord (val, allowed, caller, name, dflt)

  if (isempty (val))
    v = dflt;
    return;
  endif
  if (! (ischar (val) && isrow (val) && any (strcmpi (allowed, val))))
    error ("%s: '%s' must be one of %s.", caller, name, ...
           strjoin (strcat ("'", allowed, "'"), ', '));
  endif
  j = find (strcmpi (allowed, val), 1);
  v = allowed{j};

endfunction

## A count, or the default where none was given.
function v = limeCheckCount (val, caller, name, dflt)

  if (isempty (val))
    v = dflt;
    return;
  endif
  if (! (isnumeric (val) && isscalar (val) && isreal (val)
         && isfinite (val) && val == fix (val) && val > 0))
    error ("%s: '%s' must be a positive integer.", caller, name);
  endif
  v = double (val);

endfunction

## One observation to explain.
function q = limeCheckPoint (val, M, caller)

  if (! (isnumeric (val) && isreal (val) && isvector (val)
         && numel (val) == M))
    error (strcat ("%s: the query point must be a real numeric vector", ...
                   " of %d predictors."), caller, M);
  endif
  q = double (val(:)');

endfunction

## The explanation, measured on R2024a 2026-09-22.  Every fixture takes a
## function handle as the model: a fitted model would make the comparison
## test that model rather than this one.
%!test  # the observations are drawn when the explainer is built
%! X = [1, 10; 2, 20; 3, 30; 4, 45; 5, 50; 6, 65];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! L = lime (f, X, 'Type', 'regression', 'NumSyntheticData', 200);
%! assert_equal (size (L.SyntheticData), [200, 2]);
%! assert_equal (size (L.Fitted), [200, 1]);
%! assert_equal (isempty (L.SimpleModel), true);
%! assert_equal (isempty (L.QueryPoint), true);

%!test  # the drawn observations follow the observations they were fitted to
%! rand ('seed', 42);
%! randn ('seed', 42);
%! X = [randn(400,1) * 2 + 5, randn(400,1) * 0.5 - 1];
%! f = @(Z) Z(:,1);
%! L = lime (f, X, 'Type', 'regression', 'NumSyntheticData', 4000);
%! assert_equal (mean (L.SyntheticData), mean (X), 0.4);
%! assert_equal (std (L.SyntheticData), std (X), 0.4);

## MATLAB parity: the coefficients over five kernel widths
%!test
%! X = [0, 0, 0; 1, 0, 0; 2, 0, 0; 3, 0, 0; 4, 0, 0; 5, 0, 0; 0, 1, 0; ...
%!      0, 2, 0; 0, 3, 0; 0, 0, 1; 0, 0, 2; 0, 0, 3; 1, 1, 1; 2, 2, 2; ...
%!      3, 3, 3; 4, 4, 4; 5, 5, 5; 1, 2, 3; 3, 2, 1; 2, 3, 1];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.5 * Z(:,3) .^ 2;
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [0, 0, 0], 3, 'KernelWidth', 0.1, 'BetaTolerance', 1e-12);
%! assert_equal (a.SimpleModel.Beta', ...
%!               [1.992780427426674, -3.006957351974981, ...
%!                0.691563828807785], 1e-5);
%! assert_equal (a.SimpleModel.Bias, -0.022034064280291, 1e-5);
%! b = fit (L, [0, 0, 0], 3, 'KernelWidth', 1, 'BetaTolerance', 1e-12);
%! assert_equal (b.SimpleModel.Beta', ...
%!               [2.211849727686012, -2.799187539649991, ...
%!                1.725666244878087], 1e-5);
%! assert_equal (b.SimpleModel.Bias, -0.941458529287273, 1e-5);

## MATLAB parity: a query point outside the drawn set, which is what says
## the weights are scaled by the spread of the draw rather than by the
## distance to the furthest observation from the query
%!test
%! X = [0, 0, 0; 1, 0, 0; 2, 0, 0; 3, 0, 0; 4, 0, 0; 5, 0, 0; 0, 1, 0; ...
%!      0, 2, 0; 0, 3, 0; 0, 0, 1; 0, 0, 2; 0, 0, 3; 1, 1, 1; 2, 2, 2; ...
%!      3, 3, 3; 4, 4, 4; 5, 5, 5; 1, 2, 3; 3, 2, 1; 2, 3, 1];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.5 * Z(:,3) .^ 2;
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [0.5, 0.5, 0.5], 3, 'BetaTolerance', 1e-12);
%! assert_equal (a.SimpleModel.Beta', ...
%!               [2.197589328509646, -2.816490786994012, ...
%!                1.701570881013631], 1e-5);
%! assert_equal (a.SimpleModel.Bias, -0.887925297450557, 1e-5);

## MATLAB parity: which predictors the pursuit takes, and in what order
%!test
%! X = [0, 0, 0; 1, 0, 0; 2, 0, 0; 3, 0, 0; 4, 0, 0; 5, 0, 0; 0, 1, 0; ...
%!      0, 2, 0; 0, 3, 0; 0, 0, 1; 0, 0, 2; 0, 0, 3; 1, 1, 1; 2, 2, 2; ...
%!      3, 3, 3; 4, 4, 4; 5, 5, 5; 1, 2, 3; 3, 2, 1; 2, 3, 1];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.5 * Z(:,3) .^ 2;
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [0, 0, 0], 1, 'BetaTolerance', 1e-12);
%! assert_equal (a.ImportantPredictors, 1);
%! assert_equal (a.SimpleModel.Beta, 1.849047756469483, 1e-5);
%! b = fit (L, [0, 0, 0], 2, 'BetaTolerance', 1e-12);
%! assert_equal (b.ImportantPredictors', [1, 2]);
%! assert_equal (b.SimpleModel.Beta', ...
%!               [2.280838361112699, -1.990784316697158], 1e-5);

%!test  # the simple model is fitted on the important predictors alone
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [3, 30, 120], 2);
%! assert_equal (numel (a.SimpleModel.Beta), 2);
%! assert_equal (a.NumImportantPredictors, 2);

%!test  # a linear model is explained by its own weights
%! X = [1, 10; 2, 20; 3, 28; 4, 45; 5, 50; 6, 65; 2.5, 25; 3.5, 35];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 7;
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [3, 28], 2, 'BetaTolerance', 1e-12);
%! assert_equal (a.SimpleModel.Beta', [2, -3], 1e-6);
%! assert_equal (a.SimpleModel.Bias, 7, 1e-5);
%! assert_equal (a.SimpleModelFitted, a.BlackboxFitted, 1e-5);

%!test  # a classifier is separated from every other class at once
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! S = [meas(1:6,:); meas(51:56,:); meas(101:106,:)];
%! L = lime (Mdl, 'CustomSyntheticData', S);
%! assert_equal (L.Type, 'classification');
%! a = fit (L, meas(1,:), 2);
%! assert_equal (class (a.SimpleModel), 'ClassificationLinear');
%! assert_equal (a.SimpleModel.ClassNames', [-1, 1]);
%! assert_equal (a.BlackboxFitted, {'setosa'});

%!test  # a tree may stand in for the weighted sum
%! X = [1, 10; 2, 20; 3, 28; 4, 45; 5, 50; 6, 65; 2.5, 25; 3.5, 35];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [3, 28], 2, 'SimpleModelType', 'tree');
%! assert_equal (class (a.SimpleModel), 'RegressionTree');

%!test  # the draw may be taken around the query point rather than over all
%! load fisheriris
%! X = meas(:,2:4);
%! f = @(Z) Z(:,1);
%! L = lime (f, X, 'Type', 'regression', 'DataLocality', 'local', ...
%!           'NumNeighbors', 20, 'NumSyntheticData', 500, ...
%!           'QueryPoint', X(1,:));
%! assert_equal (L.DataLocality, 'local');
%! assert_equal (size (L.SyntheticData), [500, 3]);
%! idx = knnsearch (X, X(1,:), 'K', 20);
%! assert_equal (mean (L.SyntheticData), mean (X(idx,:)), 0.35);

%!test  # observations may be given outright rather than drawn
%! X = [1, 10; 2, 20; 3, 28; 4, 45];
%! S = [1.5, 15; 2.5, 25; 3.5, 35];
%! f = @(Z) Z(:,1);
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', S);
%! assert_equal (L.SyntheticData, S);
%! assert_equal (L.NumSyntheticData, 3);
%! assert_equal (L.Fitted, S(:,1));

%!test  # a predictor holding levels becomes one column per level beyond one
%! C = [1, 1; 1, 2; 2, 1; 2, 2; 3, 1; 3, 2; 1, 1; 2, 2];
%! f = @(Z) 5 * (Z(:,1) == 3) + Z(:,2);
%! L = lime (f, C, 'Type', 'regression', 'CategoricalPredictors', [1, 2], ...
%!           'CustomSyntheticData', C);
%! a = fit (L, [1, 1], 1);
%! assert_equal (numel (a.SimpleModel.Beta), 2);

%!test  # the two measures for levels answer differently
%! C = [1, 1, 1; 1, 1, 2; 1, 1, 3; 1, 2, 1; 1, 2, 2; 2, 1, 1; 2, 1, 2; ...
%!      2, 2, 3; 3, 1, 1; 3, 2, 2; 1, 1, 1; 1, 1, 2; 1, 2, 3; 2, 1, 1; ...
%!      1, 1, 1; 1, 1, 1; 3, 2, 3; 2, 2, 1];
%! f = @(Z) 5 * ((Z(:,2) == 1) & (Z(:,3) == 1)) + 0.5 * Z(:,1);
%! L = lime (f, C, 'Type', 'regression', ...
%!           'CategoricalPredictors', [1, 2, 3], 'CustomSyntheticData', C);
%! a = fit (L, [1, 1, 1], 3, 'Distance', 'goodall3');
%! b = fit (L, [1, 1, 1], 3, 'Distance', 'ofd');
%! assert_equal (isequal (a.SimpleModel.Beta, b.SimpleModel.Beta), false);
%! c = fit (L, [1, 1, 1], 3);
%! assert_equal (c.SimpleModel.Beta, a.SimpleModel.Beta);

%!test  # a distance given as a function handle is ours, not MATLAB's
%! X = [1, 10; 2, 20; 3, 28; 4, 45; 5, 50; 6, 65];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [3, 28], 2, 'Distance', @(q, Z) sum (abs (Z - q), 2));
%! b = fit (L, [3, 28], 2, 'Distance', 'cityblock');
%! assert_equal (a.SimpleModel.Beta, b.SimpleModel.Beta, 1e-8);

%!test  # plot draws a bar per column of the model and gives the figure
%! X = [1, 10; 2, 20; 3, 28; 4, 45; 5, 50; 6, 65];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [3, 28], 2);
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = plot (a);
%!   assert_equal (strcmp (get (h, 'type'), 'figure'), true);
%!   ax = findobj (h, 'type', 'axes');
%!   assert_equal (get (get (ax, 'title'), 'string'), ...
%!                 'LIME with Linear Model');
%!   assert_equal (get (get (ax, 'xlabel'), 'string'), 'Coefficient');
%!   assert_equal (get (get (ax, 'ylabel'), 'string'), 'Predictor');
%!   close (h);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a tree is titled after what it is, over its predictor importance
%! X = [1, 10; 2, 20; 3, 28; 4, 45; 5, 50; 6, 65];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [3, 28], 2, 'SimpleModelType', 'tree');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = plot (a);
%!   ax = findobj (h, 'type', 'axes');
%!   assert_equal (get (get (ax, 'title'), 'string'), ...
%!                 'LIME with Decision Tree Model');
%!   assert_equal (get (get (ax, 'xlabel'), 'string'), 'Predictor Importance');
%!   close (h);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # 'P' reaches the distance that takes it
%! X = [1, 10; 2, 20; 3, 28; 4, 45; 5, 50; 6, 65; 2.5, 25; 3.5, 35];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.01 * Z(:,2) .^ 2;
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [3, 28], 2, 'Distance', 'minkowski', 'P', 3);
%! same = @(q, Z) sum (abs (Z - q) .^ 3, 2) .^ (1/3);
%! b = fit (L, [3, 28], 2, 'Distance', same);
%! assert_equal (a.SimpleModel.Beta, b.SimpleModel.Beta, 1e-8);
%! c = fit (L, [3, 28], 2, 'Distance', 'minkowski', 'P', 1);
%! assert_equal (isequal (a.SimpleModel.Beta, c.SimpleModel.Beta), false);

%!test  # 'Scale' reaches the distance that takes it
%! X = [1, 10; 2, 20; 3, 28; 4, 45; 5, 50; 6, 65; 2.5, 25; 3.5, 35];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.01 * Z(:,2) .^ 2;
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! sc = [2, 20];
%! a = fit (L, [3, 28], 2, 'Distance', 'seuclidean', 'Scale', sc);
%! same = @(q, Z) sqrt (sum (((Z - q) ./ sc) .^ 2, 2));
%! b = fit (L, [3, 28], 2, 'Distance', same);
%! assert_equal (a.SimpleModel.Beta, b.SimpleModel.Beta, 1e-8);
%! c = fit (L, [3, 28], 2, 'Distance', 'seuclidean');
%! assert_equal (isequal (a.SimpleModel.Beta, c.SimpleModel.Beta), false);

%!test  # 'Cov' reaches the distance that takes it
%! X = [1, 10; 2, 20; 3, 28; 4, 45; 5, 50; 6, 65; 2.5, 25; 3.5, 35];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.01 * Z(:,2) .^ 2;
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! CV = [4, 1; 1, 400];
%! a = fit (L, [3, 28], 2, 'Distance', 'mahalanobis', 'Cov', CV);
%! same = @(q, Z) sqrt (sum (((Z - q) / CV) .* (Z - q), 2));
%! b = fit (L, [3, 28], 2, 'Distance', same);
%! assert_equal (a.SimpleModel.Beta, b.SimpleModel.Beta, 1e-8);
%! c = fit (L, [3, 28], 2, 'Distance', 'mahalanobis');
%! assert_equal (isequal (a.SimpleModel.Beta, c.SimpleModel.Beta), false);

## MATLAB parity: a predictor that adds nothing is left out, so fewer are
## used than were asked for.  The response is one no weighted sum can fit
## exactly, so the pursuit stops on the third predictor being worthless
## rather than on there being nothing left to explain.
%!test
%! X = [1, 9, 0; 2, 3, 0; 3, 7, 0; 4, 1, 0; 5, 8, 0; 6, 2, 0; 7, 5, 0; ...
%!      8, 4, 0];
%! f = @(Z) Z(:,1) .^ 2 - 0.5 * Z(:,2);
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [3, 7, 0], 3);
%! assert_equal (a.NumImportantPredictors, 3);
%! assert_equal (a.ImportantPredictors', [1, 2]);
%! assert_equal (numel (a.SimpleModel.Beta), 2);

%!test  # the pursuit stops where there is nothing left to explain
%! X = [1, 9, 0; 2, 3, 0; 3, 7, 0; 4, 1, 0; 5, 8, 0; 6, 2, 0];
%! f = @(Z) 3 * Z(:,1);
%! L = lime (f, X, 'Type', 'regression', 'CustomSyntheticData', X);
%! a = fit (L, [3, 7, 0], 3);
%! assert_equal (a.ImportantPredictors, 1);

%!test  # a classifier is drawn as its simple model, one bar per column
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! S = [meas(1:8,:); meas(51:58,:); meas(101:108,:)];
%! a = fit (lime (Mdl, 'CustomSyntheticData', S), meas(1,:), 2);
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = plot (a);
%!   ax = findobj (h, 'type', 'axes');
%!   assert_equal (get (get (ax, 'title'), 'string'), ...
%!                 'LIME with Linear Model');
%!   assert_equal (numel (get (ax, 'yticklabel')), ...
%!                 numel (a.SimpleModel.Beta));
%!   close (h);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Input validation
%!error<lime: too few input arguments.> lime ()

%!error<lime: BLACKBOX must be a fitted model that predicts, or a function handle.> ...
%! lime (42)

%!error<lime: X is required when the model is a function handle.> ...
%! lime (@(Z) Z(:,1))

%!error<lime: 'Type' is required when the model is a function handle.> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4])

%!error<lime: X must be a real numeric matrix.> ...
%! lime (@(Z) Z(:,1), {1, 2}, 'Type', 'regression')

%!error<lime: 'UseParallel' is not implemented.> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression', 'UseParallel', true)

%!error<lime: 'DataLocality' must be one of 'global', 'local'.> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression', ...
%!       'DataLocality', 'nearby')

%!error<lime: 'NumSyntheticData' must be a positive integer.> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression', ...
%!       'NumSyntheticData', 0)

%!error<lime: 'KernelWidth' must be a scalar greater than 0 and not greater than 1.> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression', 'KernelWidth', 2)

%!error<lime: 'SimpleModelType' must be one of 'linear', 'tree'.> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression', ...
%!       'SimpleModelType', 'forest')

%!error<lime: 'CustomSyntheticData' must have one column per predictor of the model.> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression', ...
%!       'CustomSyntheticData', [1, 2, 3])

%!error<lime: unknown optional argument or misplaced value.> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression', 'NoSuchThing', 1)

%!error<lime: 'Distance' must be one of 'euclidean', 'squaredeuclidean'> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression', ...
%!       'Distance', 'goodall3')

%!error<lime: 'Distance' is not a usable function handle.> ...
%! lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression', ...
%!       'Distance', @(a) a)

## The fixture needs three observations: over two, the one distance the
## trial call should answer with is a scalar, which a scalar satisfies
%!error<lime: 'Distance' must answer with one column holding one distance per observation.>
%! lime (@(Z) Z(:,1), [1, 2; 3, 4; 5, 6], 'Type', 'regression', ...
%!       'Distance', @(a, B) 5)

%!error<lime.fit: too few input arguments.> ...
%! fit (lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression'), [1, 2])

%!error<lime.fit: NUMIMPORTANTPREDICTORS must be a positive integer.> ...
%! fit (lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression'), [1, 2], 0)

%!error<lime.fit: NUMIMPORTANTPREDICTORS must not exceed the 2 predictors of the model.> ...
%! fit (lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression'), [1, 2], 5)

%!error<lime.fit: the query point must be a real numeric vector of 2 predictors.> ...
%! fit (lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression'), [1, 2, 3], 1)

%!error<lime.plot: the simple model is not fitted; use fit to compute it.> ...
%! plot (lime (@(Z) Z(:,1), [1, 2; 3, 4], 'Type', 'regression'))
