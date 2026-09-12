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

classdef CompactTreeBagger
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactTreeBagger
  ##
  ## Compact ensemble of bagged decision trees
  ##
  ## A @code{CompactTreeBagger} object carries the trees a @code{TreeBagger}
  ## ensemble grew and what prediction needs, but not the observations it was
  ## fitted on, and so none of the out-of-bag information.  It predicts new
  ## data identically to the ensemble it came from.
  ##
  ## Create one with the @code{compact} method of a @code{TreeBagger} object.
  ## Two compact ensembles fitted on the same classes are joined with
  ## @code{combine}.
  ##
  ## @seealso{TreeBagger, TreeBagger.compact, CompactClassificationTree,
  ## CompactRegressionTree}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} Method
    ##
    ## Type of the ensemble
    ##
    ## @qcode{'classification'} or @qcode{'regression'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Method = 'classification';

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} NumTrees
    ##
    ## Number of trees
    ##
    ## A nonnegative integer, the number of trees in the ensemble.  This
    ## property is read-only.
    ##
    ## @end deftp
    NumTrees = 0;

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} Trees
    ##
    ## The trees of the ensemble
    ##
    ## A column cell array holding one @code{CompactClassificationTree} or
    ## @code{CompactRegressionTree} object per tree.  This property is
    ## read-only.
    ##
    ## @end deftp
    Trees = {};

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} ClassNames
    ##
    ## Names of the classes
    ##
    ## The classes of a classification ensemble, in the type of the response
    ## it was fitted on and in the order its scores are laid out.  Empty for a
    ## regression ensemble.  This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} DefaultYfit
    ##
    ## Prediction for an observation no tree may answer for
    ##
    ## For classification, the class of greatest prior probability, in the
    ## type of the class names, or the missing label of that type after
    ## @code{setDefaultYfit} with @qcode{''}.  For regression, the weighted
    ## mean of the training response unless set otherwise.  This property is
    ## read-only; change it with @code{setDefaultYfit}.
    ##
    ## @end deftp
    DefaultYfit = [];

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors naming the columns of the predictor
    ## data.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    TreeClassIdx = {};   # columns of ClassNames each tree's scores fill
    DefaultIndex = 0;    # index of DefaultYfit into ClassNames, 0 if missing
    DefaultScore = [];   # scores of an observation no tree may answer for
  endproperties

  methods (Hidden)

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ('\n  CompactTreeBagger\n\n');
      fprintf ('%16s: %s\n', 'Method', this.Method);
      fprintf ('%16s: %d\n', 'NumTrees', this.NumTrees);
      fprintf ('%16s: %d\n', 'NumPredictors', numel (this.PredictorNames));
      if (strcmp (this.Method, 'classification'))
        fprintf ('%16s: %s\n', 'ClassNames', ...
                 classNameListing (this.ClassNames));
      endif
      fprintf ('\n');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactTreeBagger} {@var{obj} =} CompactTreeBagger (@var{B})
    ##
    ## Create a @code{CompactTreeBagger} object.
    ##
    ## @var{B} is the @code{TreeBagger} object to compact.  The documented way
    ## to reach this constructor is the @code{compact} method.
    ##
    ## @seealso{TreeBagger, TreeBagger.compact}
    ## @end deftypefn
    function this = CompactTreeBagger (B)

      ## Input validation
      if (nargin < 1)
        error ("CompactTreeBagger: too few input arguments.");
      endif
      if (! isa (B, 'TreeBagger'))
        error ("CompactTreeBagger: B must be a TreeBagger object.");
      endif

      this.Method = B.Method;
      this.NumTrees = B.NumTrees;
      this.Trees = B.Trees;
      this.ClassNames = B.ClassNames;
      this.DefaultYfit = B.DefaultYfit;
      this.PredictorNames = B.PredictorNames;
      this.TreeClassIdx = B.TreeClassIdx;
      this.DefaultIndex = B.DefaultIndex;
      this.DefaultScore = B.DefaultScore;

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {@var{label} =} predict (@var{obj}, @var{X})
    ## @deftypefnx {CompactTreeBagger} {[@var{label}, @var{scores}] =} predict (@dots{})
    ## @deftypefnx {CompactTreeBagger} {[@var{label}, @var{scores}, @var{stdevs}] =} predict (@dots{})
    ## @deftypefnx {CompactTreeBagger} {[@var{Yfit}, @var{stdevs}] =} predict (@dots{})
    ## @deftypefnx {CompactTreeBagger} {@dots{} =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict responses with a compact bagged ensemble.
    ##
    ## For a classification ensemble, @var{label} holds the predicted class of
    ## each row of @var{X}, in the type of @code{ClassNames}.  @var{scores} is
    ## an @math{NxK} matrix, the weighted average over the trees of the class
    ## probability each tree gives, and @var{stdevs} holds their standard
    ## deviations over the trees.  The label is the class of highest score,
    ## whatever cost matrix the trees were grown with, as MATLAB documents.
    ##
    ## For a regression ensemble, @var{Yfit} is the weighted average of the
    ## trees' predictions and @var{stdevs} their standard deviation.
    ##
    ## The standard deviations are population deviations, taken over the
    ## trees that answer for the observation with their weights.  An
    ## observation no tree may answer for is given @code{DefaultYfit}, with
    ## the prior as its scores and @code{NaN} as its deviations.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Trees'} @tab @tab @qcode{'all'} (default) or a vector of
    ## indices of the trees to use.
    ## @item @qcode{'TreeWeights'} @tab @tab A nonnegative vector with one
    ## weight per tree used.  The default weighs them equally.
    ## @item @qcode{'UseInstanceForTree'} @tab @tab An @math{NxNumTrees}
    ## logical matrix saying which tree may answer for which observation.  The
    ## default lets every tree answer for every observation.
    ## @end multitable
    ##
    ## MATLAB returns the labels as a cell array of character vectors whatever
    ## the type of the response; they are returned here in the type of
    ## @code{ClassNames}, as by every other classifier in this package.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.error, TreeBagger.predict}
    ## @end deftypefn
    function [Yfit, scores, stdevs] = predict (this, X, varargin)

      if (nargin < 2)
        error ("CompactTreeBagger.predict: too few input arguments.");
      endif
      [Yfit, scores, stdevs] = bagPredict (this, X, varargin, ...
                                           'CompactTreeBagger.predict', []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {@var{err} =} error (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactTreeBagger} {@var{err} =} error (@dots{}, @var{name}, @var{value})
    ##
    ## Misclassification probability or mean squared error of the ensemble.
    ##
    ## For a classification ensemble @var{err} is the weighted share of the
    ## rows of @var{X} whose predicted class differs from @var{Y}; for a
    ## regression ensemble it is the weighted mean squared error.  By default
    ## @var{err} is a column with one element per tree, the error of the
    ## first tree, then of the first two, and so on.
    ##
    ## An observation no tree may answer for is predicted as
    ## @code{DefaultYfit} and counted; when @code{DefaultYfit} is the missing
    ## label it has no prediction and is left out.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Mode'} @tab @tab @qcode{'cumulative'} (default),
    ## @qcode{'individual'} for the error of each tree on its own, or
    ## @qcode{'ensemble'} for a single error over every tree used.
    ## @item @qcode{'Trees'} @tab @tab @qcode{'all'} (default) or a vector of
    ## indices of the trees to use, in the order they are accumulated.
    ## @item @qcode{'TreeWeights'} @tab @tab A nonnegative vector with one
    ## weight per tree used.  It may not be given in @qcode{'individual'} mode.
    ## @item @qcode{'UseInstanceForTree'} @tab @tab An @math{NxNumTrees}
    ## logical matrix saying which tree may answer for which observation.
    ## @item @qcode{'Weights'} @tab @tab A nonnegative vector with one weight
    ## per observation.  The default is uniform.
    ## @end multitable
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.predict,
    ## CompactTreeBagger.meanMargin, TreeBagger.oobError}
    ## @end deftypefn
    function err = error (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactTreeBagger.error: too few input arguments.");
      endif
      err = bagLoss ('error', this, X, Y, varargin, ...
                     'CompactTreeBagger.error', [], []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactTreeBagger} {@var{m} =} margin (@dots{}, @var{name}, @var{value})
    ##
    ## Classification margin of each observation.
    ##
    ## The margin is the score of the true class less the highest score among
    ## the other classes.  By default @var{m} has one row per observation and
    ## one column per tree, the margin of the first tree, then of the first
    ## two, and so on.  An observation left without a prediction, as described
    ## under @code{CompactTreeBagger.error}, has a @code{NaN} margin.
    ##
    ## @qcode{'Mode'}, @qcode{'Trees'}, @qcode{'TreeWeights'} and
    ## @qcode{'UseInstanceForTree'} are taken as by
    ## @code{CompactTreeBagger.error}.  A regression ensemble has no margins.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.meanMargin,
    ## CompactTreeBagger.error}
    ## @end deftypefn
    function m = margin (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactTreeBagger.margin: too few input arguments.");
      endif
      m = bagLoss ('margin', this, X, Y, varargin, ...
                   'CompactTreeBagger.margin', [], []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {@var{mm} =} meanMargin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactTreeBagger} {@var{mm} =} meanMargin (@dots{}, @var{name}, @var{value})
    ##
    ## Weighted mean classification margin.
    ##
    ## @var{mm} is the weighted mean over the observations of the margins
    ## @code{CompactTreeBagger.margin} returns, a row with one element per
    ## tree by default.  The Name-Value arguments are those of
    ## @code{CompactTreeBagger.error}, @qcode{'Weights'} included.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.margin}
    ## @end deftypefn
    function mm = meanMargin (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactTreeBagger.meanMargin: too few input arguments.");
      endif
      mm = bagLoss ('meanMargin', this, X, Y, varargin, ...
                    'CompactTreeBagger.meanMargin', [], []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactTreeBagger} {@var{C} =} combine (@var{C1}, @var{C2})
    ##
    ## Join two compact ensembles.
    ##
    ## @var{C} holds the trees of @var{C1} followed by those of @var{C2}.  The
    ## two must be of the same type and, for classification, have the same
    ## class names and the same default scores.  @var{C} keeps the
    ## @code{DefaultYfit} of @var{C1}.
    ##
    ## @seealso{CompactTreeBagger, TreeBagger.append}
    ## @end deftypefn
    function this = combine (this, other)

      if (nargin < 2)
        error ("CompactTreeBagger.combine: too few input arguments.");
      endif
      if (! isa (other, 'CompactTreeBagger'))
        error (strcat ("CompactTreeBagger.combine: C2 must be a", ...
                       " CompactTreeBagger object."));
      endif
      if (! strcmp (this.Method, other.Method))
        error (strcat ("CompactTreeBagger.combine: the two ensembles must", ...
                       " be of the same type."));
      endif
      if (! (isequal (this.ClassNames, other.ClassNames)
             && isequal (this.DefaultScore, other.DefaultScore)))
        error (strcat ("CompactTreeBagger.combine: the two ensembles must", ...
                       " have the same classes and priors."));
      endif
      if (! isequal (this.PredictorNames, other.PredictorNames))
        error (strcat ("CompactTreeBagger.combine: the two ensembles must", ...
                       " have the same predictors."));
      endif
      this.Trees = [this.Trees; other.Trees];
      this.TreeClassIdx = [this.TreeClassIdx; other.TreeClassIdx];
      this.NumTrees = numel (this.Trees);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactTreeBagger} {@var{C} =} setDefaultYfit (@var{obj}, @var{Yfit})
    ##
    ## Set the prediction for an observation no tree may answer for.
    ##
    ## For a classification ensemble @var{Yfit} is @qcode{'MostPopular'}, the
    ## class of greatest prior probability, or @qcode{''}, the missing label,
    ## which leaves such an observation without a prediction and out of the
    ## error and the mean margin.  A logical response has no missing label.
    ## For a regression ensemble @var{Yfit} is a numeric scalar.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.predict}
    ## @end deftypefn
    function this = setDefaultYfit (this, Yfit)

      if (nargin < 2)
        error ("CompactTreeBagger.setDefaultYfit: too few input arguments.");
      endif
      if (strcmp (this.Method, 'classification'))
        if (ischar (Yfit) && isempty (Yfit))
          if (islogical (this.ClassNames))
            error (strcat ("CompactTreeBagger.setDefaultYfit: a logical", ...
                           " response has no missing label, so YFIT must", ...
                           " be 'MostPopular'."));
          endif
          this.DefaultIndex = 0;
          this.DefaultYfit = missingLabels (this.ClassNames, 1);
        elseif (ischar (Yfit) && strcmpi (Yfit, 'MostPopular'))
          [~, this.DefaultIndex] = max (this.DefaultScore);
          this.DefaultYfit = labelsFromIndex (this.ClassNames, ...
                                              this.DefaultIndex);
        else
          error (strcat ("CompactTreeBagger.setDefaultYfit: YFIT must be", ...
                         " '' or 'MostPopular' for a classification", ...
                         " ensemble."));
        endif
      else
        if (! (isnumeric (Yfit) && isscalar (Yfit) && isreal (Yfit)))
          error (strcat ("CompactTreeBagger.setDefaultYfit: YFIT must be", ...
                         " a real numeric scalar for a regression", ...
                         " ensemble."));
        endif
        this.DefaultYfit = double (Yfit);
      endif

    endfunction

  endmethods

endclassdef

%!demo
%! ## Compact a random forest and classify new flowers with it.  The compact
%! ## ensemble keeps the trees but not the training data.
%! load fisheriris
%! rng (42);
%! C = compact (TreeBagger (30, meas, species));
%! label = predict (C, [5.0, 3.4, 1.5, 0.2; 6.7, 3.0, 5.2, 2.3])

%!test  # a compact ensemble predicts as the ensemble it came from
%! load fisheriris
%! rng (1);
%! B = TreeBagger (5, meas, species);
%! C = compact (B);
%! assert_equal (class (C), 'CompactTreeBagger');
%! [la, sa] = predict (B, meas);
%! [lc, sc] = predict (C, meas);
%! assert_equal (lc, la);
%! assert_equal (sc, sa);
%! assert_equal (C.ClassNames, B.ClassNames);
%! assert_equal (C.DefaultYfit, B.DefaultYfit);

%!test  # MATLAB parity: combining keeps the first ensemble's default
%! load fisheriris
%! rng (1);
%! C1 = setDefaultYfit (compact (TreeBagger (3, meas, species)), '');
%! C2 = compact (TreeBagger (4, meas, species));
%! C = combine (C1, C2);
%! assert_equal (C.NumTrees, 7);
%! assert_equal (C.DefaultYfit, {''});

%!test  # a missing default leaves an observation without a prediction
%! load fisheriris
%! rng (1);
%! C = setDefaultYfit (compact (TreeBagger (4, meas, species)), '');
%! U = true (150, 4);
%! U(1,:) = false;
%! label = predict (C, meas(1:2,:), 'UseInstanceForTree', U(1:2,:));
%! assert_equal (label(1), {''});
%! m = margin (C, meas(1:2,:), species(1:2), 'Mode', 'ensemble', ...
%!             'UseInstanceForTree', U(1:2,:));
%! assert_equal (isnan (m), [true; false]);
%! miss = ! strcmp (predict (C, meas(2:end,:)), species(2:end));
%! e = error (C, meas, species, 'Mode', 'ensemble', 'UseInstanceForTree', U);
%! assert_equal (e, sum (miss) / 149, 1e-15);

%!test  # 'MostPopular' restores the class of greatest prior probability
%! load fisheriris
%! rng (1);
%! C = setDefaultYfit (compact (TreeBagger (2, meas, species)), '');
%! C = setDefaultYfit (C, 'MostPopular');
%! assert_equal (C.DefaultYfit, {'setosa'});

%!test  # a regression default is any numeric scalar
%! load fisheriris
%! rng (1);
%! C = compact (TreeBagger (2, meas(:,2:4), meas(:,1), 'Method', 'regression'));
%! C = setDefaultYfit (C, 3);
%! assert_equal (C.DefaultYfit, 3);

## Test input validation
%!shared x, y, C, R
%! load fisheriris
%! x = meas;
%! y = species;
%! C = compact (TreeBagger (2, x, y));
%! R = compact (TreeBagger (2, x(:,2:4), x(:,1), 'Method', 'regression'));
%!error<CompactTreeBagger: too few input arguments.> CompactTreeBagger ()
%!error<CompactTreeBagger: B must be a TreeBagger object.> CompactTreeBagger (1)
%!error<CompactTreeBagger.predict: too few input arguments.> predict (C)
%!error<CompactTreeBagger.margin: too few input arguments.> margin (C, x)
%!error<CompactTreeBagger.meanMargin: too few input arguments.> ...
%! meanMargin (C, x)
## Octave's test.m cuts a message through its first 'error:' before matching
## it, so these patterns hold the message after the method's prefix.
%!error<too few input arguments.> error (C, x)
%!error<CompactTreeBagger.combine: too few input arguments.> combine (C)
%!error<CompactTreeBagger.combine: C2 must be a CompactTreeBagger object.> ...
%! combine (C, 1)
%!error<CompactTreeBagger.combine: the two ensembles must be of the same type.> ...
%! combine (C, R)
%!error<CompactTreeBagger.combine: the two ensembles must have the same classes and priors.> ...
%! combine (C, compact (TreeBagger (1, x, y, 'Prior', 'uniform', ...
%!          'ClassNames', {'virginica'; 'versicolor'; 'setosa'})))
%!error<CompactTreeBagger.combine: the two ensembles must have the same predictors.> ...
%! combine (C, compact (TreeBagger (1, x, y, 'PredictorNames', ...
%!                                  {'a', 'b', 'c', 'd'})))
%!error<CompactTreeBagger.setDefaultYfit: too few input arguments.> ...
%! setDefaultYfit (C)
%!error<CompactTreeBagger.setDefaultYfit: a logical response has no missing label, so YFIT must be 'MostPopular'.> ...
%! setDefaultYfit (compact (TreeBagger (1, x, strcmp (y, 'setosa'))), '')
%!error<CompactTreeBagger.setDefaultYfit: YFIT must be '' or 'MostPopular' for a classification ensemble.> ...
%! setDefaultYfit (C, 'setosa')
%!error<CompactTreeBagger.setDefaultYfit: YFIT must be a real numeric scalar for a regression ensemble.> ...
%! setDefaultYfit (R, 'mean')
