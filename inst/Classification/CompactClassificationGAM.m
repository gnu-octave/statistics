## Copyright (C) 2024-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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

classdef CompactClassificationGAM
## -*- texinfo -*-
## @deftp {statistics} CompactClassificationGAM
##
## Compact generalized additive model classification
##
## The @code{CompactClassificationGAM} class is a compact version of a
## Generalized Additive Model classifier, @code{ClassificationGAM}.  It does
## not include the training data, resulting in a smaller classifier size that
## can be used for making predictions from new data, but not for tasks such as
## cross validation.
##
## A @code{CompactClassificationGAM} object can only be created from a
## @code{ClassificationGAM} model by using the @code{compact} method.
##
## @seealso{ClassificationGAM, fitcgam}
## @end deftp

  properties (Access = public)
    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the ClassificationGAM model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors   = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames  = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName    = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} ClassNames
    ##
    ## Names of classes in the response variable
    ##
    ## An array of unique values of the response variable @var{Y}, which has the
    ## same data types as the data in @var{Y}.  This property is read-only.
    ## @qcode{ClassNames} can have any of the following datatypes:
    ##
    ## @itemize
    ## @item Cell array of character vectors
    ## @item Character array
    ## @item Logical vector
    ## @item Numeric vector
    ## @end itemize
    ##
    ## @end deftp
    ClassNames      = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Prior
    ##
    ## Prior probability for each class
    ##
    ## A 2-element numeric vector specifying the prior probabilities for each
    ## class.  The order of the elements in @qcode{Prior} corresponds to the
    ## order of the classes in @qcode{ClassNames}.  This property is read-only.
    ##
    ## @end deftp
    Prior           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Cost
    ##
    ## Cost of Misclassification
    ##
    ## A square matrix specifying the cost of misclassification of a point.
    ## @qcode{Cost(i,j)} is the cost of classifying a point into class @qcode{j}
    ## if its true class is @qcode{i} (that is, the rows correspond to the true
    ## class and the columns correspond to the predicted class).  The order of
    ## the rows and columns in @qcode{Cost} corresponds to the order of the
    ## classes in @qcode{ClassNames}.  The number of rows and columns in
    ## @qcode{Cost} is the number of unique classes in the response.  By
    ## default, @qcode{Cost(i,j) = 1} if @qcode{i != j}, and
    ## @qcode{Cost(i,j) = 0} if @qcode{i = j}.  In other words, the cost is 0
    ## for correct classification and 1 for incorrect classification.
    ##
    ## Add or change the @qcode{Cost} property using dot notation as in:
    ## @itemize
    ## @item @qcode{@var{obj}.Cost = @var{costMatrix}}
    ## @end itemize
    ##
    ## @end deftp
    Cost            = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} ScoreTransform
    ##
    ## Transformation function for classification scores
    ##
    ## Specified as a function handle for transforming the classification
    ## scores.  Add or change the @qcode{ScoreTransform} property using dot
    ## notation as in:
    ##
    ## @itemize
    ## @item @qcode{@var{obj}.ScoreTransform = 'function_name'}
    ## @item @qcode{@var{obj}.ScoreTransform = @@function_handle}
    ## @end itemize
    ##
    ## When specified as a character vector, it can be any of the following
    ## built-in functions.  Nevertheless, the @qcode{ScoreTransform} property
    ## always stores their function handle equivalent.
    ##
    ## @multitable @columnfractions 0.2 0.05 0.75
    ## @headitem @var{Value} @tab @tab @var{Description}
    ## @item @qcode{"doublelogit"} @tab @tab @math{1 ./ (1 + exp (-2 * x))}
    ## @item @qcode{"invlogit"} @tab @tab @math{1 ./ (1 + exp (-x))}
    ## @item @qcode{"ismax"} @tab @tab Sets the score for the class with the
    ## largest score to 1, and for all other classes to 0
    ## @item @qcode{"logit"} @tab @tab @math{log (x ./ (1 - x))}
    ## @item @qcode{"none"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"identity"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"sign"} @tab @tab @math{-1 for x < 0, 0 for x = 0, 1 for x > 0}
    ## @item @qcode{"symmetric"} @tab @tab @math{2 * x - 1}
    ## @item @qcode{"symmetricismax"} @tab @tab Sets the score for the class
    ## with the largest score to 1, and for all other classes to -1
    ## @item @qcode{"symmetriclogit"} @tab @tab @math{2 ./ (1 + exp (-x)) - 1}
    ## @end multitable
    ##
    ## @end deftp
    ScoreTransform  = @(x) x;

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Formula
    ##
    ## Model specification formula
    ##
    ## A character vector specifying the model formula in the form
    ## @qcode{"Y ~ terms"} where @qcode{Y} represents the response variable and
    ## @qcode{terms} specifies the predictor variables and interaction terms.
    ## This property is read-only.
    ##
    ## @end deftp
    Formula         = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Interactions
    ##
    ## Interaction terms specification
    ##
    ## A logical matrix, positive integer scalar, or character vector
    ## @qcode{"all"} specifying the interaction terms between predictor
    ## variables.  This property is read-only.
    ##
    ## @end deftp
    Interactions    = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Knots
    ##
    ## Knots for spline fitting
    ##
    ## A scalar or row vector specifying the number of knots for each predictor
    ## variable in the spline fitting.  This property is read-only.
    ##
    ## @end deftp
    Knots           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Order
    ##
    ## Order of spline fitting
    ##
    ## A scalar or row vector specifying the order of the spline for each
    ## predictor variable.  This property is read-only.
    ##
    ## @end deftp
    Order           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} DoF
    ##
    ## Degrees of freedom for spline fitting
    ##
    ## A scalar or row vector specifying the degrees of freedom for each
    ## predictor variable in the spline fitting.  This property is read-only.
    ##
    ## @end deftp
    DoF             = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} LearningRate
    ##
    ## Learning rate for gradient boosting
    ##
    ## A scalar value between 0 and 1 specifying the learning rate used in the
    ## gradient boosting algorithm.  This property is read-only.
    ##
    ## @end deftp
    LearningRate    = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} NumIterations
    ##
    ## Maximum number of iterations
    ##
    ## A positive integer specifying the maximum number of iterations for the
    ## gradient boosting algorithm.  This property is read-only.
    ##
    ## @end deftp
    NumIterations   = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} BaseModel
    ##
    ## Base model parameters
    ##
    ## A structure containing the parameters of the base model without any
    ## interaction terms.  The base model represents the generalized additive
    ## model with only the main effects (predictor terms) included.
    ## This property is read-only.
    ##
    ## @end deftp
    BaseModel       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} ModelwInt
    ##
    ## Model parameters with interactions
    ##
    ## A structure containing the parameters of the model that includes
    ## interaction terms.  This model extends the base model by adding
    ## interaction terms between predictors.  This property is read-only.
    ##
    ## @end deftp
    ModelwInt       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} IntMatrix
    ##
    ## Interaction matrix
    ##
    ## A logical matrix or matrix of column indices describing the interaction
    ## terms applied to the predictor data.  This property is read-only.
    ##
    ## @end deftp
    IntMatrix       = [];
  endproperties

  properties (Access = private, Hidden)
    STname = 'none';
  endproperties

  methods (Hidden)

    ## Class object constructor
    function this = CompactClassificationGAM (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), "ClassificationGAM"))
        error ("CompactClassificationGAM: invalid classification object.");
      endif

      ## Save properties to compact model
      this.NumPredictors   = Mdl.NumPredictors;
      this.PredictorNames  = Mdl.PredictorNames;
      this.ResponseName    = Mdl.ResponseName;
      this.ClassNames      = Mdl.ClassNames;

      this.Cost            = Mdl.Cost;
      this.Prior           = Mdl.Prior;
      this.ScoreTransform  = Mdl.ScoreTransform;
      this.STname          = Mdl.STname;

      this.Formula         = Mdl.Formula;
      this.Interactions    = Mdl.Interactions;
      this.Knots           = Mdl.Knots;
      this.Order           = Mdl.Order;
      this.DoF             = Mdl.DoF;
      this.LearningRate    = Mdl.LearningRate;
      this.NumIterations   = Mdl.NumIterations;
      this.BaseModel       = Mdl.BaseModel;
      this.ModelwInt       = Mdl.ModelwInt;
      this.IntMatrix       = Mdl.IntMatrix;

    endfunction

    ## Custom display
    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    ## Custom display
    function disp (this)
      fprintf ("\n  CompactClassificationGAM\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      if (iscellstr (this.ClassNames))
        str = repmat ({"'%s'"}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      elseif (ischar (this.ClassNames))
        str = repmat ({"'%s'"}, 1, rows (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, cellstr (this.ClassNames){:});
      else # single, double, logical
        str = repmat ({"%d"}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: %s\n", 'ClassNames', str);
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.STname);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      if (! isempty (this.Formula))
        fprintf ("%+25s: '%s'\n", 'Formula', this.Formula);
      endif
      if (! isempty (this.Interactions))
        if (ischar (this.Interactions))
          fprintf ("%+25s: '%s'\n", 'Interactions', this.Interactions);
        else
          fprintf ("%+25s: [%dx%d %s]\n", 'Interactions', ...
                   size (this.Interactions, 1), size (this.Interactions, 2), ...
                   class (this.Interactions));
        endif
      endif
    endfunction

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("Invalid () indexing for referencing values", ...
                         " in a CompactClassificationGAM object."));
        case '{}'
          error (strcat ("Invalid {} indexing for referencing values", ...
                         " in a CompactClassificationGAM object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("CompactClassificationGAM.subsref: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("CompactClassificationGAM.subsref:", ...
                           " unrecognized property: '%s'"), s.subs);
          end_try_catch
      endswitch
      ## Chained references
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

    ## Class specific subscripted assignment
    function this = subsasgn (this, s, val)
      if (numel (s) > 1)
        error (strcat ("CompactClassificationGAM.subsasgn:", ...
                       " chained subscripts not allowed."));
      endif
      switch s.type
        case '()'
          error (strcat ("Invalid () indexing for assigning values", ...
                         " to a CompactClassificationGAM object."));
        case '{}'
          error (strcat ("Invalid {} indexing for assigning values", ...
                         " to a CompactClassificationGAM object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("CompactClassificationGAM.subsasgn: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          switch (s.subs)
            case 'Cost'
              this = setCost (this, val);
            case 'ScoreTransform'
              name = "CompactClassificationGAM";
              [this.ScoreTransform, this.STname] = parseScoreTransform ...
                                                   (varargin{2}, name);
            otherwise
              error (strcat ("CompactClassificationGAM.subsasgn:", ...
                             " unrecognized or read-only property: '%s'"), ...
                             s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationGAM} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationGAM} {[@var{label}, @var{score}] =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationGAM} {[@var{label}, @var{score}] =} predict (@dots{}, @qcode{'IncludeInteractions'}, @var{includeInteractions})
    ##
    ## Predict labels for new data using the Generalized Additive Model (GAM)
    ## stored in a CompactClassificationGAM object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the predicted
    ## labels for the data in @var{XC} based on the model stored in the
    ## CompactClassificationGAM object, @var{obj}.
    ##
    ## @code{[@var{label}, @var{score}] = predict (@var{obj}, @var{XC})} also
    ## returns @var{score}, which contains the predicted class scores or
    ## posterior probabilities for each observation.
    ##
    ## @code{[@var{label}, @var{score}] = predict (@var{obj}, @var{XC},
    ## 'IncludeInteractions', @var{includeInteractions})} allows you to specify
    ## whether interaction terms should be included when making predictions.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactClassificationGAM} class object.
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix where each row is an
    ## observation and each column corresponds to a predictor variable.
    ## @item
    ## @var{includeInteractions} is a logical scalar indicating whether to
    ## include interaction terms in the predictions.
    ## @end itemize
    ##
    ## @seealso{CompactClassificationGAM, ClassificationGAM, fitcgam}
    ## @end deftypefn
    function [labels, scores] = predict (this, XC, varargin)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("CompactClassificationGAM.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactClassificationGAM.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat (["CompactClassificationGAM.predict: XC must have"], ...
                       [" the same number of features as the trained model."]));
      endif

      ## Clean XC data
      notnansf  = ! logical (sum (isnan (XC), 2));
      XC        = XC (notnansf, :);

      ## Default values for Name-Value Pairs
      incInt = ! isempty (this.IntMatrix);
      Cost = this.Cost;

      ## Parse optional arguments
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "includeinteractions"
            tmpInt = varargin{2};
            if (! islogical (tmpInt) || (tmpInt != 0 && tmpInt != 1))
              error (strcat (["CompactClassificationGAM.predict:"], ...
                     [" includeinteractions must be a logical value."]));
            endif
            ## Check model for interactions
            if (tmpInt && isempty (this.IntMatrix))
              error (strcat (["CompactClassificationGAM.predict: trained"], ...
                             [" model does not include any interactions."]));
            endif
            incInt = tmpInt;

          otherwise
            error (strcat (["CompactClassificationGAM.predict: invalid"], ...
                           [" NAME in optional pairs of arguments."]));
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Choose whether interactions must be included
      if (incInt)
        if (! isempty (this.Interactions))
          ## Append interaction terms to the predictor matrix
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = XC(:,tindex);
            Xinter = ones (rows (XC), 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append interaction terms
            XC = [XC, Xinter];
          endfor
        else
          ## Add selected predictors and interaction terms
          XN = [];
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = XC(:,tindex);
            Xinter = ones (rows (XC), 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append selected predictors and interaction terms
            XN = [XN, Xinter];
          endfor
          XC = XN;
        endif
        ## Get parameters and intercept vectors from model with interactions
        params = this.ModelwInt.Parameters;
        Interc = this.ModelwInt.Intercept;
      else
        ## Get parameters and intercept vectors from base model
        params = this.BaseModel.Parameters;
        Interc = this.BaseModel.Intercept;
      endif

      ## Predict probabilities from testing data
      scores = predict_val (params, XC, Interc);

      ## Compute the expected misclassification cost matrix
      numObservations = size (XC, 1);
      CE = zeros (numObservations, 2);

      for k = 1:2
        for i = 1:2
          CE(:, k) = CE(:, k) + scores(:, i) * Cost(k, i);
        endfor
      endfor

      ## Select the class with the minimum expected misclassification cost
      [~, minIdx] = min (CE, [], 2);
      labels = this.ClassNames (minIdx);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationGAM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactClassificationGAM object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactClassificationGAM object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable,
    ## which defines the type classification object these variables constitute.
    ## Use @code{loadmodel} in order to load a classification object into
    ## Octave's workspace.
    ##
    ## @seealso{loadmodel, fitcgam, ClassificationGAM, CompactClassificationGAM}
    ## @end deftypefn
    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "CompactClassificationGAM";

      ## Create variables from model properties
      NumPredictors   = this.NumPredictors;
      PredictorNames  = this.PredictorNames;
      ResponseName    = this.ResponseName;
      ClassNames      = this.ClassNames;
      Prior           = this.Prior;
      Cost            = this.Cost;
      ScoreTransform  = this.ScoreTransform;
      STname          = this.STname;
      Formula         = this.Formula;
      Interactions    = this.Interactions;
      Knots           = this.Knots;
      Order           = this.Order;
      DoF             = this.DoF;
      LearningRate    = this.LearningRate;
      NumIterations   = this.NumIterations;
      BaseModel       = this.BaseModel;
      ModelwInt       = this.ModelwInt;
      IntMatrix       = this.IntMatrix;

      ## Save classdef name and all model properties as individual variables
      save ("-binary", fname, "classdef_name", "NumPredictors", ...
            "PredictorNames", "ResponseName", "ClassNames", "Prior", "Cost", ...
            "ScoreTransform", "STname", "Formula", "Interactions", "Knots", ...
            "Order", "DoF", "BaseModel", "ModelwInt", "IntMatrix", ...
            "LearningRate", "NumIterations");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a CompactClassificationGAM object
      mdl = CompactClassificationGAM ();

      ## Check that fieldnames in DATA match properties in CompactClassificationGAM
      names = fieldnames (data);
      props = fieldnames (mdl);
      if (! isequal (sort (names), sort (props)))
        error ("CompactClassificationGAM.load_model: invalid model in '%s'.", ...
               filename)
      endif

      ## Copy data into object
      for i = 1:numel (props)
        mdl.(props{i}) = data.(props{i});
      endfor
    endfunction

  endmethods

  methods (Access = private)

    ## Set cost
    function this = setCost (this, Cost, gnY = [])
      if (isempty (gnY))
        [~, gnY, gY] = unique (this.Y(this.RowsUsed));
      endif
      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), "double");
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat ("CompactClassificationGAM: the number", ...
                         " of rows and columns in 'Cost' must", ...
                         " correspond to selected classes in Y."));
        endif
        this.Cost = Cost;
      endif
    endfunction

  endmethods

endclassdef

## Helper function
function scores = predict_val (params, XC, intercept)
  [nsample, ndims_X] = size (XC);
  ypred = ones (nsample, 1) * intercept;

  ## Add the remaining terms
  for j = 1:ndims_X
    ypred = ypred + ppval (params(j), XC (:,j));
  endfor

  ## Apply the sigmoid function to get probabilities
  pos_prob = 1 ./ (1 + exp (-ypred));
  neg_prob = 1 - pos_prob;

  scores = [neg_prob, pos_prob];
endfunction

%!demo
%! ## Create a generalized additive model classifier and its compact version
%! # and compare their size
%!
%! load fisheriris
%! X = meas;
%! Y = species;
%!
%! Mdl = fitcdiscr (X, Y, 'ClassNames', unique (species))
%! CMdl = crossval (Mdl)

## Test constructor
%!test
%! Mdl = CompactClassificationGAM ();
%! assert (class (Mdl), "CompactClassificationGAM")
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [0; 0; 1; 1];
%! PredictorNames = {'Feature1', 'Feature2', 'Feature3'};
%! Mdl = fitcgam (x, y, "PredictorNames", PredictorNames);
%! CMdl = compact (Mdl);
%! assert (class (CMdl), "CompactClassificationGAM");
%! assert ({CMdl.NumPredictors, CMdl.ResponseName}, {3, "Y"})
%! assert (CMdl.ClassNames, {'0'; '1'})
%! assert (CMdl.PredictorNames, PredictorNames)
%! assert (CMdl.BaseModel.Intercept, 0)
%!test
%! load fisheriris
%! inds = strcmp (species,'versicolor') | strcmp (species,'virginica');
%! X = meas(inds, :);
%! Y = species(inds, :)';
%! Y = strcmp (Y, 'virginica')';
%! Mdl = fitcgam (X, Y, 'Formula', 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3');
%! CMdl = compact (Mdl);
%! assert (class (CMdl), "CompactClassificationGAM");
%! assert ({CMdl.NumPredictors, CMdl.ResponseName}, {4, "Y"})
%! assert (CMdl.ClassNames, {'0'; '1'})
%! assert (CMdl.Formula, 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3')
%! assert (CMdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'})
%! assert (CMdl.ModelwInt.Intercept, 0)
%!test
%! X = [2, 3, 5; 4, 6, 8; 1, 2, 3; 7, 8, 9; 5, 4, 3];
%! Y = [0; 1; 0; 1; 1];
%! Mdl = fitcgam (X, Y, 'Knots', [4, 4, 4], 'Order', [3, 3, 3]);
%! CMdl = compact (Mdl);
%! assert (class (CMdl), "CompactClassificationGAM");
%! assert ({CMdl.NumPredictors, CMdl.ResponseName}, {3, "Y"})
%! assert (CMdl.ClassNames, {'0'; '1'})
%! assert (CMdl.PredictorNames, {'x1', 'x2', 'x3'})
%! assert (CMdl.Knots, [4, 4, 4])
%! assert (CMdl.Order, [3, 3, 3])
%! assert (CMdl.DoF, [7, 7, 7])
%! assert (CMdl.BaseModel.Intercept, 0.4055, 1e-1)

## Test input validation for constructor
%!error<CompactClassificationGAM: invalid classification object.> ...
%! CompactClassificationGAM (1)

## Test predict method
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10];
%! y = [1; 0; 1; 0; 1];
%! Mdl = fitcgam (x, y, "interactions", "all");
%! CMdl = compact (Mdl);
%! l = {'0'; '0'; '0'; '0'; '0'};
%! s = [0.3760, 0.6240; 0.4259, 0.5741; 0.3760, 0.6240; ...
%!      0.4259, 0.5741; 0.3760, 0.6240];
%! [labels, scores] = predict (CMdl, x);
%! assert (class (CMdl), "CompactClassificationGAM");
%! assert ({CMdl.NumPredictors, CMdl.ResponseName}, {2, "Y"})
%! assert (CMdl.ClassNames, {'1'; '0'})
%! assert (CMdl.PredictorNames, {'x1', 'x2'})
%! assert (CMdl.ModelwInt.Intercept, 0.4055, 1e-1)
%! assert (labels, l)
%! assert (scores, s, 1e-1)
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [0; 0; 1; 1];
%! interactions = [false, true, false; true, false, true; false, true, false];
%! Mdl = fitcgam (x, y, "learningrate", 0.2, "interactions", interactions);
%! CMdl = compact (Mdl);
%! [label, score] = predict (CMdl, x, "includeinteractions", true);
%! l = {'0'; '0'; '1'; '1'};
%! s = [0.5106, 0.4894; 0.5135, 0.4865; 0.4864, 0.5136; 0.4847, 0.5153];
%! assert (class (CMdl), "CompactClassificationGAM");
%! assert ({CMdl.NumPredictors, CMdl.ResponseName}, {3, "Y"})
%! assert (CMdl.ClassNames, {'0'; '1'})
%! assert (CMdl.PredictorNames, {'x1', 'x2', 'x3'})
%! assert (CMdl.ModelwInt.Intercept, 0)
%! assert (label, l)
%! assert (score, s, 1e-1)

## Test input validation for predict method
%!shared CMdl
%! Mdl = fitcgam (ones (4,2), ones (4,1));
%! CMdl = compact (Mdl);
%!error<CompactClassificationGAM.predict: too few input arguments.> ...
%! predict (CMdl)
%!error<CompactClassificationGAM.predict: XC is empty.> ...
%! predict (CMdl, [])
%!error<CompactClassificationGAM.predict: XC must have the same number of features as the trained model.> ...
%! predict (CMdl, 1)
