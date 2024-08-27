## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {statistics} CompactClassificationGAM
##
## A @qcode{CompactClassificationGAM} object is a compact version of a
## Generalized Additive Model, @qcode{ClassificationGAM}.
##
## The @qcode{CompactClassificationDiscriminant} does not include the training
## data resulting to a smaller classifier size, which can be used for making
## predictions from new data, but not for tasks such as cross validation.  It
## can only be created from a @qcode{ClassificationDiscriminant} model by using
## the @code{compact} object method.
##
## The available methods for a @qcode{CompactClassificationGAM} object
## are:
## @itemize
## @item
## @code{predict}
## @item
## @code{savemodel}
## @end itemize
##
## @seealso{fitcgam, compact, ClassificationGAM}
## @end deftypefn

  properties (Access = public)

    NumPredictors   = [];   # Number of predictors
    PredictorNames  = [];   # Predictor variable names
    ResponseName    = [];   # Response variable name
    ClassNames      = [];   # Names of classes in Y
    Prior           = [];   # Prior probability for each class
    Cost            = [];   # Cost of Misclassification

    ScoreTransform  = [];   # Transformation for classification scores

    Formula         = [];   # Formula for GAM model
    Interactions    = [];   # Number or matrix of interaction terms

    Knots           = [];   # Knots of spline fitting
    Order           = [];   # Order of spline fitting
    DoF             = [];   # Degrees of freedom for fitting spline

    LearningRate    = [];   # Learning rate for training
    NumIterations   = [];   # Max number of iterations for training

    BaseModel       = [];   # Base model parameters (no interactions)
    ModelwInt       = [];   # Model parameters with interactions
    IntMatrix       = [];   # Interactions matrix applied to predictor data

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
      this.Prior           = Mdl.Prior;
      this.Cost            = Mdl.Cost;

      this.ScoreTransform  = Mdl.ScoreTransform;

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

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationGAM} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationGAM} {@var{label} =} predict (@dots{}, @qcode{'IncludeInteractions'}, @var{includeInteractions})
    ## @deftypefnx {CompactClassificationGAM} {[@var{label}, @var{score}] =} predict (@dots{})
    ##
    ## Predict labels for new data using the Generalized Additive Model (GAM)
    ## stored in a CompactClassificationGAM object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the predicted
    ## labels for the data in @var{XC} based on the model stored in the
    ## CompactClassificationGAM object, @var{obj}.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC}, 'IncludeInteractions',
    ## @var{includeInteractions})} allows you to specify whether interaction
    ## terms should be included when making predictions.
    ##
    ## @code{[@var{label}, @var{score}] = predict (@dots{})} also returns
    ## @var{score}, which contains the predicted class scores or posterior
    ## probabilities for each observation.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactClassificationGAM} class object.
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix where each row is an
    ## observation and each column corresponds to a predictor variable.
    ## @item
    ## @var{includeInteractions} is a 'true' or 'false' indicating whether to
    ## include interaction terms in the predictions.
    ## @end itemize
    ##
    ## @seealso{CompactClassificationGAM, fitcgam}
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
    ## @deftypefn  {ClassificationGAM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationGAM object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves a ClassificationGAM
    ## object into a file defined by @var{filename}.
    ##
    ## @seealso{loadmodel, fitcgam, ClassificationGAM, cvpartition,
    ## ClassificationPartitionedModel}
    ## @end deftypefn

    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "ClassificationGAM";

      ## Create variables from model properties
      NumPredictors   = this.NumPredictors;
      PredictorNames  = this.PredictorNames;
      ResponseName    = this.ResponseName;
      ClassNames      = this.ClassNames;
      Prior           = this.Prior;
      Cost            = this.Cost;
      ScoreTransform  = this.ScoreTransform;
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
      save (fname, "classdef_name", "NumPredictors", "PredictorNames", ...
            "ResponseName", "ClassNames", "Prior", "Cost", "ScoreTransform", ...
            "Formula", "Interactions", "Knots", "Order", "DoF", "BaseModel", ...
            "ModelwInt", "IntMatrix", "LearningRate", "NumIterations");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationGAM object
      mdl = CompactClassificationGAM ();

      ## Check that fieldnames in DATA match properties in ClassificationGAM
      names = fieldnames (data);
      props = fieldnames (mdl);
      if (! isequal (sort (names), sort (props)))
        error ("ClassificationGAM.load_model: invalid model in '%s'.", filename)
      endif

      ## Copy data into object
      for i = 1:numel (props)
        mdl.(props{i}) = data.(props{i});
      endfor
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
%! CMdl = crossval (Mdl);
%!
%! whos ('Mdl', 'CMdl')

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
