## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
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

classdef ClassificationGAM
  properties
    X = [];                   # Predictor data
    Y = [];                   # Class labels

    NumObservations = [];     # Number of observations in training dataset
    RowsUsed        = [];     # Rows used in fitting
    PredictorNames  = [];     # Predictor variables names
    ResponseName    = [];     # Response variable name
    ClassNames      = [];     # Names of classes in Y
    Prior           = [];     # Prior probability for each class
    Cost            = [];     # Cost of Misclassification

    BinEdges        = [];     # Bin edges for each predictor
    Interactions    = [];     # Interactions between predictors
    Intercept       = [];     # Intercept term
    PairDetectionBinEdges = []; # Pair-detection bin edges
    ReasonForTermination  = []; # Reason for termination
  endproperties

  methods (Access = public)
    ## class object constructor
    function this = ClassificationGAM (X, Y, varargin)

      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationGAM: too few input arguments.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error ("ClassificationGAM: number of rows in X and Y must be equal.");
      endif

      ## Assign original X and Y data to the ClassificationGAM object
      this.X = X;
      this.Y = Y;
      
      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      ## Set default values before parsing optional parameters
      PredictorNames  = [];     # Predictor variables names
      ResponseName    = [];     # Response variable name
      ClassNames      = [];     # Names of classes in Y
      Prior           = [];     # Prior probability for each class
      Cost            = [];     # Cost of misclassification
      
      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))
          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat (["ClassificationKNN: PredictorNames must"], ...
                             [" be supplied as a cellstring array."]));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat (["ClassificationKNN: PredictorNames must"], ...
                             [" have the same number of columns as X."]));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error ("ClassificationKNN: ResponseName must be a character array.");
            endif

          case "classnames"
            ClassNames = varargin{2};
            if (! iscellstr (ClassNames))
              error (strcat (["ClassificationKNN: ClassNames must"], ...
                               [" be a cellstring array."]));
            endif
            ## Check that all class names are available in gnY
            if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                 ClassNames, "UniformOutput", false))))
              error ("ClassificationKNN: not all ClassNames are present in Y.");
            endif
          
          case "prior"
            Prior = varargin{2};
            if (! ((isnumeric (Prior) && isvector (Prior)) ||
                  (strcmpi (Prior, "empirical") || strcmpi (Prior, "uniform"))))
              error (strcat (["ClassificationKNN: Prior must be either a"], ...
                             [" numeric vector or a string."]));
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error ("ClassificationKNN: Cost must be a numeric square matrix.");
            endif
          
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Get number of variables in training data
      ndims_X = columns (X);

      ## Assign the number of predictors to the ClassificationKNN object
      this.NumPredictors = ndims_X;

      ## Generate default predictors and response variabe names (if necessary)
      if (isempty (PredictorNames))
        for i = 1:ndims_X
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = "Y";
      endif

      ## Assign predictors and response variable names
      this.PredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;

      ## Handle class names
      if (isempty (ClassNames))
        ClassNames = gnY;
      else
        ru = logical (zeros (size (Y)));
        for i = 1:numel (ClassNames)
          ac = find (strcmp (gnY, ClassNames{i}));
          ru = ru | ac;
        endfor
        X = X(ru);
        Y = Y(ru);
        gY = gY(ru);
      endif

      ## Remove missing values from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, gY]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      ## Renew groups in Y
      [gY, gnY, glY] = grp2idx (Y);
      this.ClassNames = gnY;

      ## Check X contains valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("ClassificationKNN: invalid values in X.");
      endif

      ## Assign the number of observations and their correspoding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationKNN object
      this.NumObservations = rows (X);
      this.RowsUsed = cast (RowsUsed, "double");

      ## Handle Prior and Cost
      if (strcmpi ("uniform", Prior))
        this.Prior = ones (size (gnY)) ./ numel (gnY);
      elseif (isempty (Prior) || strcmpi ("empirical", Prior))
        pr = [];
        for i = 1:numel (gnY)
          pr = [pr; sum(gY==i)];
        endfor
        this.Prior = pr ./ sum (pr);
      elseif (isnumeric (Prior))
        if (numel (gnY) != numel (Prior))
          error (strcat (["ClassificationKNN: the elements in Prior do not"], ...
                         [" correspond to selected classes in Y."]));
        endif
        this.Prior = Prior ./ sum (Prior);
      endif
      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), "double");
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat (["ClassificationKNN: the number of rows and"], ...
                         [" columns in Cost must correspond to selected"], ...
                         [" classes in Y."]));
        endif
        this.Cost = Cost;
      endif


    endfunction
  endmethods
endclassdef

## demo

## tests for constructor
## tests for input validation of constructor
## tests for predict
## tests for input validation of predict
