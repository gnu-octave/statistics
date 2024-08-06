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
    ## constructor
    function this = ClassificationGAM (X, Y, varargin)

      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationGAM: too few input arguments.");
      endif

      ## Validate X
      if (! isnumeric (X))
        error ("ClassificationGAM: X must be a numeric matrix.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error (["ClassificationGAM: number of rows ", ...
                "in X and Y must be equal."]);
      endif

      ## Assign original X and Y data
      this.X = X;
      this.Y = Y;

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      ## Set default values before parsing optional parameters
      ClassNames           = [];
      Cost                 = [];
      PredictorNames       = [];
      ResponseName         = 'Y';
      Prior                = "empirical";
      InitialLearnRateForInteractions = 1;
      InitialLearnRateForPredictors   = 1;
      NumTreesPerPredictor = 300;
      NumTreesPerInteraction = 100;

      ## Parse optional parameters
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (["ClassificationGAM: PredictorNames ", ...
                      "must be supplied as a cellstring array."]);
            elseif (columns (PredictorNames) != columns (X))
              error (["ClassificationGAM: PredictorNames ", ...
                      "must have the same number of columns as X."]);
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (["ClassificationGAM: ResponseName", ...
                       " must be a character vector."]);
            endif

          case "classnames"
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames)
                                          || islogical (ClassNames)))
              error (["ClassificationGAM: ClassNames must be a", ...
                      " cellstring, logical or numeric vector."]);
            endif
            ## Check that all class names are available in gnY
            if (iscellstr (ClassNames))
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, "UniformOutput", false))))
                error (["ClassificationGAM: not all ClassNames", ...
                        " are present in Y."]);
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, "UniformOutput", false))))
                error (["ClassificationGAM: not all ClassNames", ...
                        " are present in Y."]);
              endif
            endif

          case "prior"
            Prior = varargin{2};
            if (! ((isnumeric (Prior) && isvector (Prior)) ||
                  (strcmpi (Prior, "empirical") || strcmpi (Prior, "uniform"))))
              error (["ClassificationGAM: Prior must be either", ...
                      " a numeric vector or a character vector."]);
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error (["ClassificationGAM: Cost must be", ...
                      " a numeric square matrix."]);
            endif

          case "initiallearnrateforinteractions"
            InitialLearnRateForInteractions = varargin{2};
            if (InitialLearnRateForInteractions <= 0 && InitialLearnRateForInteractions > 1)
              error ("ClassificationGAM: InitialLearnRateForInteractions must be between 0 and 1.")
            endif

          case "initiallearnrateforpredictors"
            InitialLearnRateForPredictors = varargin{2};
            if (InitialLearnRateForPredictors <= 0 && InitialLearnRateForPredictors > 1)
              error ("ClassificationGAM: InitialLearnRateForPredictors must be between 0 and 1.")
            endif

          case "numtreesperpredictor"
            NumTreesPerPredictor = varargin{2};
            if (NumTreesPerPredictor <= 0)
              error ("ClassificationGAM: NumTreesPerPredictor must be positive integer scalar.")
            endif
          
          case "numtreesperinteraction"
            NumTreesPerInteraction = varargin{2};
            if (NumTreesPerInteraction <= 0)
              error ("ClassificationGAM: NumTreesPerInteraction must be positive integer scalar.")
            endif

          case "interactions"
            Interactions = varargin{2};
            if (isnumeric (Interactions) && isscalar (Interactions) && Interactions >= 0)
              Interactions = Interactions;
            elseif (islogical (Interactions) && ismatrix (Interactions) && columns (Interactions) == columns (X))
              Interactions = Interactions;
            elseif (ischar (Interactions) && strcmpi (Interactions, 'all'))
              Interactions = 'all';
            else
              error ("ClassificationGAM: Interactions must be a nonnegative integer, logical matrix, or 'all'.");
            endif


          otherwise
            error ("ClassificationGAM: invalid name-value arguments.");
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Generate default predictors and response variabe names (if necessary)
      if (isempty (PredictorNames))
        for i = 1:columns (X)
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
      if (! isempty (ClassNames))
        if (iscellstr (ClassNames))
          ru = find (! ismember (gnY, ClassNames));
        else
          ru = find (! ismember (glY, ClassNames));
        endif
        for i = 1:numel (ru)
          gY(gY == ru(i)) = NaN;
        endfor
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
        error ("ClassificationGAM: invalid values in X.");
      endif

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
          error (["ClassificationGAM: the elements in Prior", ...
                  " do not correspond to selected classes in Y."]);
        endif
        this.Prior = Prior ./ sum (Prior);
      endif
      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), "double");
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (["ClassificationGAM: the number of rows", ...
                    " and columns in Cost must correspond", ...
                    " to selected classes in Y."]);
        endif
        this.Cost = Cost;
      endif
    endfunction
  endmethods
endclassdef