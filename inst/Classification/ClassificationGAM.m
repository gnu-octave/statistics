classdef ClassificationGAM
  properties (Access = public)

    X         = [];         # Predictor data
    Y         = [];         # Response data
    BaseModel = [];         # Base model parameters (no interactions)
    ModelwInt = [];         # Model parameters with interactions
    IntMatrix = [];         # Interactions matrix applied to predictor data

    NumObservations = [];       # Number of observations in training dataset
    RowsUsed        = [];       # Rows used in fitting
    NumPredictors   = [];       # Number of predictors
    PredictorNames  = [];       # Predictor variable names
    ResponseName    = [];       # Response variable name

    Formula         = [];       # Formula for GAM model
    Interactions    = [];       # Number or matrix of interaction terms

    Knots           = [];       # Knots of spline fitting
    Order           = [];       # Order of spline fitting
    DoF             = [];       # Degrees of freedom for fitting spline

  endproperties

  methods (Access = public)

    function this = ClassificationGAM (X, Y, varargin)
      nsample = rows (X);
      ndims_X = columns (X);

      ## Set default values before parsing optional parameters
      PredictorNames = {};                    # Predictor variable names
      ResponseName   = [];                    # Response variable name
      Formula        = [];                    # Formula for GAM model
      Interactions   = [];                    # Interaction terms
      DoF            = ones (1, ndims_X) * 8; # Degrees of freedom
      Order          = ones (1, ndims_X) * 3; # Order of spline
      Knots          = ones (1, ndims_X) * 5; # Knots
      learning_rate  = 0.1;
      num_iterations = 100;
      Cost            = [];

      ## Number of parameters for Knots, DoF, Order (maximum 2 allowed)
      KOD = 0;
      ## Number of parameters for Formula, Ineractions (maximum 1 allowed)
      F_I = 0;

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "predictors"
            PredictorNames = varargin{2};

          case "responsename"
            ResponseName = varargin{2};

          case "formula"
            if (F_I < 1)
              Formula = varargin{2};
              if (! ischar (Formula) && ! islogical (Formula))
                error ("ClassificationGAM: Formula must be a string.");
              endif
              F_I += 1;
            else
              error ("ClassificationGAM: Interactions have been already defined.");
            endif

          case "interactions"
            if (F_I < 1)
              tmp = varargin{2};
              if (isnumeric (tmp) && isscalar (tmp)
                                  && tmp == fix (tmp) && tmp >= 0)
                Interactions = tmp;
              elseif (islogical (tmp))
                Interactions = tmp;
              elseif (ischar (tmp) && strcmpi (tmp, "all"))
                Interactions = tmp;
              else
                error ("ClassificationGAM: invalid Interactions parameter.");
              endif
              F_I += 1;
            else
              error ("ClassificationGAM: Formula has been already defined.");
            endif

          case "knots"
            if (KOD < 2)
              Knots = varargin{2};
              if (! isnumeric (Knots) || ! (isscalar (Knots) ||
                  isequal (size (Knots), [1, ndims_X])))
                error ("ClassificationGAM: invalid value for Knots.");
              endif
              DoF = Knots + Order;
              Order = DoF - Knots;
              KOD += 1;
            else
              error ("ClassificationGAM: DoF and Order have been set already.");
            endif

          case "order"
            if (KOD < 2)
              Order = varargin{2};
              if (! isnumeric (Order) || ! (isscalar (Order) ||
                  isequal (size (Order), [1, ndims_X])))
                error ("ClassificationGAM: invalid value for Order.");
              endif
              DoF = Knots + Order;
              Knots = DoF - Order;
              KOD += 1;
            else
              error ("ClassificationGAM: DoF and Knots have been set already.");
            endif

          case "dof"
            if (KOD < 2)
              DoF = varargin{2};
              if (! isnumeric (DoF) ||
                  ! (isscalar (DoF) || isequal (size (DoF), [1, ndims_X])))
                error ("ClassificationGAM: invalid value for DoF.");
              endif
              Knots = DoF - Order;
              Order = DoF - Knots;
              KOD += 1;
            else
              error ("ClassificationGAM: Knots and Order have been set already.");
            endif

          case "learning_rate"
            learning_rate = varargin{2};
          
          case "num_iterations"
            num_iterations = varargin{2};

          otherwise
            error (strcat (["ClassificationGAM: invalid parameter name"],...
                           [" in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile

      ## Assign original X and Y data to the ClassificationGAM object
      this.X = X;
      this.Y = Y;

      ## Remove nans from X and Y
      RowsUsed  = ! logical (sum (isnan ([Y, X]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      this.NumObservations = rows (X);
      this.RowsUsed = cast (RowsUsed, "double");

      ## Assign the number of original predictors to the ClassificationGAM object
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

      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), "double");
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat (["ClassificationKNN: the number of rows"], ...
                         [" and columns in 'Cost' must correspond"], ...
                         [" to selected classes in Y."]));
        endif
        this.Cost = Cost;
      endif

      ## Assign remaining optional parameters
      this.Formula      = Formula;
      this.Interactions = Interactions;
      this.Knots        = Knots;
      this.Order        = Order;
      this.DoF          = DoF;

      ## Fit the basic model
      y = (Y > min(Y)) ;
      Inter = mean (y);
      [iter, param, res, RSS, intercept] = this.fitGAM (X, Y, Inter, Knots, Order, learning_rate, num_iterations);
      this.BaseModel.Intercept  = intercept;
      this.BaseModel.Parameters = param;
      this.BaseModel.Iterations = iter;
      this.BaseModel.Residuals  = res;
      this.BaseModel.RSS        = RSS;

      ## Handle interaction terms (if given)
      if (F_I > 0)
        if (isempty (this.Formula))
          ## Analyze Interactions optional parameter
          this.IntMatrix = this.parseInteractions ();
          ## Append interaction terms to the predictor matrix
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = X(:,tindex);
            Xinter = ones (this.NumObservations, 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append interaction terms
            X = [X, Xinter];
          endfor

        else
          ## Analyze Formula optional parameter
          this.IntMatrix = this.parseFormula ();
          ## Add selected predictors and interaction terms
          XN = [];
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = X(:,tindex);
            Xinter = ones (this.NumObservations, 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append selected predictors and interaction terms
            XN = [XN, Xinter];
          endfor
          X = XN;
        endif

        ## Update length of Knots, Order, and DoF vectors to match
        ## the columns of X with the interaction terms
        Knots = ones (1, columns (X)) * Knots(1); # Knots
        Order = ones (1, columns (X)) * Order(1); # Order of spline
        DoF   = ones (1, columns (X)) * DoF(1);   # Degrees of freedom

        ## Fit the model with interactions
        [iter, param, res, RSS, intercept] = this.fitGAM (X, Y, Inter, Knots, Order, learning_rate, num_iterations);
        this.ModelwInt.Intercept  = intercept;
        this.ModelwInt.Parameters = param;
        this.ModelwInt.Iterations = iter;
        this.ModelwInt.Residuals  = res;
        this.ModelwInt.RSS        = RSS;
      endif

    endfunction

    function [label, score] = predict (this, Xfit, varargin)

      ## Clean Xfit data
      notnansf  = ! logical (sum (isnan (Xfit), 2));
      Xfit      = Xfit (notnansf, :);

      ## Default values for Name-Value Pairs
      incInt = ! isempty (this.IntMatrix);
      Cost = this.Cost;

      ## Parse optional arguments
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "includeinteractions"    
            incInt = varargin{2};           

          otherwise
            error (strcat(["ClassificationGAM.predict: invalid NAME in"], ...
                          [" optional pairs of arguments."]));
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Choose whether interactions must be included
      if (incInt)
        if (! isempty (this.Interactions))
          ## Append interaction terms to the predictor matrix
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = Xfit(:,tindex);
            Xinter = ones (rows (Xfit), 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append interaction terms
            Xfit = [Xfit, Xinter];
          endfor
        else
          ## Add selected predictors and interaction terms
          XN = [];
          for i = 1:rows (this.IntMat)
            tindex = logical (this.IntMat(i,:));
            Xterms = Xfit(:,tindex);
            Xinter = ones (rows (Xfit), 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append selected predictors and interaction terms
            XN = [XN, Xinter];
          endfor
          Xfit = XN;
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
      scores = predict_val(params, Xfit, Interc);

      ## Compute the expected misclassification cost matrix
      numClasses = size (Cost, 1);
      numObservations = size (Xfit, 1);
      CE = zeros (numObservations, numClasses);

      for k = 1:numClasses
        for i = 1:numClasses
          CE(:, k) = CE(:, k) + scores(:, i) * Cost(k, i);
        endfor
      endfor

      ## Select the class with the minimum expected misclassification cost
      [~, labels] = min (CE, [], 2);

    endfunction
  endmethods
  
  ## Helper functions
  methods (Access = private)

    ## Determine interactions from Interactions optional parameter
    function intMat = parseInteractions (this)
      if (islogical (this.Interactions))
        ## Check that interaction matrix corresponds to predictors
        if (numel (this.PredictorNames) != columns (this.Interactions))
          error (strcat (["ClassificationGAM: columns in Interactions logical"], ...
                         [" matrix must equal to the number of predictors."]));
        endif
        intMat = this.Interactions
      elseif (isnumeric (this.Interactions))
        ## Need to measure the effect of all interactions to keep the best
        ## performing. Just check that the given number is not higher than
        ## p*(p-1)/2, where p is the number of predictors.
        p = this.NumPredictors;
        if (this.Interactions > p * (p - 1) / 2)
          error (strcat (["ClassificationGAM: number of interaction terms"], ...
                         [" requested is larger than all possible"], ...
                         [" combinations of predictors in X."]));
        endif
        ## Get all combinations except all zeros
        allMat = flip (fullfact(p)([2:end],:), 2);
        ## Only keep interaction terms
        iterms = find (sum (allMat, 2) != 1);
        intMat = allMat(iterms);
      elseif (strcmpi (this.Interactions, "all"))
        ## Calculate all p*(p-1)/2 interaction terms
        allMat = flip (fullfact(p)([2:end],:), 2);
        ## Only keep interaction terms
        iterms = find (sum (allMat, 2) != 1);
        intMat = allMat(iterms);
      endif
    endfunction

    ## Determine interactions from formula
    function intMat = parseFormula (this)
      intMat = [];
      ## Check formula for syntax
      if (isempty (strfind (this.Formula, '~')))
        error ("ClassificationGAM: invalid syntax in Formula.");
      endif
      ## Split formula and keep predictor terms
      formulaParts = strsplit (this.Formula, '~');
      ## Check there is some string after '~'
      if (numel (formulaParts) < 2)
        error ("ClassificationGAM: no predictor terms in Formula.");
      endif
      predictorString = strtrim (formulaParts{2});
      if (isempty (predictorString))
        error ("ClassificationGAM: no predictor terms in Formula.");
      endif
      ## Spit additive terms (between + sign)
      aterms = strtrim (strsplit (predictorString, '+'));
      ## Process all terms
      for i = 1:numel (aterms)
        ## Find individual terms (string missing ':')
        if (isempty (strfind (aterms(i), ':'){:}))
          ## Search PredictorNames to associate with column in X
          sterms = strcmp (this.PredictorNames, aterms(i));
          ## Append to interactions matrix
          intMat = [intMat; sterms];
        else
          ## Split interaction terms (string contains ':')
          mterms = strsplit (aterms{i}, ':');
          ## Add each individual predictor to interaction term vector
          iterms = logical (zeros (1, this.NumPredictors));
          for t = 1:numel (mterms)
            iterms = iterms | strcmp (this.PredictorNames, mterms(t));
          endfor
          ## Check that all predictors have been identified
          if (sum (iterms) != t)
            error ("ClassificationGAM: some predictors have not been identified.");
          endif
          ## Append to interactions matrix
          intMat = [intMat; iterms];
        endif
      endfor
      ## Check that all terms have been identified
      if (! all (sum (intMat, 2) > 0))
        error ("ClassificationGAM: some terms have not been identified.");
      endif
    endfunction

    ## Fit the model
    function [iter, param, res, RSS, intercept] = fitGAM(this, X, Y, Inter, Knots, Order, learning_rate, num_iterations)
      % Initialize variables
      [n_samples, n_features] = size(X);
      RSS = zeros(1, n_features);
      iter = num_iterations;
      
      % Initialize model predictions with the intercept (log-odds)
      p = Inter;
      intercept = log(p / (1 - p));
      f = intercept * ones(n_samples, 1);
      
      % Initialize parameters to store weak learners for each feature
      #param = cell(1, n_features);
      
      % Start boosting iterations
      for iter = 1:num_iterations
        % Compute the gradient (negative gradient of log-loss)
        y_pred = 1 ./ (1 + exp(-f));  % Sigmoid function
        gradient = Y - y_pred;        % Negative gradient of log-loss
        
        % Initialize a variable to store predictions for this iteration
        f_new = zeros(n_samples, 1);
        
        for j = 1:n_features
          % Fit a spline to the gradient for feature X_j
          spline_model = splinefit(X(:, j), gradient, Knots(j), "order", Order(j));
          
          % Predict using the fitted spline
          spline_pred = ppval(spline_model, X(:, j));
          
          % Store the spline model parameters
          param(j) = spline_model;
          
          % Update the model predictions
          f_new = f_new + learning_rate * spline_pred;
        endfor
        
        % Update the overall model predictions
        f = f + f_new;
      endfor
      
      % Final residuals and RSS calculation
      res = Y - 1 ./ (1 + exp(-f));
      RSS = sum(res .^ 2);
    endfunction

  endmethods
endclassdef

## Helper function for making prediction of new data based on GAM model
function ypred = predict_val (params, X, intercept)
  [nsample, ndims_X] = size (X);
  ypred = ones (nsample, 1) * intercept;
  ## Add the remaining terms
  for j = 1:ndims_X
    ypred = ypred + ppval (params(j), X (:,j));
  endfor
endfunction