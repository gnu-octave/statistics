## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef RegressionGAM
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{obj} =} RegressionGAM (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} RegressionGAM (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{RegressionGAM} class object containing a Generalised Additive
## Model (GAM) for regression.
##
## A @qcode{RegressionGAM} class object can store the predictors and response
## data along with various parameters for the GAM model.  It is recommended to
## use the @code{fitrgam} function to create a @qcode{RegressionGAM} object.
##
## @code{@var{obj} = RegressionGAM (@var{X}, @var{Y})} returns an object of
## class RegressionGAM, with matrix @var{X} containing the predictor data and
## vector @var{Y} containing the continuous response data.
##
## @itemize
## @item
## @var{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the GAM model.
## @item
## @var{Y} must be @math{Nx1} numeric vector containing the response data
## corresponding to the predictor data in @var{X}. @var{Y} must have same
## number of rows as @var{X}.
## @end itemize
##
## @code{@var{obj} = RegressionGAM (@dots{}, @var{name}, @var{value})} returns
## an object of class RegressionGAM with additional properties specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"predictors"} @tab Predictor Variable names, specified as
## a row vector cell of strings with the same length as the columns in @var{X}.
## If omitted, the program will generate default variable names
## @qcode{(x1, x2, ..., xn)} for each column in @var{X}.
##
## @item @tab @qcode{"responsename"} @tab Response Variable Name, specified as
## a string.  If omitted, the default value is @qcode{"Y"}.
##
## @item @tab @qcode{"formula"} @tab a model specification given as a string in
## the form @qcode{"Y ~ terms"} where @qcode{Y} represents the reponse variable
## and @qcode{terms} the predictor variables.  The formula can be used to
## specify a subset of variables for training model.  For example:
## @qcode{"Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3"} specifies four linear terms
## for the first four columns of for predictor data, and @qcode{x1:x2} and
## @qcode{x2:x3} specify the two interaction terms for 1st-2nd and 3rd-4th
## columns respectively.  Only these terms will be used for training the model,
## but @var{X} must have at least as many columns as referenced in the formula.
## If Predictor Variable names have been defined, then the terms in the formula
## must reference to those.  When @qcode{"formula"} is specified, all terms used
## for training the model are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object as a matrix containing the column indexes for each
## term including both the predictors and the interactions used.
##
## @item @tab @qcode{"interactions"} @tab a logical matrix, a positive integer
## scalar, or the string @qcode{"all"} for defining the interactions between
## predictor variables.  When given a logical matrix, it must have the same
## number of columns as @var{X} and each row corresponds to a different
## interaction term combining the predictors indexed as @qcode{true}.  Each
## interaction term is appended as a column vector after the available predictor
## column in @var{X}.  When @qcode{"all"} is defined, then all possible
## combinations of interactions are appended in @var{X} before training.  At the
## moment, parsing a positive integer has the same effect as the @qcode{"all"}
## option.  When @qcode{"interactions"} is specified, only the interaction terms
## appended to @var{X} are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object.
##
## @item @tab @qcode{"knots"} @tab a scalar or a row vector with the same
## columns as @var{X}.  It defines the knots for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 5, hence expanded to @qcode{ones (1, columns (X)) * 5}.  You can
## parse a row vector with different number of knots for each predictor
## variable to be fitted with, although not recommended.
##
## @item @tab @qcode{"order"} @tab a scalar or a row vector with the same
## columns as @var{X}.  It defines the order of the polynomial when training the
## GAM.  As a scalar, it is expanded to a row vector.  The default values is 3,
## hence expanded to @qcode{ones (1, columns (X)) * 3}.  You can parse a row
## vector with different number of polynomial order for each predictor variable
## to be fitted with, although not recommended.
##
## @item @tab @qcode{"dof"} @tab a scalar or a row vector with the same columns
## as @var{X}.  It defines the degrees of freedom for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 8, hence expanded to @qcode{ones (1, columns (X)) * 8}.  You can
## parse a row vector with different degrees of freedom for each predictor
## variable to be fitted with, although not recommended.
##
## @item @tab @qcode{"tol"} @tab a positive scalar to set the tolerance for
## covergence during training. By defaul, it is set to @qcode{1e-3}.
## @end multitable
##
## You can parse either a @qcode{"formula"} or an @qcode{"interactions"}
## optional parameter.  Parsing both parameters will result an error.
## Accordingly, you can only pass up to two parameters among @qcode{"knots"},
## @qcode{"order"}, and @qcode{"dof"} to define the required polynomial for
## training the GAM model.
##
## @seealso{fitrgam, @@RegressionGAM/predict}
## @end deftypefn

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

    Tol             = [];       # Tolerence for convergence
  endproperties

  properties (Access = private)
    ySD  = [];
    yInt = [];
    yFit = [];
  endproperties


  methods (Access = public)
    ## Class object contructor
    function this = RegressionGAM (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("RegressionGAM: too few input arguments.");
      endif

      ## Get training sample size and number of variables in training data
      nsample = rows (X);
      ndims_X = columns (X);

      ## Check correspodence between predictors and response
      if (nsample != rows (Y))
        error ("RegressionGAM: number of rows in X and Y must be equal.");
      endif

      ## Set default values before parsing optional parameters
      PredictorNames = {};                    # Predictor variable names
      ResponseName   = [];                    # Response variable name
      Formula        = [];                    # Formula for GAM model
      Interactions   = [];                    # Interaction terms
      DoF            = ones (1, ndims_X) * 8; # Degrees of freedom
      Order          = ones (1, ndims_X) * 3; # Order of spline
      Knots          = ones (1, ndims_X) * 5; # Knots
      Tol            = 1e-3;                  # Tolerence for convergence

      ## Number of parameters for Knots, DoF, Order (maximum 2 allowed)
      KOD = 0;
      ## Number of parameters for Formula, Ineractions (maximum 1 allowed)
      F_I = 0;

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "predictors"
            PredictorNames = varargin{2};
            if (! isempty (PredictorNames))
              if (! iscellstr (PredictorNames))
                error (strcat (["RegressionGAM: PredictorNames must"], ...
                               [" be a cellstring array."]));
              elseif (columns (PredictorNames) != columns (X))
                error (strcat (["RegressionGAM: PredictorNames must"], ...
                               [" have same number of columns as X."]));
              endif
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error ("RegressionGAM: ResponseName must be a char string.");
            endif

          case "formula"
            if (F_I < 1)
              Formula = varargin{2};
              if (! ischar (Formula) && ! islogical (Formula))
                error ("RegressionGAM: Formula must be a string.");
              endif
              F_I += 1;
            else
              error ("RegressionGAM: Interactions have been already defined.");
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
                error ("RegressionGAM: invalid Interactions parameter.");
              endif
              F_I += 1;
            else
              error ("RegressionGAM: Formula has been already defined.");
            endif

          case "knots"
            if (KOD < 2)
              Knots = varargin{2};
              if (! isnumeric (Knots) || ! (isscalar (Knots) ||
                  isequal (size (Knots), [1, ndims_X])))
                error ("RegressionGAM: invalid value for Knots.");
              endif
              DoF = Knots + Order;
              Order = DoF - Knots;
              KOD += 1;
            else
              error ("RegressionGAM: DoF and Order have been set already.");
            endif

          case "order"
            if (KOD < 2)
              Order = varargin{2};
              if (! isnumeric (Order) || ! (isscalar (Order) ||
                  isequal (size (Order), [1, ndims_X])))
                error ("RegressionGAM: invalid value for Order.");
              endif
              DoF = Knots + Order;
              Knots = DoF - Order;
              KOD += 1;
            else
              error ("RegressionGAM: DoF and Knots have been set already.");
            endif

          case "dof"
            if (KOD < 2)
              DoF = varargin{2};
              if (! isnumeric (DoF) ||
                  ! (isscalar (DoF) || isequal (size (DoF), [1, ndims_X])))
                error ("RegressionGAM: invalid value for DoF.");
              endif
              Knots = DoF - Order;
              Order = DoF - Knots;
              KOD += 1;
            else
              error ("RegressionGAM: Knots and Order have been set already.");
            endif

          case "tol"
            Tol = varargin{2};
            if (! (isnumeric (Tol) && isscalar (Tol) && (Tol > 0)))
              error ("RegressionGAM: Tolerance must be a Positive scalar.");
            endif

          otherwise
            error (strcat (["RegressionGAM: invalid parameter name"],...
                           [" in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile

      ## Assign original X and Y data to the RegressionGAM object
      this.X = X;
      this.Y = Y;

      ## Remove nans from X and Y
      RowsUsed  = ! logical (sum (isnan ([Y, X]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      ## Check X and Y contain valid data
      if (! isnumeric (X) || ! isfinite (X))
        error ("RegressionGAM: invalid values in X.");
      endif
      if (! isnumeric (Y) || ! isfinite (Y))
        error ("RegressionGAM: invalid values in Y.");
      endif

      ## Assign the number of observations and their correspoding indices
      ## on the original data, which will be used for training the model,
      ## to the RegressionGAM object
      this.NumObservations = rows (X);
      this.RowsUsed = cast (RowsUsed, "double");

      ## Assign the number of original predictors to the RegressionGAM object
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

      ## Assign remaining optional parameters
      this.Formula      = Formula;
      this.Interactions = Interactions;
      this.Knots        = Knots;
      this.Order        = Order;
      this.DoF          = DoF;
      this.Tol          = Tol;

      ## Fit the basic model
      Inter = mean (Y);
      [iter, param, res, RSS] = this.fitGAM (X, Y, Inter, Knots, Order);
      this.BaseModel.Intercept  = Inter;
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
        [iter, param, res, RSS] = this.fitGAM (X, Y, Inter, Knots, Order);
        this.ModelwInt.Intercept  = Inter;
        this.ModelwInt.Parameters = param;
        this.ModelwInt.Iterations = iter;
        this.ModelwInt.Residuals  = res;
        this.ModelwInt.RSS        = RSS;
      endif

    endfunction
  endmethods

  ## Helper functions
  methods (Access = private)

    ## Determine interactions from Interactions optional parameter
    function intMat = parseInteractions (this)
      if (islogical (this.Interactions))
        ## Check that interaction matrix corresponds to predictors
        if (numel (this.PredictorNames) != columns (this.Interactions))
          error (strcat (["RegressionGAM: columns in Interactions logical"], ...
                         [" matrix must equal to the number of predictors."]));
        endif
        intMat = this.Interactions
      elseif (isnumeric (this.Interactions))
        ## Need to measure the effect of all interactions to keep the best
        ## performing. Just check that the given number is not higher than
        ## p*(p-1)/2, where p is the number of predictors.
        p = this.NumPredictors;
        if (this.Interactions > p * (p - 1) / 2)
          error (strcat (["RegressionGAM: number of interaction terms"], ...
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
        error ("RegressionGAM: invalid syntax in Formula.");
      endif
      ## Split formula and keep predictor terms
      formulaParts = strsplit (this.Formula, '~');
      ## Check there is some string after '~'
      if (numel (formulaParts) < 2)
        error ("RegressionGAM: no predictor terms in Formula.");
      endif
      predictorString = strtrim (formulaParts{2});
      if (isempty (predictorString))
        error ("RegressionGAM: no predictor terms in Formula.");
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
            error ("RegressionGAM: some predictors have not been identified.");
          endif
          ## Append to interactions matrix
          intMat = [intMat; iterms];
        endif
      endfor
      ## Check that all terms have been identified
      if (! all (sum (intMat, 2) > 0))
        error ("RegressionGAM: some terms have not been identified.");
      endif
    endfunction

    ## Fit the model
    function [iter, param, res, RSS] = fitGAM (this, X, Y, Inter, Knots, Order)
      ## Initialize variables
      converged = false;
      iter      = 0;
      RSS       = zeros (1, columns (X));
      res       = Y - Inter;
      ns        = rows (X);
      Tol       = this.Tol;
      ## Start training
      while (! converged)
        iter += 1;
        ## Single cycle of backfit
        for j = 1:columns (X)

          ## Calculate residuals to fit spline
          if (iter > 1)
            res = res + ppval (param(j), X(:, j));
          endif

          ## Fit an spline to the data
          gk = splinefit (X(:,j), res, Knots(j), "order", Order(j));

          ## This might be wrong! We need to check this out
          RSSk(j) = abs (sum (abs (Y - ppval (gk, X(:,j)) - Inter)) .^ 2) / ns;
          param(j) = gk;
          res = res - ppval (param(j), X(:,j));
        endfor

        ## Check if RSS is less than the tolerence
        if (all (abs (RSS - RSSk) <= Tol))
          converged = true;
        endif

        ## Update RSS
        RSS = RSSk;
      endwhile

    endfunction

  endmethods

endclassdef


%!demo
%! ## Train a RegressionGAM Model for synthetic values
%! f1 = @(x) cos (3 * x);
%! f2 = @(x) x .^ 3;
%! x1 = 2 * rand (50, 1) - 1;
%! x2 = 2 * rand (50, 1) - 1;
%! y = f1(x1) + f2(x2);
%! y = y + y .* 0.2 .* rand (50,1);
%! X = [x1, x2];
%! a = RegressionGAM (X, y, "tol", 1e-3)

## Test constructor
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [1; 2; 3; 4];
%! a = RegressionGAM (x, y);
%! assert ({a.X, a.Y}, {x, y})
%! assert ({a.BaseModel.Intercept}, {2.5000})
%! assert ({a.Knots, a.Order, a.DoF}, {[5, 5, 5], [3, 3, 3], [8, 8, 8]})
%! assert ({a.NumObservations, a.NumPredictors}, {4, 3})
%! assert ({a.ResponseName, a.PredictorNames}, {"Y", {"x1", "x2", "x3"}})
%! assert ({a.Formula}, {[]})
%!test
%! x = [1, 2, 3, 4; 4, 5, 6, 7; 7, 8, 9, 1; 3, 2, 1, 2];
%! y = [1; 2; 3; 4];
%! pnames = {"A", "B", "C", "D"};
%! formula = "Y ~ A + B + C + D + A:C";
%! intMat = logical ([1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1;1,0,1,0]);
%! a = RegressionGAM (x, y, "predictors", pnames, "formula", formula);
%! assert ({a.IntMatrix}, {intMat})
%! assert ({a.ResponseName, a.PredictorNames}, {"Y", pnames})
%! assert ({a.Formula}, {formula})

## Test input validation
%!error<RegressionGAM: too few input arguments.> RegressionGAM ()
%!error<RegressionGAM: too few input arguments.> RegressionGAM (ones(10,2))
%!error<RegressionGAM: number of rows in X and Y must be equal.> ...
%! RegressionGAM (ones(10,2), ones (5,1))
%!error<RegressionGAM: invalid values in X.> ...
%! RegressionGAM ([1;2;3;"a";4], ones (5,1))
%!error<RegressionGAM: invalid parameter name in optional pair arguments.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "some", "some")
%!error<RegressionGAM: Formula must be a string.>
%! RegressionGAM (ones(10,2), ones (10,1), "formula", {"y~x1+x2"})
%!error<RegressionGAM: Formula must be a string.>
%! RegressionGAM (ones(10,2), ones (10,1), "formula", [0, 1, 0])
%!error<RegressionGAM: invalid syntax in Formula.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "formula", "something")
%!error<RegressionGAM: no predictor terms in Formula.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "formula", "something~")
%!error<RegressionGAM: no predictor terms in Formula.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "formula", "something~")
%!error<RegressionGAM: some predictors have not been identified> ...
%! RegressionGAM (ones(10,2), ones (10,1), "formula", "something~x1:")
%!error<RegressionGAM: invalid Interactions parameter.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "interactions", "some")
%!error<RegressionGAM: invalid Interactions parameter.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "interactions", -1)
%!error<RegressionGAM: invalid Interactions parameter.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "interactions", [1 2 3 4])
%!error<RegressionGAM: number of interaction terms requested is larger than> ...
%! RegressionGAM (ones(10,2), ones (10,1), "interactions", 3)
%!error<RegressionGAM: Formula has been already defined.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "formula", "y ~ x1 + x2", "interactions", 1)
%!error<RegressionGAM: Interactions have been already defined.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "interactions", 1, "formula", "y ~ x1 + x2")
%!error<RegressionGAM: invalid value for Knots.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "knots", "a")
%!error<RegressionGAM: DoF and Order have been set already.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "order", 3, "dof", 2, "knots", 5)
%!error<RegressionGAM: invalid value for DoF.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "dof", 'a')
%!error<RegressionGAM: Knots and Order have been set already.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "knots", 5, "order", 3, "dof", 2)
%!error<RegressionGAM: invalid value for Order.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "order", 'a')
%!error<RegressionGAM: DoF and Knots have been set already.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "knots", 5, "dof", 2, "order", 2)
%!error<RegressionGAM: Tolerance must be a Positive scalar.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "tol", -1)
%!error<RegressionGAM: ResponseName must be a char string.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "responsename", -1)
%!error<RegressionGAM: PredictorNames must be a cellstring array.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "predictors", -1)
%!error<RegressionGAM: PredictorNames must be a cellstring array.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "predictors", ['a','b','c'])
%!error<RegressionGAM: PredictorNames must have same number of columns as X.> ...
%! RegressionGAM (ones(10,2), ones (10,1), "predictors", {'a','b','c'})
