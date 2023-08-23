## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
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
## @deftypefn  {statistics} {@var{P} =} RegressionGAM
## @deftypefnx {statistics} {@var{P} =} RegressionGAM (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{P} =} RegressionGAM (@var{X}, @var{Y}, @var{name}, @var{value})
##
## Create a Generalised additive model (GAM) Regression Model,
## RegressionGAM Object using @var{X}, @var{Y} and other additional 
## Name-Value pairs. Object of RegressionGAM can be used to store the 
## training data and values can be altered via Object to predict New Values
## and other data using built-in functions for new observations.
## 
## 
## An object of @qcode{RegressionGAM} contains the following properties :
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the GAM model.
## @item
## @code{Y} is @math{Nx1} numeric matrix containing the Response data for
## corresponding predictor data in @var{X}.
## @var{Y} must have same numbers of Rows as @var{X}.
##
## @emph{Additional parameters can be passed as name value pairs :}
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"Xfit"} @tab  must be a @math{MxP} numeric matrix of 
## query/new points to predict the response value.
## @var{Xfit} must have same numbers of columns as @var{X}.
##
## @item @tab @qcode{"formula"} @tab  model specification specified as a
## string of the form @qcode{'Y ~ terms'} where 'Y' represents the reponse
## variable and 'terms' the predictor variables. formula is used to specify
## the subset of variables for training model.
## for example @qcode{"Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3"} specifies the
## linear terms x1, x2, x3, and x4 for predictor variables, x1:x2 and x2:x3
## specifies the interaction term for x1, x2 and x2, x3 respectively.
##
## @item @tab @qcode{"responsename"} @tab Response Variable Name specified as
## a string. default value is 'Y'.
##
## @item @tab @qcode{"predictors"} @tab Predictor Variable names, specified as
## cell of string(s). the length or columns of @qcode{"predictors"} must be
## same as @qcode{"X"}. If not supplied the program will generate default
## variable names (x1, x2, ... xn) for each column in @qcode{"X"}.
##
## @item @tab @qcode{"fitstd"} @tab Logical Value 0(false) or 1(true) to
## specify flag to fit model for the standard deviation of the response
## variable.
##
## @item @tab @qcode{"interactions"} @tab
##
## @item @tab @qcode{"maxpval"} @tab Maximum p-value for detecting interaction
## terms, must be a numeric scalar between 0 to 1. Interaction terms with
## p-vale less than maxpval will be used. default value is set to 0.05.
##
## @item @tab @qcode{"catpredictors"} @tab List of categorical predictors in
## predictor data @qcode{"X"}, specified as the index of column in @qcode{"X"}.
##
## @item @tab @qcode{"weights"} @tab Observational weights specified as a
## numeric matrix with each row correspoding to the observations in @qcode{"X"}.
## @qcode{"weights"} must have same number of rows as @qcode{"X"}. Default is
## ones (size (X,1),1).
##
## @item @tab @qcode{"alpha"} @tab Significance level of the prediction
## intervals @qcode{"yInt"}. Specified as scalar in range [0,1]. This argument
## is only valid when @qcode{"fitstd"} is set true. default value is 0.05.
## for example 'alpha',0.05 return 95% prediction intervals.
##
## @item @tab @qcode{"dof"} @tab Degree of freedom to fit a third order spline.
## for fitting a spline @qcode{"dof = knots + order"}, for fitting a GAM a
## polynomial spline of degree '3' is used hence the number of knots can be
## controlled by degree of freedom, degree of freedom can be used to adjust the
## fit of the each variable. the length of @qcode{"dof"} must be same as the
## columns of @qcode{"X"}. default value is 8 for each predictor variable.
##
## @end multitable
## 
##
##
## @end itemize
## for demo use demo RegressionGAM
##
## @seealso{regress}
## @end deftypefn

  properties (Access = public)

    X         = [];
    Y         = [];
    Intercept = [];
    ppfk      = [];
    
    fitstd    = [];             ## flag to fit model for std dev
    MaxPval   = [];             ## MAx p-value for detecting interaction terms
    catpred   = [];             ## list of categorical predictors
    NosObservations = [];       ## Number of observations in training dataset
    NosPredictors   = [];       ## Number of predictors 
    Interactions    = [];       ## Nos or list of interaction term
    PredictorNames  = [];       ## predictor variable names
    ResponseName    = [];       ## response variable name
    weights         = [];       ## observational weights
 
    Formula         = [];       ## formula for GAM model
    Order           = [];
    Knots           = [];
    DoF             = [];       ## degree of freedom for fitting spline
    NosIterations   = [];
    Tol             = [];       ## positive tolerence to decide convergence
  endproperties
  
  properties (Access =  private)
    ySD  = [];
    yInt = [];
    yFit = [];
  endproperties


  methods (Access = public)
    function this = RegressionGAM (X, Y, varargin)
      if (nargin < 2 && nargin != 0) ## return empty object when no arguments
        error ("RegressionGAM: too few arguments.");
      endif

      if (nargin >= 2) ## add upper condtion for nargin
        
        ## Get training sample size and number of variables in training data
        nsample = rows (X);
        ndims_X = columns (X);
        
        if (nsample != rows (Y))
          error ("RegressionGAM: number of rows in X and Y must be equal.");
        endif

        ## remove nans from X and y
        notnans  = ! logical (sum (isnan ([Y, X]), 2));
        Y        = Y (notnans);
        X        = X (notnans, :);

        if (! isnumeric (X) || ! isfinite (X))
          error ("RegressionGAM: Invalid values in X.");
        else
          this.X = X;
          NumObsv = size (X, 1);
        endif

        ## assign Y
        this.Y  = Y;

        ## adding default values
        fitstd  = false;        ## flag to fit model for std dev
        MaxPval = 1;            ## MAx p-value for detecting interaction terms
        Formula        = [];                   ## Formula for GAM model
        catpred        = {};                   ## list of categorical predictors
        Interactions   = [];                   ## Nos or list of interaction
        DoF            = ones (1, ndims_X) * 8;## Degrees of freedom
        Order          = ones (1, ndims_X) * 3;## Order of spline
        Knots          = ones (1, ndims_X) * 5;## Knots
        Tol            = 1e-2;                 ## tolerence for convergence
        PredictorNames = [];                   ## predictor variable names
        ResponseName   = "Y";                  ## response variable name
        weights        = ones (size (X,1), 1); ## observational weights

        KDO            = 0;

        while (numel (varargin) > 0)
          switch (tolower (varargin {1}))
            case "alpha"
              Alpha = varargin{2};
              if (! (isnumeric (Alpha) && isscalar (Alpha) 
                                       && Alpha > 0 && Alpha < 1))
                error (strcat (["RegressionGAM: alpha must be numeric"], ...
                               [" and in between 0 and 1."]));
              endif
              
            case "formula"
              Formula = varargin{2};
              if (! ischar (Formula) && !islogical(Formula))
                error ("RegressionGAM: Formula must be a string.");
              endif
            case "interactions"
              tmp = varargin{2};
              if (isnumeric (tmp) && isscalar (tmp) 
                                  && tmp == fix (tmp) && tmp >= 0)
                Interactions = tmp;
              elseif (ismatrix (tmp) && islogical (tmp))
                Interactions = tmp;
              elseif (ischar (tmp) && strcmpi (tmp, "all"))
                Interactions = tmp;
              else
                error (strcat (["RegressionGAM: wrong value for "],...
                               [" interactions parameter."]));
              endif
              
            case "maxpval"
                MaxPValue = varargin{2};
                if (! isnumeric (MaxPValue) ||
                    ! isscalar (MaxPValue)  || MaxPValue < 0 || MaxPValue > 1)
                  error ("RegressionGAM: MaxPValue must range from 0 to 1.");
                endif
              
            case "knots"
                if (KDO < 2)
                  Knots = varargin{2};
                  if (! isnumeric (Knots) ||
                      ! (isscalar (Knots) || 
                         isequal (size (Knots), [1, ndims_X])))
                    error ("RegressionGAM: invalid value for Knots.");
                  endif
                  DoF = Knots + Order;
                  Order = DoF - Knots;
                  KDO += 1;
                else
                  error (strcat (["RegressionGAM: DoF and Order "],...
                                 [" have been set already."]));
                endif

            case "dof"
                if (KDO < 2)
                  DoF = varargin{2};
                  if (! isnumeric (DoF) ||
                      ! (isscalar (DoF) || isequal (size (DoF), [1, ndims_X])))
                    error ("RegressionGAM: invalid value for DoF.");
                  endif
                  Knots = DoF - Order;
                  Order = DoF - Knots;
                  KDO += 1;
                else
                  error (strcat (["RegressionGAM: Knots and Order "],...
                                 [" have been set already."]));
                endif

            case "order"
                if (KDO < 2)
                  Order = varargin{2};
                  if (! isnumeric (Order) ||
                      ! (isscalar (Order) || 
                      isequal (size (Order), [1, ndims_X])))
                    error ("RegressionGAM: invalid value for Order.");
                  endif
                  DoF = Knots + Order;
                  Knots = DoF - Order;
                  KDO += 1;
                else
                  error (strcat (["RegressionGAM: DoF and Knots "],...
                                 [" have been set already."]));
                endif
              
            case "tol"
                Tol = varargin{2};
                if (! isnumeric (Tol) || ! isscalar (Tol) || !(Tol > 0))
                  error (strcat (["RegressionGAM: Tolerence (tol)"],...
                                 [" must be a Positive scalar."]));
                endif
                
            case "responsename"
                ResponseName = varargin{2};
                if (! ischar (ResponseName))
                  error (strcat (["RegressionGAM: ResponseName must "],...
                                 [" be a string or char."]));
                endif
                
            case "predictors"
                  PredictorNames = varargin{2};
                  if (! isempty (PredictorNames))
                    if (! iscellstr (PredictorNames))
                      error (strcat (["RegressionGAM: PredictorNames must"],...
                                     [" be supplied as cellstr."]));
                    elseif (columns (PredictorNames) != columns (X))
                      error (strcat (["RegressionGAM: PredictorNames must"],...
                                     [" have same number of columns as X."]));
                    endif
                  endif

            case "fitstd"
              fitstd = varargin{2};
              if (! islogical (fitstd) && fitstd != 1 && fitstd != 0)
                error (strcat (["RegressionGAM: Value in fitstd must be"],...
                               [" logical or boolean."]));
              endif
            
            case "catpredictors"
              catpred = varargin{2};

            case "weights"
              weights = varargin{2};
              if (! isvector (weights) || rows (weights) != nsample)
                error (strcat (["RegressionGAM: Weights must be a vector"],...
                               [" with same No. of rows as X."]));
              endif

            otherwise
              error (strcat (["RegressionGAM: invalid NAME in"],...
                             [" optional pairs of arguments."]));
          endswitch
          varargin (1:2) = [];
        endwhile

        ##------assign optional parameters------##
        
        
        this.NosObservations = rows (X);     ## after removing NaNs
        this.NosPredictors   = ndims_X;
        this.ResponseName    = ResponseName;
        this.Interactions    = Interactions;
        this.Intercept       = mean (Y);
        this.fitstd  = fitstd;
        this.Formula = Formula;
        this.MaxPval = MaxPval;
        this.catpred = catpred;
        this.weights = weights;
        this.Knots   = Knots;
        this.Order   = Order;
        this.DoF     = DoF;
        

        ## Adding interaction terms to the predictor matrix for training
        if (isempty (PredictorNames))
          ## empty predictornames generate default
          for i = 1:ndims_X
            PredictorNames {i} = strcat ('x', num2str (i));
          endfor
          this.PredictorNames = PredictorNames;
        endif
        IntMat = [];
        
        if (! isempty (Formula))
          IntMat = this.parseInteractions ();
          if (isempty (IntMat))
            ## user has not Specified any Ineractions in the Formula
            ## check if specified in interactions
            if (! isempty (Interactions))
              IntMat = this.parseInteractions ();
            endif
          endif
        else
        ## Formula is empty train model for given terms and check Interactions
          if (! isempty (Interactions))
            IntMat = this.parseInteractions ();
          endif
        endif

        if (! isempty (IntMat))
          for i = 1:rows (IntMatat)
            Xint   = this.X (:,IntMat (i, 1)) .* this.X (:,IntMat (i, 1));
            this.X = [this.X, Xint]; ## adding interaction term column
          endfor
        endif
        
        ## Fit the model 
        
        converged   = false;
        num_itn     = 0;
        RSS         = zeros (1, columns (X));
        res         = Y - this.Intercept;
        
        while (! converged)
          num_itn += 1;
          ## Single cycle of backfit
            for j = 1:columns (X)

              ## Calculate residuals to fit spline
              if (num_itn > 1)
                res = res + ppval (ppfk (j), X(:, j));
              endif

              ## Fit an spline to the data
              gk = splinefit (X(:,j), res, Knots(j), 'order', Order(j));

              ## This might be wrong! We need to check this out
              RSSk (j) = abs (sum (abs (Y - ppval (gk, X(:,j)) ...
                                  - this.Intercept )) .^ 2) / nsample;
              ppfk (j) = gk;
              res = res - ppval (ppfk (j), X (:,j));
            endfor

            ## check if RSS is less than the tolerence
            if (all (abs (RSS - RSSk) <= Tol))
              converged = true;
            endif

           ## update RSS
           RSS = RSSk;
        endwhile    ## loop end
        
        this.NosIterations = num_itn;
        this.ppfk          = ppfk;
      endif ## of nargin => 2

    endfunction
  endmethods
  
  ## Helper functions
  methods (Access = private)
    
    ## Function parseInteractions to detect user given interaction
    function intMat = parseInteractions (this)
      intMat = [];
      if (islogical (this.Formula))
        if (numel (this.PredictorNames) != columns (this.Formula))
          error (strcat (["RegressionGAM: Columns in interactions must"],...
                         [" Be equal to No. of predictors"]));
        endif
        for i = 1:rows (this.Formula)
          ind = find (this.Formula(i, :) == 1);
          intMat = [intMat; ind];
        endfor
      elseif (ischar (this.Formula))
        if (strcmpi (this.Formula, "all"))
        #calculate all p*(p-1)/2 interaction terms
      
        else
          #calculate interaction matrix by Formula
          if (isempty( strfind (this.Formula, '~')))
            error (strcat (["RegressionGAM: Formula Incomplete or"],...
                           [" Invalid Syntax of Formula."]));
          endif
          
          formulaParts = strsplit(this.Formula, '~');
          responseVar = strtrim(formulaParts{1});
          predictorString = strtrim(formulaParts{2});
          if (isempty (predictorString))
            error ("RegressionGAM: Incomplete Formula.");
          endif
          
          ## Now the Formula is of Valid syntax and has atleast one predictor
          intterms = strsplit(predictorString, '+');
          in = [];
          if (!isempty (strfind (this.Formula, ':')))
            for i = 1:numel (intterms)
              if (cell2mat (strfind (intterms (i), ':')))
                inters = strtrim (strsplit (cell2mat (intterms (i)), ':'));
                if (isempty (cell2mat (inters (1))) || 
                    isempty (cell2mat (inters (2))))
                  error (strcat (["RegressionGAM: Invalid syntax or"],...
                                 [" Incomplete Values in Interaction"], ...
                                 [" Term in Formula."]));
                endif
                in = [in; (index (this.PredictorNames, inters(1)) + ...
                           index (this.PredictorNames, inters(2)))];
              endif
            endfor
          endif
          ## pass recursively to get intMat
          ## cant pass recursively coz of inheritence hence just decoding
          for i = 1:rows (this.Formula)
            ind = find (this.Formula(i, :) == 1);
            intMat = [intMat; ind];
          endfor
          #this.Formula = logical (in)  ## Bad
          #in
          #intMat = this.parseInteractions ();
        endif
      else
        error ("RegressionGAM: Invalid value in Interactions.");
      endif
    endfunction
  

    ##------ gendefault to generate default predictornames-----##
    function defs = gendefault (obj)
      p = obj.NosPredictors;
      for i = 1:p
        defs {i} = strcat ('x', num2str (i));
      endfor
    endfunction ## gendefault
    
  endmethods ## private helper fns
endclassdef


%!demo
%! # Train a RegressionGAM Model for synthetic values
%! f1 = @(x) cos (3 *x);
%! f2 = @(x) x .^ 3;
%! x1 = 2 * rand (50, 1) - 1;
%! x2 = 2 * rand (50, 1) - 1;
%! y = f1(x1) + f2(x2);
%! y = y + y .* 0.2 .* rand (50,1);
%! X = [x1, x2];
%! a = RegressionGAM (X, y, "tol", 1e-3)

## Test constructor
%!test
%! a = RegressionGAM ();
%! assert (class (a), "RegressionGAM");
%! assert ({a.X, a.Y, a.Intercept, a.ppfk}, {[],[],[],[]})
%! assert ({a.fitstd, a.MaxPval, a.catpred, a.NosObservations}, {[],[],[],[]})
%! assert ({a.NosPredictors, a.Interactions, a.PredictorNames}, {[],[],[]})
%! assert ({a.ResponseName, a.weights, a.Formula, a.Order}, {[],[],[],[]})
%! assert ({a.Knots, a.DoF, a.NosIterations, a.Tol}, {[],[],[],[]})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = [1; 2; 3; 4];
%! a = RegressionGAM (x, y);
%! assert ({a.X, a.Y}, {[1,2,3;4,5,6;7,8,9;3,2,1], [1; 2; 3; 4]})
%! assert ({a.Intercept, a.Knots, a.Order}, {2.5000, [5, 5, 5], [3, 3, 3]})
%! assert ({a.DoF, a.NosIterations, a.fitstd}, {[8, 8, 8], 2, 0})
%! assert ({a.Tol, a.NosObservations, a.NosPredictors}, {[], 4, 3})
%! assert (class (a.ppfk), "struct")
%! assert ({a.ResponseName,a.PredictorNames,a.Formula},{'Y',{'x1','x2','x3'},[]})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = [1; 2; 3; 4];
%! weights = zeros (4, 1);
%! a = RegressionGAM (x, y, "weights", weights);
%! assert ({a.X, a.Y}, {[1,2,3;4,5,6;7,8,9;3,2,1], [1; 2; 3; 4]})
%! assert ({a.Intercept, a.Knots, a.Order}, {2.5000, [5, 5, 5], [3, 3, 3]})
%! assert ({a.DoF, a.NosIterations, a.fitstd}, {[8, 8, 8], 2, 0})
%! assert ({a.Tol, a.NosObservations, a.NosPredictors}, {[], 4, 3})
%! assert (class (a.ppfk), "struct")
%! assert ( a.weights, [0;0;0;0])
%! assert ({a.ResponseName,a.PredictorNames,a.Formula},{'Y',{'x1','x2','x3'},[]})


## Test input validation
%!error<RegressionGAM: too few arguments.> RegressionGAM( ones(10,2))
%!error<RegressionGAM: number of rows in X and Y must be equal.> ...
%! RegressionGAM( ones(10,2), ones (5,1))
%!error<RegressionGAM: Invalid values in X.> ...
%! RegressionGAM ([1;2;3;'a';4], ones (5,1))
%!error<RegressionGAM: invalid NAME in optional pairs of arguments.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "some", "some")
%!error<RegressionGAM: Formula must be a string.>
%! RegressionGAM( ones(10,2), ones (10,1), "formula", {"y~x1+x2"})
%!error<RegressionGAM: Formula must be a string.>
%! RegressionGAM( ones(10,2), ones (10,1), "formula", [0, 1, 0])
%!error<RegressionGAM: Formula Incomplete or Invalid Syntax of Formula.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "formula", "something")
%!error<RegressionGAM: Incomplete Formula.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "formula", "something~")
%!error<RegressionGAM: Invalid syntax or Incomplete Values in Interaction Term in Formula.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "formula", "something~x1:")
%!error<RegressionGAM: wrong value for interactions parameter.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "interactions", "some")
%!error<RegressionGAM: wrong value for interactions parameter.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "interactions", -1)
%!error<RegressionGAM: wrong value for interactions parameter.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "interactions", [1 2 3 4])
%!error<RegressionGAM: MaxPValue must range from 0 to 1.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "maxpval", [1 2 3])
%!error<RegressionGAM: MaxPValue must range from 0 to 1.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "maxpval", -1)
%!error<RegressionGAM: MaxPValue must range from 0 to 1.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "maxpval", 2)
%!error<RegressionGAM: invalid value for Knots.> ... 
%! RegressionGAM( ones(10,2), ones (10,1), "knots", 'a')
%!error<RegressionGAM: DoF and Order have been set already.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "order", 3, "dof", 2, "knots", 5)
%!error<RegressionGAM: invalid value for DoF.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "dof", 'a')
%!error<RegressionGAM: Knots and Order have been set already.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "knots", 5, "order", 3, "dof", 2)
%!error<RegressionGAM: invalid value for Order.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "order", 'a')
%!error<RegressionGAM: DoF and Knots have been set already.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "knots", 5, "dof", 2, "order", 2)
#%!error<RegressionGAM: Tolerence (tol) must be a Positive scalar.> ...
#%! RegressionGAM( ones(10,2), ones (10,1), "tol", -1)
%!error<RegressionGAM: ResponseName must be a string or char.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "responsename", -1)
%!error<RegressionGAM: PredictorNames must be supplied as cellstr.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "predictors", -1)
%!error<RegressionGAM: PredictorNames must be supplied as cellstr.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "predictors", ['a','b','c'])
%!error<RegressionGAM: PredictorNames must have same number of columns as X.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "predictors", {'a','b','c'})
%!error<RegressionGAM: Value in fitstd must be logical or boolean.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "fitstd", 5)
%!error<RegressionGAM: Value in fitstd must be logical or boolean.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "fitstd", 'a')
%!error<RegressionGAM: Weights must be a vector with same No. of rows as X.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "weights", ones (5,1))
%!error<RegressionGAM: Weights must be a vector with same No. of rows as X.> ...
%! RegressionGAM( ones(10,2), ones (10,1), "weights", ones (5,2))
