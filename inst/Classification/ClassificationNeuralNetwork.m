## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
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

classdef ClassificationNeuralNetwork
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{obj} =} ClassificationNeuralNetwork (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} ClassificationNeuralNetwork (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{ClassificationNeuralNetwork} class object containing a Neural
## Network classifier model.
##
## @code{@var{obj} = ClassificationNeuralNetwork (@var{X}, @var{Y})} returns a
## ClassificationNeuralNetwork object, with @var{X} as the predictor data and @var{Y}
## containing the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the Neural Network model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have same numbers of Rows as @var{X}.
## @item
## @end itemize
##
## @code{@var{obj} = ClassificationNeuralNetwork (@dots{}, @var{name}, @var{value})}
## returns a ClassificationNeuralNetwork object with parameters specified by
## @qcode{Name-Value} pair arguments.  Type @code{help fitcnet} for more info.
##
## A @qcode{ClassificationNeuralNetwork} object, @var{obj}, stores the labelled training
## data and various parameters for the Neural Network Classifier
## model, which can be accessed in the following fields:
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Field} @tab @tab @var{Description}
##
## @item @qcode{obj.X} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Y} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.NumObservations} @tab @tab Class labels, specified as a logical or
## numeric vector, or cell array of character vectors.  Each value in @var{Y} is
## the observed class label for the corresponding row in @var{X}.
##
## @end multitable
##
## @seealso{fitcnet}
## @end deftypefn

  properties (Access = public)

    X                                 = [];     # Predictor data
    Y                                 = [];     # Class labels
    W                                 = [];     # Weights of observations used to train this model

    NumObservations                   = [];     # Number of observations in training dataset
    PredictorNames                    = [];     # Predictor variables names
    ResponseName                      = [];     # Response variable name
    RowsUsed                          = [];     # Rows used in fitting
    NumIterations                     = [];     # Number of iterations taken by optimization

  endproperties


  methods (Access = public)

    ## Class object constructor
    function this = ClassificationNeuralNetwork (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationNeuralNetwork: too few input arguments.");
      endif

      ## Get training sample size and number of variables in training data
      nsample = rows (X);                    #Number of samples in X
      ndims_X = columns (X);                 #Number of dimensions in X

      ## Check if X is a table
      if (istable(X))
        if (ischar(Y) && ismember(Y, X.Properties.VariableNames))
          ## Use the variable Y from the table as the response
          Y = X.(Y);
          X(:, Y) = [];
        elseif (isstring(Y))          ## if formula is given as input
          parts = strsplit(Y, '~');
          if (numel(parts) != 2)
              error("ClassificationNeuralNetwork: Formula must be of the form 'y ~ x1 + x2 + ...'");
          endif
          responseVar = strtrim(parts{1});
          predictorStr = strtrim(parts{2});
          predictorVars = strsplit(predictorStr, '+');

          if (!ismember(responseVar, X.Properties.VariableNames))
             error("ClassificationNeuralNetwork: Response variable not found in table.");
          endif
          for i = 1:numel(predictorVars)
              if (!ismember(predictorVars{i}, X.Properties.VariableNames))
                  error("ClassificationNeuralNetwork: Predictor variable not found in table.");
              end
          endif
          ## Extract response variable
          Y = X.(responseVar);

          ## Extract predictor variables
          X = X(:, predictorVars);

        else
          error('ClassificationNeuralNetwork: Invalid Y.');
      endif

      ## Check correspodence between predictors and response
      if (nsample != rows (Y))
        error ("ClassificationNeuralNetwork: number of rows in X and Y must be equal.");
      endif


      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "weights"

          otherwise
            error (strcat (["ClassificationNeuralNetwork: invalid parameter name"],...
                           [" in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile
    endfunction
   endmethods
endclassdef


## Test input validation for constructor


## Test input validation for predict method






