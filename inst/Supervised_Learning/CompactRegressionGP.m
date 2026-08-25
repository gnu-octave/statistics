## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {statistics} {@var{CMdl} =} CompactRegressionGP (@var{Mdl})
##
## Create a @qcode{CompactRegressionGP} object containing a Gaussian process
## regression model without its training data.
##
## @code{@var{CMdl} = CompactRegressionGP (@var{Mdl})} returns the compact form
## of the @qcode{RegressionGP} model @var{Mdl}, keeping what is needed to
## predict and dropping the rest.  It is what @code{compact} returns, and it is
## not usually constructed directly.
##
## A compact model keeps the active set it predicts from, the prediction
## weights, the covariance function and its parameters, the explicit basis and
## its coefficients, the noise standard deviation and the standardizing
## location and scale.  It drops the response, the observation weights, the
## rows used, the count of observations and the maximized log likelihood, so
## it can predict but cannot be cross validated, refitted, or asked for its
## resubstitution loss or its post-fit statistics.
##
## The standard deviation and the prediction intervals remain available,
## because the active set of an exactly fitted model is the whole of the
## training predictors and the factorization can be rebuilt from it.
##
## @seealso{RegressionGP, fitrgp}
## @end deftypefn

classdef CompactRegressionGP

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} PredictorNames
    ##
    ## Predictor variable names
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} ExpandedPredictorNames
    ##
    ## Expanded predictor variable names
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = 'Y';

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A vector of positive integers, or empty.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} FitMethod
    ##
    ## Method used to estimate the parameters
    ##
    ## @qcode{'Exact'} or @qcode{'None'}.  This property is read-only.
    ##
    ## @end deftp
    FitMethod             = 'Exact';

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} BasisFunction
    ##
    ## Explicit basis of the model
    ##
    ## A character vector or a function handle.  This property is read-only.
    ##
    ## @end deftp
    BasisFunction         = 'Constant';

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} Beta
    ##
    ## Estimated coefficients of the explicit basis
    ##
    ## A numeric vector, empty when the basis is @qcode{'None'}.  This property
    ## is read-only.
    ##
    ## @end deftp
    Beta                  = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} Sigma
    ##
    ## Estimated noise standard deviation
    ##
    ## A positive scalar.  This property is read-only.
    ##
    ## @end deftp
    Sigma                 = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} KernelFunction
    ##
    ## Form of the covariance function
    ##
    ## A character vector or a function handle.  This property is read-only.
    ##
    ## @end deftp
    KernelFunction        = 'SquaredExponential';

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} KernelInformation
    ##
    ## Covariance function and its parameters
    ##
    ## A structure with fields @qcode{Name}, @qcode{KernelParameters} and
    ## @qcode{KernelParameterNames}.  This property is read-only.
    ##
    ## @end deftp
    KernelInformation     = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} PredictMethod
    ##
    ## Method used to make predictions
    ##
    ## @qcode{'Exact'}.  This property is read-only.
    ##
    ## @end deftp
    PredictMethod         = 'Exact';

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} Alpha
    ##
    ## Weights the predictions are made from
    ##
    ## A numeric vector with one weight per active set vector.  This property
    ## is read-only.
    ##
    ## @end deftp
    Alpha                 = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} ActiveSetVectors
    ##
    ## Subset of the training data used for predictions
    ##
    ## An @math{MxP} numeric matrix, standardized where the model standardized
    ## its predictors.  This property is read-only.
    ##
    ## @end deftp
    ActiveSetVectors      = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} ActiveSetMethod
    ##
    ## Method used to select the active set
    ##
    ## @qcode{'Random'}.  This property is read-only.
    ##
    ## @end deftp
    ActiveSetMethod       = 'Random';

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} ActiveSetSize
    ##
    ## Size of the active set
    ##
    ## A positive integer scalar.  This property is read-only.
    ##
    ## @end deftp
    ActiveSetSize         = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} PredictorLocation
    ##
    ## Means the predictors were centred by
    ##
    ## A @math{1xP} numeric vector, or empty.  This property is read-only.
    ##
    ## @end deftp
    PredictorLocation     = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} PredictorScale
    ##
    ## Standard deviations the predictors were scaled by
    ##
    ## A @math{1xP} numeric vector, or empty.  This property is read-only.
    ##
    ## @end deftp
    PredictorScale        = [];

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGP} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A character vector, or the text of the function handle that was
    ## supplied.  Assigning to it accepts either.
    ##
    ## @end deftp
    ResponseTransform     = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The callable behind ResponseTransform.
    RTfun                 = @(y) y;

  endproperties

  methods (Access = public)

    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = ...
                     parseResponseTransform (val, 'CompactRegressionGP');
    endfunction

    ## Class constructor
    function this = CompactRegressionGP (Mdl)

      if (nargin < 1)
        error ("CompactRegressionGP: too few input arguments.");
      endif
      if (! isa (Mdl, 'RegressionGP'))
        error (strcat ("CompactRegressionGP: MDL must be a RegressionGP", ...
                       " object."));
      endif

      this.PredictorNames = Mdl.PredictorNames;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;
      this.ResponseName = Mdl.ResponseName;
      this.CategoricalPredictors = Mdl.CategoricalPredictors;
      this.FitMethod = Mdl.FitMethod;
      this.BasisFunction = Mdl.BasisFunction;
      this.Beta = Mdl.Beta;
      this.Sigma = Mdl.Sigma;
      this.KernelFunction = Mdl.KernelFunction;
      this.KernelInformation = Mdl.KernelInformation;
      this.PredictMethod = Mdl.PredictMethod;
      this.Alpha = Mdl.Alpha;
      this.ActiveSetVectors = Mdl.ActiveSetVectors;
      this.ActiveSetMethod = Mdl.ActiveSetMethod;
      this.ActiveSetSize = Mdl.ActiveSetSize;
      this.PredictorLocation = Mdl.PredictorLocation;
      this.PredictorScale = Mdl.PredictorScale;
      this.ResponseTransform = Mdl.ResponseTransform;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionGP} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactRegressionGP} {[@var{yFit}, @var{ySD}, @var{yInt}] =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactRegressionGP} {[@dots{}] =} predict (@dots{}, @qcode{'Alpha'}, @var{alpha})
    ##
    ## Predict the response for new data with a compact Gaussian process model.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{XC})} returns the predicted
    ## response of the @qcode{CompactRegressionGP} model @var{obj} at the
    ## points in @var{XC}, and the further outputs are the standard deviation
    ## of each predicted response and the prediction intervals, exactly as the
    ## full model returns them.
    ##
    ## @end deftypefn
    function [yFit, ySD, yInt] = predict (this, XC, varargin)

      if (nargin < 2)
        error ("CompactRegressionGP.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("CompactRegressionGP.predict: XC is empty.");
      endif
      if (columns (XC) != columns (this.ActiveSetVectors))
        error (strcat ("CompactRegressionGP.predict: XC must have the same", ...
                       " number of predictors as the trained model."));
      endif

      CIAlpha = 0.05;
      while (numel (varargin) > 0)
        if (numel (varargin) < 2)
          error (strcat ("CompactRegressionGP.predict: optional arguments", ...
                         " must be given in Name-Value pairs."));
        endif
        switch (lower (varargin{1}))
          case 'alpha'
            CIAlpha = varargin{2};
            if (! (isnumeric (CIAlpha) && isscalar (CIAlpha) && ...
                   CIAlpha >= 0 && CIAlpha <= 1))
              error (strcat ("CompactRegressionGP.predict: 'Alpha' must", ...
                             " be a scalar between 0 and 1."));
            endif
          otherwise
            error (strcat ("CompactRegressionGP.predict: invalid NAME in", ...
                           " optional pairs of arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      M = struct ('X', this.ActiveSetVectors, 'Alpha', this.Alpha, ...
                  'KernelFunction', this.KernelFunction, ...
                  'Theta', this.KernelInformation.KernelParameters, ...
                  'BasisFunction', this.BasisFunction, 'Beta', this.Beta, ...
                  'Sigma', this.Sigma, 'Location', this.PredictorLocation, ...
                  'Scale', this.PredictorScale, 'CIAlpha', CIAlpha);

      if (nargout < 2)
        yFit = this.RTfun (gpPredict (XC, M));
      elseif (nargout < 3)
        [yFit, ySD] = gpPredict (XC, M);
        yFit = this.RTfun (yFit);
      else
        [yFit, ySD, yInt] = gpPredict (XC, M);
        yFit = this.RTfun (yFit);
        yInt = this.RTfun (yInt);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionGP} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactRegressionGP} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the regression loss of a compact Gaussian process model.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the mean
    ## squared error of the model @var{obj} on the data @var{X} and @var{Y},
    ## and accepts the same @qcode{'LossFun'} and @qcode{'Weights'} pairs the
    ## full model accepts.
    ##
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactRegressionGP.loss: too few input arguments.");
      endif
      if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
        error ("CompactRegressionGP.loss: invalid values in X.");
      endif
      if (! (isnumeric (Y) && isreal (Y) && isvector (Y)))
        error ("CompactRegressionGP.loss: invalid values in Y.");
      endif
      Y = Y(:);
      if (rows (X) != rows (Y))
        error (strcat ("CompactRegressionGP.loss: number of rows in X and", ...
                       " Y must be equal."));
      endif

      LossFun = 'mse';
      Weights = ones (rows (X), 1);
      Epsilon = 0;
      while (numel (varargin) > 0)
        if (numel (varargin) < 2)
          error (strcat ("CompactRegressionGP.loss: optional arguments", ...
                         " must be given in Name-Value pairs."));
        endif
        switch (lower (varargin{1}))
          case 'lossfun'
            LossFun = varargin{2};
            if (! (ischar (LossFun) || is_function_handle (LossFun)))
              error (strcat ("CompactRegressionGP.loss: 'LossFun' must be", ...
                             " a character vector or a function handle."));
            endif
            if (ischar (LossFun) && ...
                ! any (strcmpi (LossFun, {'mse', 'mae', ...
                                          'epsiloninsensitive'})))
              error (strcat ("CompactRegressionGP.loss: unsupported", ...
                             " 'LossFun' value."));
            endif
          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isvector (Weights) && ...
                   numel (Weights) == rows (X) && all (Weights >= 0)))
              error (strcat ("CompactRegressionGP.loss: 'Weights' must be", ...
                             " a vector of non-negative values with one", ...
                             " element per observation."));
            endif
            Weights = Weights(:);
          case 'epsilon'
            Epsilon = varargin{2};
          otherwise
            error (strcat ("CompactRegressionGP.loss: invalid NAME in", ...
                           " optional pairs of arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      yFit = this.predict (X);
      if (is_function_handle (LossFun))
        L = LossFun (Y, yFit);
        return;
      endif

      switch (lower (LossFun))
        case 'mse'
          L = sum (Weights .* (Y - yFit) .^ 2) / sum (Weights);
        case 'mae'
          L = sum (Weights .* abs (Y - yFit)) / sum (Weights);
        case 'epsiloninsensitive'
          e = max (0, abs (Y - yFit) - Epsilon);
          L = sum (Weights .* e) / sum (Weights);
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionGP} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a compact Gaussian process model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves the model @var{obj}
    ## into @var{filename} in a form @code{loadmodel} can read back.
    ##
    ## @end deftypefn
    function savemodel (obj, fname)

      classdef_name = 'CompactRegressionGP';
      PredictorNames = obj.PredictorNames;
      ExpandedPredictorNames = obj.ExpandedPredictorNames;
      ResponseName = obj.ResponseName;
      CategoricalPredictors = obj.CategoricalPredictors;
      FitMethod = obj.FitMethod;
      BasisFunction = obj.BasisFunction;
      Beta = obj.Beta;
      Sigma = obj.Sigma;
      KernelFunction = obj.KernelFunction;
      KernelInformation = obj.KernelInformation;
      PredictMethod = obj.PredictMethod;
      Alpha = obj.Alpha;
      ActiveSetVectors = obj.ActiveSetVectors;
      ActiveSetMethod = obj.ActiveSetMethod;
      ActiveSetSize = obj.ActiveSetSize;
      PredictorLocation = obj.PredictorLocation;
      PredictorScale = obj.PredictorScale;
      ResponseTransform = obj.ResponseTransform;

      save ('-binary', fname, 'classdef_name', 'PredictorNames', ...
            'ExpandedPredictorNames', 'ResponseName', ...
            'CategoricalPredictors', 'FitMethod', 'BasisFunction', 'Beta', ...
            'Sigma', 'KernelFunction', ...
            'KernelInformation', 'PredictMethod', 'Alpha', ...
            'ActiveSetVectors', 'ActiveSetMethod', 'ActiveSetSize', ...
            'PredictorLocation', 'PredictorScale', 'ResponseTransform');

    endfunction

  endmethods

  methods (Access = public, Hidden)

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        printf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      printf ("\n  CompactRegressionGP\n\n");
      printf ("%25s: '%s'\n", 'ResponseName', this.ResponseName);
      printf ("%25s: %d\n", 'NumPredictors', ...
              columns (this.ActiveSetVectors));
      if (is_function_handle (this.KernelFunction))
        printf ("%25s: '%s'\n", 'KernelFunction', ...
                func2str (this.KernelFunction));
      else
        printf ("%25s: '%s'\n", 'KernelFunction', this.KernelFunction);
      endif
      printf ("%25s: '%s'\n", 'PredictMethod', this.PredictMethod);
      printf ("%25s: %g\n", 'Sigma', this.Sigma);
      printf ("\n");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      ## The compact model is rebuilt field by field: it has no training data
      ## to construct itself from, so an empty shell is filled directly.
      mdl = CompactRegressionGP.empty_model ();
      fields = fieldnames (data);
      for k = 1:numel (fields)
        mdl.(fields{k}) = data.(fields{k});
      endfor

    endfunction

    ## An unfitted shell, which only load_model needs.
    function mdl = empty_model ()
      Mdl = RegressionGP ([0; 1], [0; 1], 'FitMethod', 'none');
      mdl = CompactRegressionGP (Mdl);
    endfunction

  endmethods

endclassdef

%!demo
%! ## A compact model predicts what the full model predicts, and carries none
%! ## of the training data.
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.05 * cos (9*x);
%! Mdl = fitrgp (x, y);
%! CMdl = compact (Mdl)
%! xq = [0.15; 0.55; 0.85];
%! [yq, ysd] = predict (CMdl, xq)

%!test
%! ## The compact model carries the fitted surface and nothing that describes
%! ## the training data
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! Mdl = RegressionGP (x, y);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), 'CompactRegressionGP');
%! assert_equal (CMdl.Beta, Mdl.Beta);
%! assert_equal (CMdl.Sigma, Mdl.Sigma);
%! assert_equal (CMdl.KernelFunction, Mdl.KernelFunction);
%! assert_equal (CMdl.KernelInformation, Mdl.KernelInformation);
%! assert_equal (CMdl.Alpha, Mdl.Alpha);
%! assert_equal (CMdl.ActiveSetVectors, Mdl.ActiveSetVectors);
%! assert_equal (CMdl.ResponseName, Mdl.ResponseName);
%! assert_equal (CMdl.PredictorNames, Mdl.PredictorNames);

%!test
%! ## It predicts exactly what the full model predicts, standard deviation
%! ## and interval included, because it predicts through the same code
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! Mdl = RegressionGP (x, y);
%! CMdl = compact (Mdl);
%! xq = [0.05; 0.33; 0.5; 0.77; 0.95];
%! [y1, s1, i1] = predict (Mdl, xq);
%! [y2, s2, i2] = predict (CMdl, xq);
%! assert_equal (y2, y1, 1e-14);
%! assert_equal (s2, s1, 1e-14);
%! assert_equal (i2, i1, 1e-14);

%!test
%! ## The R2024a values, reached through the compact model
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! CMdl = compact (RegressionGP (x, y));
%! xq = [0.05; 0.33; 0.5; 0.77; 0.95];
%! assert_equal (predict (CMdl, xq), ...
%!               [0.404401855406521; 0.809355242891053; ...
%!                -0.093708552144171; -0.928420803239529; ...
%!                -0.217104024154071], 1e-7);

%!test
%! ## The interval level is settable here as it is on the full model
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! CMdl = compact (RegressionGP (x, y));
%! xq = [0.2; 0.6];
%! [yp, ysd, yint] = predict (CMdl, xq);
%! assert_equal (yint(:,2) - yp, norminv (0.975) * ysd, 1e-12);
%! [~, ~, yint90] = predict (CMdl, xq, 'Alpha', 0.10);
%! assert (all (yint90(:,2) - yint90(:,1) < yint(:,2) - yint(:,1)));

%!test
%! ## Standardization is carried over, so the compact model transforms new
%! ## data the way the full model did
%! X = [linspace(0, 10, 20)', linspace(-5, 5, 20)'];
%! y = 0.3 * X(:,1) - 0.2 * X(:,2);
%! Mdl = RegressionGP (X, y, 'Standardize', true);
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.PredictorLocation, Mdl.PredictorLocation);
%! assert_equal (CMdl.PredictorScale, Mdl.PredictorScale);
%! assert_equal (predict (CMdl, X), predict (Mdl, X), 1e-14);

%!test
%! ## loss on the compact model agrees with loss on the full one
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! Mdl = RegressionGP (x, y);
%! CMdl = compact (Mdl);
%! assert_equal (loss (CMdl, x, y), loss (Mdl, x, y), 1e-14);
%! assert_equal (loss (CMdl, x, y, 'LossFun', 'mae'), ...
%!               loss (Mdl, x, y, 'LossFun', 'mae'), 1e-14);
%! w = linspace (1, 2, 15)';
%! assert_equal (loss (CMdl, x, y, 'Weights', w), ...
%!               loss (Mdl, x, y, 'Weights', w), 1e-14);

%!test
%! ## A compact model saved and loaded predicts what it predicted before
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! CMdl = compact (RegressionGP (x, y));
%! fname = tempname ();
%! savemodel (CMdl, fname);
%! CMdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (CMdl2), 'CompactRegressionGP');
%! assert_equal (CMdl2.Beta, CMdl.Beta);
%! assert_equal (CMdl2.Sigma, CMdl.Sigma);
%! assert_equal (CMdl2.Alpha, CMdl.Alpha);
%! assert_equal (CMdl2.ActiveSetVectors, CMdl.ActiveSetVectors);
%! assert_equal (predict (CMdl2, x), predict (CMdl, x), 1e-14);

%!test
%! ## A response transform survives compacting
%! x = linspace (0, 1, 12)';
%! y = cos (3*x);
%! Mdl = RegressionGP (x, y, 'ResponseTransform', 'exp');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.ResponseTransform, 'exp');
%! assert_equal (predict (CMdl, x), predict (Mdl, x), 1e-14);

## Test input validation for the constructor
%!error<CompactRegressionGP: too few input arguments.> CompactRegressionGP ()
%!error<CompactRegressionGP: MDL must be a RegressionGP object.> ...
%! CompactRegressionGP (5)

## Test input validation for the predict method
%!error<CompactRegressionGP.predict: too few input arguments.> ...
%! predict (compact (RegressionGP (ones (5, 2), ones (5, 1))))
%!error<CompactRegressionGP.predict: XC is empty.> ...
%! predict (compact (RegressionGP (ones (5, 2), ones (5, 1))), [])
%!error<CompactRegressionGP.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (compact (RegressionGP (ones (5, 2), ones (5, 1))), ones (3, 3))
%!error<CompactRegressionGP.predict: 'Alpha' must be a scalar between 0 and 1.> ...
%! predict (compact (RegressionGP (ones (5, 2), ones (5, 1))), ...
%!          ones (3, 2), 'Alpha', 2)
%!error<CompactRegressionGP.predict: invalid NAME in optional pairs of arguments.> ...
%! predict (compact (RegressionGP (ones (5, 2), ones (5, 1))), ...
%!          ones (3, 2), 'bogus', 1)

## Test input validation for the loss method
%!error<CompactRegressionGP.loss: too few input arguments.> ...
%! loss (compact (RegressionGP (ones (5, 2), ones (5, 1))), ones (3, 2))
%!error<CompactRegressionGP.loss: invalid values in X.> ...
%! loss (compact (RegressionGP (ones (5, 2), ones (5, 1))), 'a', ones (3, 1))
%!error<CompactRegressionGP.loss: number of rows in X and Y must be equal.> ...
%! loss (compact (RegressionGP (ones (5, 2), ones (5, 1))), ...
%!       ones (3, 2), ones (2, 1))
%!error<CompactRegressionGP.loss: unsupported 'LossFun' value.> ...
%! loss (compact (RegressionGP (ones (5, 2), ones (5, 1))), ones (3, 2), ...
%!       ones (3, 1), 'LossFun', 'bogus')
