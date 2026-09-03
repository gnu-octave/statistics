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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftp {statistics} {} ReconstructionICA
##
## Reconstruction independent component analysis (RICA) feature-extraction model.
##
## A @qcode{ReconstructionICA} object stores the transformation learned by
## @code{rica} for extracting features from data.  Create one with @code{rica};
## apply it to data with the @code{transform} method.
##
## @seealso{rica, sparsefilt}
## @end deftp

classdef ReconstructionICA

  properties (SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ReconstructionICA} {property} ModelParameters
    ##
    ## Options the fit used
    ##
    ## A scalar structure holding the options the fit ran with:
    ## @qcode{IterationLimit}, @qcode{Lambda}, @qcode{Standardize},
    ## @qcode{ContrastFcn}, @qcode{InitialTransformWeights},
    ## @qcode{GradientTolerance}, @qcode{StepTolerance}, @qcode{Solver} and
    ## @qcode{NonGaussianityIndicator}.
    ## This property is read-only.
    ##
    ## @end deftp
    ModelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {ReconstructionICA} {property} NumPredictors
    ##
    ## Number of input predictors
    ##
    ## A positive integer @math{P}, the number of columns of the training
    ## data.  This property is read-only.
    ##
    ## @end deftp
    NumPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ReconstructionICA} {property} NumLearnedFeatures
    ##
    ## Number of learned features
    ##
    ## A positive integer @math{Q}, the number of features the learned
    ## transformation produces.  This property is read-only.
    ##
    ## @end deftp
    NumLearnedFeatures = [];

    ## -*- texinfo -*-
    ## @deftp {ReconstructionICA} {property} Mu
    ##
    ## Predictor means used when standardizing
    ##
    ## A column vector with one entry per predictor, the mean of each
    ## column of the training data.  It is empty unless @qcode{'Standardize'}
    ## was true.  This property is read-only.
    ##
    ## @end deftp
    Mu = [];

    ## -*- texinfo -*-
    ## @deftp {ReconstructionICA} {property} Sigma
    ##
    ## Predictor standard deviations used when standardizing
    ##
    ## A column vector with one entry per predictor, the standard deviation
    ## of each column of the training data.  It is empty unless
    ## @qcode{'Standardize'} was true.  This property is read-only.
    ##
    ## @end deftp
    Sigma = [];

    ## -*- texinfo -*-
    ## @deftp {ReconstructionICA} {property} FitInfo
    ##
    ## History of the fit
    ##
    ## A scalar structure with the fields @qcode{Iteration} and
    ## @qcode{Objective}, both column vectors of the same length.
    ## @qcode{Iteration} counts from zero and @qcode{Objective(1)} is the
    ## objective at the starting weights, so the last entry of each is the
    ## solution the fit returned.  This property is read-only.
    ##
    ## The trajectory is this implementation's own.  The default
    ## @qcode{'quasinewton'} solver minimises through Octave's @code{fminunc},
    ## and @qcode{'Solver', 'lbfgs'} selects the limited-memory BFGS solver
    ## MATLAB uses.  Either way the steps taken from the same starting weights
    ## differ from MATLAB's, so the length of the history and the iteration
    ## counts differ, and on an objective this far from convex the optimum
    ## reached need not be MATLAB's either.
    ##
    ## @end deftp
    FitInfo = [];

    ## -*- texinfo -*-
    ## @deftp {ReconstructionICA} {property} TransformWeights
    ##
    ## Learned feature transformation weights
    ##
    ## A @math{P}-by-@math{Q} matrix of learned weights, its columns of unit
    ## length.  The @code{transform} method applies it to data.  This property
    ## is read-only.
    ##
    ## @end deftp
    TransformWeights = [];

    ## -*- texinfo -*-
    ## @deftp {ReconstructionICA} {property} InitialTransformWeights
    ##
    ## Starting feature transformation weights
    ##
    ## A @math{P}-by-@math{Q} matrix, the weights the fit started from.  It
    ## is the matrix given as @qcode{'InitialTransformWeights'} when one was
    ## given, and the random start the fit drew otherwise.  This property is
    ## read-only.
    ##
    ## @end deftp
    InitialTransformWeights = [];

    ## -*- texinfo -*-
    ## @deftp {ReconstructionICA} {property} NonGaussianityIndicator
    ##
    ## Non-Gaussianity of each learned feature
    ##
    ## A @math{Q}-by-1 vector of @math{+1} and @math{-1}, one per learned
    ## feature: @math{+1} where the feature is taken to be super-Gaussian and
    ## @math{-1} where it is taken to be sub-Gaussian.  The entry sets the
    ## sign its feature's contrast term carries in the objective, so the fit
    ## seeks a sparse feature where the entry is @math{+1} and a spread one
    ## where it is @math{-1}.  The default is all @math{+1}.  This property is
    ## read-only.
    ##
    ## @end deftp
    NonGaussianityIndicator = [];

  endproperties

  methods

    ## -*- texinfo -*-
    ## @deftypefn {statistics} {@var{Mdl} =} ReconstructionICA (@var{X}, @var{Q}, @dots{})
    ## Fit a reconstruction ICA model.  This constructor is invoked by
    ## @code{rica}; see @code{help rica} for the arguments.
    ## @end deftypefn
    function this = ReconstructionICA (X, Q, varargin)

      if (nargin == 0)
        return;
      endif

      [n, p] = size (X);

      ## Options
      opts = struct ("IterationLimit", 1000, "Lambda", 1, ...
                     "Standardize", false, "ContrastFcn", "logcosh", ...
                     "InitialTransformWeights", [], ...
                     "GradientTolerance", 1e-6, "StepTolerance", 1e-6, ...
                     "Solver", "quasinewton", ...
                     "NonGaussianityIndicator", []);
      for k = 1:2:numel (varargin)
        name = varargin{k};
        val = varargin{k+1};
        switch (lower (name))
          case 'iterationlimit'
            opts.IterationLimit = val;
          case 'lambda'
            opts.Lambda = val;
          case 'standardize'
            opts.Standardize = val;
          case 'contrastfcn'
            opts.ContrastFcn = lower (val);
          case 'initialtransformweights'
            opts.InitialTransformWeights = val;
          case 'gradienttolerance'
            opts.GradientTolerance = val;
          case 'steptolerance'
            opts.StepTolerance = val;
          case 'nongaussianityindicator'
            opts.NonGaussianityIndicator = val;
          case 'solver'
            if (! ischar (val) || ...
                ! any (strcmpi (val, {'quasinewton', 'lbfgs'})))
              error ("rica: 'Solver' must be 'quasinewton' or 'lbfgs'.");
            endif
            opts.Solver = lower (val);
          otherwise
            error ("rica: unknown parameter name '%s'.", name);
        endswitch
      endfor

      if (! any (strcmp (opts.ContrastFcn, {'logcosh', 'exp', 'sqrt'})))
        error ("rica: 'ContrastFcn' must be 'logcosh', 'exp', or 'sqrt'.");
      endif

      ## One sign per learned feature, all of them super-Gaussian by default
      if (isempty (opts.NonGaussianityIndicator))
        opts.NonGaussianityIndicator = ones (Q, 1);
      else
        sigma = opts.NonGaussianityIndicator;
        if (! (isnumeric (sigma) && isreal (sigma) && isvector (sigma)
               && numel (sigma) == Q && all (abs (sigma(:)) == 1)))
          error (strcat ("rica: 'NonGaussianityIndicator' must be a real", ...
                         " vector of Q elements, each +1 or -1."));
        endif
        opts.NonGaussianityIndicator = sigma(:);
      endif

      ## Standardize (center and scale) the data if requested
      if (opts.Standardize)
        ## MU and SIGMA are held as columns, one entry per predictor, as
        ## MATLAB holds them, and transposed where the data is scaled.
        this.Mu = mean (X)';
        this.Sigma = std (X)';
        sig = this.Sigma';
        sig(sig == 0) = 1;
        X = (X - this.Mu') ./ sig;
      endif

      ## Initial weights
      if (isempty (opts.InitialTransformWeights))
        W0 = randn (p, Q);
      else
        W0 = opts.InitialTransformWeights;
      endif
      this.InitialTransformWeights = W0;

      ## Minimize the RICA objective
      lambda = opts.Lambda;
      sigma = opts.NonGaussianityIndicator;
      ofun = @(wv) __rica_objective__ (wv, X, p, Q, lambda, ...
                                       opts.ContrastFcn, sigma);
      if (strcmpi (opts.Solver, 'lbfgs'))
        ## LossTolerance is switched off because MATLAB has no such option
        ## here: the fit stops on the gradient or the step, never on the
        ## objective's own value.
        ## The history is sized to the problem rather than left at the
        ## engine's default ten.  A transform carries only p * Q parameters,
        ## and once the stored pairs reach that count the recursion is full
        ## BFGS, which converges further here and in fewer iterations.  The
        ## cap bounds the memory when Q is large.
        lbopts = struct ("IterationLimit", opts.IterationLimit, ...
                         "GradientTolerance", opts.GradientTolerance, ...
                         "StepTolerance", opts.StepTolerance, ...
                         "LossTolerance", -Inf, ...
                         "HistorySize", min (max (10, numel (W0)), 50));
        [wv, lbinfo] = __lbfgs__ (ofun, W0(:), lbopts);
        fval = lbinfo.Fval;
        ## The engine records from the first step, where MATLAB's trajectory
        ## starts at the initial weights.
        f0 = ofun (W0(:));
        obj = [f0; lbinfo.History.Fval];
      else
        objective_history ("reset");
        fmopts = optimset ("GradObj", "on", "MaxIter", opts.IterationLimit, ...
                           "TolFun", 1e-10, "TolX", 1e-10, "Display", "off", ...
                           "OutputFcn", @objective_history_fcn);
        [wv, fval, ~, output] = fminunc (ofun, W0(:), fmopts);
        obj = objective_history ("get");
      endif

      W = reshape (wv, p, Q);
      W = W ./ sqrt (sum (W .^ 2, 1));   ## unit-length columns
      this.TransformWeights = W;

      this.NumPredictors = p;
      this.NumLearnedFeatures = Q;
      this.NonGaussianityIndicator = opts.NonGaussianityIndicator;
      ## MATLAB reports the whole trajectory: a column with the objective at
      ## the starting weights first, one entry per iteration after it, and
      ## Iteration the matching 0-based index.  fminunc does not call its
      ## OutputFcn for the step it returns on, so the final value is appended.
      if (isempty (obj) || obj(end) != fval)
        obj(end+1, 1) = fval;
      endif
      this.FitInfo = struct ("Iteration", (0:numel (obj) - 1)', ...
                             "Objective", obj);
      this.ModelParameters = opts;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {statistics} {@var{Z} =} transform (@var{Mdl}, @var{X})
    ## Transform data @var{X} into the learned feature space, returning the
    ## @math{N * Q} matrix @var{Z} of features.
    ## @end deftypefn
    function Z = transform (this, X)
      if (nargin != 2)
        print_usage ();
      endif
      if (! isempty (this.Mu))
        sig = this.Sigma';
        sig(sig == 0) = 1;
        X = (X - this.Mu') ./ sig;
      endif
      Z = X * this.TransformWeights;
    endfunction

  endmethods

endclassdef

## RICA objective and gradient with respect to the raw (un-normalized) weights.
## SQRT_EPS smooths abs (z) at the origin for the 'sqrt' contrast; the value is
## MATLAB's own, measured on R2024a against a fixture carrying a zero, which is
## the only place it is visible.
function [f, g] = __rica_objective__ (wv, X, p, Q, lambda, contrast, sigma)

  SQRT_EPS = 1e-8;

  W = reshape (wv, p, Q);
  nrm = sqrt (sum (W .^ 2, 1));
  Wn = W ./ nrm;
  Z = X * Wn;
  E = X * (Wn * Wn') - X;

  switch (contrast)
    case 'logcosh'
      C = 0.5 * log (cosh (2 * Z));
    case 'exp'
      C = -exp (-Z .^ 2 / 2);
    case 'sqrt'
      C = sqrt (Z .^ 2 + SQRT_EPS);
  endswitch
  f = lambda * sum (E(:) .^ 2) + sum (C * sigma);

  if (nargout > 1)
    switch (contrast)
      case 'logcosh'
        dC = tanh (2 * Z);
      case 'exp'
        dC = Z .* exp (-Z .^ 2 / 2);
      case 'sqrt'
        dC = Z ./ sqrt (Z .^ 2 + SQRT_EPS);
    endswitch
    gWn = lambda * 2 * (X' * E * Wn + E' * X * Wn) + X' * (dC .* sigma');
    g = zeros (p, Q);
    for j = 1:Q
      wnj = Wn(:,j);
      g(:,j) = (gWn(:,j) - wnj * (wnj' * gWn(:,j))) / nrm(j);
    endfor
    g = g(:);
  endif

endfunction

%!shared X, W0
%! X = reshape (mod ((1:60)*7, 13), 12, 5) - 6;
%! W0 = reshape (mod ((1:15)*3, 7), 5, 3) - 3;

%!test
%! ## Construct the object directly and check its properties.
%! Mdl = ReconstructionICA (X, 3, "InitialTransformWeights", W0, ...
%!                          "Lambda", 1, "IterationLimit", 500);
%! assert_equal (isa (Mdl, "ReconstructionICA"), true);
%! assert_equal (Mdl.NumPredictors, 5);
%! assert_equal (Mdl.NumLearnedFeatures, 3);
%! assert_equal (size (Mdl.TransformWeights), [5, 3]);
%! assert_equal (Mdl.InitialTransformWeights, W0);
%! assert_equal (isfield (Mdl.FitInfo, "Objective"), true);

%!test
%! ## The transform method projects onto the (unit-column) weights.
%! Mdl = ReconstructionICA (X, 3, "InitialTransformWeights", W0, ...
%!                          "IterationLimit", 500);
%! Z = transform (Mdl, X);
%! assert_equal (Z, X * Mdl.TransformWeights, 1e-12);
%! assert_equal (sqrt (sum (Mdl.TransformWeights .^ 2, 1)), [1, 1, 1], 1e-10);

%!test
%! ## The transform method standardizes new data when the model does.
%! Mdl = ReconstructionICA (X, 2, "Standardize", true, "IterationLimit", 100, ...
%!                          "InitialTransformWeights", W0(:,1:2));
%! Z = transform (Mdl, X);
%! Xs = (X - Mdl.Mu') ./ Mdl.Sigma';
%! assert_equal (Z, Xs * Mdl.TransformWeights, 1e-12);
%!test
%! ## Mu and Sigma are columns, one entry per predictor, as MATLAB reports them
%! Mdl = ReconstructionICA (X, 2, "Standardize", true, "IterationLimit", 100, ...
%!                          "InitialTransformWeights", W0(:,1:2));
%! assert_equal (size (Mdl.Mu), [5, 1]);
%! assert_equal (size (Mdl.Sigma), [5, 1]);
%! assert_equal (Mdl.Mu, mean (X)', 1e-12);
%! assert_equal (Mdl.Sigma, std (X)', 1e-12);

%!test
%! ## The default constructor returns an empty model.
%! Mdl = ReconstructionICA ();
%! assert_equal (isempty (Mdl.TransformWeights), true);

%!test
%! ## NonGaussianityIndicator defaults to one per feature, all +1.
%! Mdl = ReconstructionICA (X, 3, "InitialTransformWeights", W0, ...
%!                          "IterationLimit", 100);
%! assert_equal (Mdl.NonGaussianityIndicator, [1; 1; 1]);

%!test
%! ## A given indicator is held as a column, whatever shape it arrives in.
%! Mdl = ReconstructionICA (X, 3, "InitialTransformWeights", W0, ...
%!                          "IterationLimit", 100, ...
%!                          "NonGaussianityIndicator", [1, -1, 1]);
%! assert_equal (Mdl.NonGaussianityIndicator, [1; -1; 1]);

%!test
%! ## The indicator signs each feature's contrast term, so it moves the fit.
%! args = {"InitialTransformWeights", W0, "IterationLimit", 200};
%! Mdl = ReconstructionICA (X, 3, args{:});
%! Neg = ReconstructionICA (X, 3, args{:}, ...
%!                          "NonGaussianityIndicator", [-1; 1; 1]);
%! assert_equal (isequal (Mdl.TransformWeights, Neg.TransformWeights), false);

%!test
%! ## All +1 is the default, so it must reproduce the default fit exactly.
%! args = {"InitialTransformWeights", W0, "IterationLimit", 200};
%! Mdl = ReconstructionICA (X, 3, args{:});
%! Pos = ReconstructionICA (X, 3, args{:}, ...
%!                          "NonGaussianityIndicator", [1; 1; 1]);
%! assert_equal (Pos.TransformWeights, Mdl.TransformWeights);

%!test
%! ## Objective against MATLAB R2024a: 801.1816174574 and 197.5588279632.
%! t = (1:80)';
%! Xr = double ([mod(t*7,11)-5, mod(t*13,17)-8, mod(t*5,7)-3]);
%! args = {"InitialTransformWeights", [1, 0; 0, 1; 1, 1], "Lambda", 1, ...
%!         "ContrastFcn", "logcosh", "IterationLimit", 1000, ...
%!         "Solver", "lbfgs", "GradientTolerance", 1e-10, ...
%!         "StepTolerance", 1e-10};
%! Mdl = ReconstructionICA (Xr, 2, args{:});
%! Neg = ReconstructionICA (Xr, 2, args{:}, ...
%!                          "NonGaussianityIndicator", [-1; 1]);
%! assert_equal (Mdl.FitInfo.Objective(end), 801.1816174574, 1e-6);
%! assert_equal (Neg.FitInfo.Objective(end), 197.5588279632, 1e-6);

%!error<rica: 'NonGaussianityIndicator' must be a real vector of Q elements, each \+1 or -1.> ...
%! ReconstructionICA (ones (6, 5), 2, "NonGaussianityIndicator", [1; 1; 1])
%!error<rica: 'NonGaussianityIndicator' must be a real vector of Q elements, each \+1 or -1.> ...
%! ReconstructionICA (ones (6, 5), 2, "NonGaussianityIndicator", [1; 0])
%!error<rica: 'NonGaussianityIndicator' must be a real vector of Q elements, each \+1 or -1.> ...
%! ReconstructionICA (ones (6, 5), 2, "NonGaussianityIndicator", {1, -1})

%!error<rica: 'ContrastFcn' must be 'logcosh', 'exp', or 'sqrt'.> ...
%! ReconstructionICA (ones (6, 5), 2, "ContrastFcn", "bogus")
%!error<rica: unknown parameter name 'bogus'.> ...
%! ReconstructionICA (ones (6, 5), 2, "bogus", 1)
