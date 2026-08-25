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
## @deftp {statistics} {} SparseFiltering
##
## Sparse filtering feature-extraction model.
##
## A @qcode{SparseFiltering} object stores the transformation learned by
## @code{sparsefilt} for extracting features from data.  Create one with
## @code{sparsefilt}; apply it to data with the @code{transform} method.
##
## An object of this class has the following properties:
##
## @itemize
## @item @qcode{TransformWeights}: the @math{P * Q} matrix of learned
## transformation weights.
## @item @qcode{Mu}, @qcode{Sigma}: the per-predictor mean and standard
## deviation used when @qcode{'Standardize'} is @qcode{true} (empty otherwise),
## each a column vector with one entry per predictor.
## @item @qcode{FitInfo}: a structure with the @qcode{Iteration} indices of
## the fit and the @qcode{Objective} value at each of them, both column
## vectors of the same length.  @qcode{Iteration} counts from zero and
## @qcode{Objective(1)} is the objective at the starting weights, so the last
## entry of each is the solution the fit returned.
##
## The values along that trajectory are this implementation's own.  Under the
## default @qcode{'quasinewton'} solver the minimisation runs through Octave's
## @code{fminunc}; @qcode{'Solver', 'lbfgs'} selects the limited-memory BFGS
## solver MATLAB uses.  Either way the steps taken from the same starting
## weights differ from MATLAB's, so the length of the history and the
## iteration counts differ, and on an objective this far from convex the
## optimum reached need not be MATLAB's either.
## @item @qcode{ModelParameters}: a structure of the options used for the fit.
## @item @qcode{NumPredictors}, @qcode{NumLearnedFeatures}: the input dimension
## @math{P} and the number of learned features @math{Q}.
## @item @qcode{InitialTransformWeights}: the starting weights used by the fit.
## @end itemize
##
## @seealso{sparsefilt, rica}
## @end deftp

classdef SparseFiltering

  properties (SetAccess = protected)
    ModelParameters = [];
    NumPredictors = [];
    NumLearnedFeatures = [];
    Mu = [];
    Sigma = [];
    FitInfo = [];
    TransformWeights = [];
    InitialTransformWeights = [];
  endproperties

  methods

    ## -*- texinfo -*-
    ## @deftypefn {statistics} {@var{Mdl} =} SparseFiltering (@var{X}, @var{Q}, @dots{})
    ## Fit a sparse filtering model.  This constructor is invoked by
    ## @code{sparsefilt}; see @code{help sparsefilt} for the arguments.
    ## @end deftypefn
    function this = SparseFiltering (X, Q, varargin)

      if (nargin == 0)
        return;
      endif

      [n, p] = size (X);

      opts = struct ("IterationLimit", 1000, "Lambda", 1, ...
                     "Standardize", false, ...
                     "InitialTransformWeights", [], ...
                     "GradientTolerance", 1e-6, "StepTolerance", 1e-6, ...
                     "Solver", "quasinewton");
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
          case 'initialtransformweights'
            opts.InitialTransformWeights = val;
          case 'gradienttolerance'
            opts.GradientTolerance = val;
          case 'steptolerance'
            opts.StepTolerance = val;
          case 'solver'
            if (! ischar (val) || ...
                ! any (strcmpi (val, {'quasinewton', 'lbfgs'})))
              error ("sparsefilt: 'Solver' must be 'quasinewton' or 'lbfgs'.");
            endif
            opts.Solver = lower (val);
          otherwise
            error ("sparsefilt: unknown parameter name '%s'.", name);
        endswitch
      endfor

      if (opts.Standardize)
        ## MU and SIGMA are held as columns, one entry per predictor, as
        ## MATLAB holds them, and transposed where the data is scaled.
        this.Mu = mean (X)';
        this.Sigma = std (X)';
        sig = this.Sigma';
        sig(sig == 0) = 1;
        X = (X - this.Mu') ./ sig;
      endif

      if (isempty (opts.InitialTransformWeights))
        W0 = randn (p, Q);
      else
        W0 = opts.InitialTransformWeights;
      endif
      this.InitialTransformWeights = W0;

      lambda = opts.Lambda;
      ofun = @(wv) __sparsefilt_objective__ (wv, X, p, Q, lambda);
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

      this.TransformWeights = reshape (wv, p, Q);
      this.NumPredictors = p;
      this.NumLearnedFeatures = Q;
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
    ## @math{N * Q} matrix @var{Z} of sparse features.
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
      Z = __sparsefilt_features__ (X, this.TransformWeights);
    endfunction

  endmethods

endclassdef

## Soft-absolute features, normalized across examples then across features.
function F = __sparsefilt_features__ (X, W)
  A = X * W;
  F = sqrt (A .^ 2 + 1e-8);
  F = F ./ sqrt (sum (F .^ 2, 1) + 1e-8);   ## per-feature (column) normalization
  F = F ./ sqrt (sum (F .^ 2, 2) + 1e-8);   ## per-example (row) normalization
endfunction

## Sparse filtering objective (sum of features + L2 penalty) and its gradient.
function [f, g] = __sparsefilt_objective__ (wv, X, p, Q, lambda)

  W = reshape (wv, p, Q);
  A = X * W;
  F = sqrt (A .^ 2 + 1e-8);
  cn = sqrt (sum (F .^ 2, 1) + 1e-8);
  Ft = F ./ cn;
  rn = sqrt (sum (Ft .^ 2, 2) + 1e-8);
  Fh = Ft ./ rn;
  f = sum (Fh(:)) + lambda * sum (W(:) .^ 2);

  if (nargout > 1)
    dFh = ones (size (Fh));
    dFt = (dFh - sum (dFh .* Fh, 2) .* Fh) ./ rn;
    dF = (dFt - sum (dFt .* Ft, 1) .* Ft) ./ cn;
    dA = dF .* (A ./ F);
    g = X' * dA + 2 * lambda * W;
    g = g(:);
  endif

endfunction

%!shared X, W0
%! X = reshape (mod ((1:60)*7, 13), 12, 5) - 6;
%! W0 = reshape (mod ((1:15)*3, 7), 5, 3) - 3;

%!test
%! ## Construct the object directly and check its properties.
%! Mdl = SparseFiltering (X, 3, "InitialTransformWeights", W0, ...
%!                        "Lambda", 1, "IterationLimit", 500);
%! assert_equal (isa (Mdl, "SparseFiltering"), true);
%! assert_equal (Mdl.NumPredictors, 5);
%! assert_equal (Mdl.NumLearnedFeatures, 3);
%! assert_equal (size (Mdl.TransformWeights), [5, 3]);
%! assert_equal (Mdl.InitialTransformWeights, W0);
%! assert_equal (isfield (Mdl.FitInfo, "Objective"), true);

%!test
%! ## The transform method returns nonnegative features in [0, 1].
%! Mdl = SparseFiltering (X, 3, "InitialTransformWeights", W0, ...
%!                        "IterationLimit", 500);
%! Z = transform (Mdl, X);
%! assert_equal (size (Z), [12, 3]);
%! assert_equal (all (Z(:) >= 0) && all (Z(:) <= 1 + 1e-12), true);

%!test
%! ## The transform method standardizes new data when the model does.
%! Mdl = SparseFiltering (X, 2, "Standardize", true, "IterationLimit", 100, ...
%!                        "InitialTransformWeights", W0(:,1:2));
%! Z = transform (Mdl, X);
%! assert_equal (size (Z), [12, 2]);
%! assert_equal (isempty (Mdl.Mu), false);
%! assert_equal (size (Mdl.Mu), [5, 1]);
%! assert_equal (size (Mdl.Sigma), [5, 1]);

%!test
%! ## The default constructor returns an empty model.
%! Mdl = SparseFiltering ();
%! assert_equal (isempty (Mdl.TransformWeights), true);

%!error<sparsefilt: unknown parameter name 'bogus'.> ...
%! SparseFiltering (ones (6, 5), 2, "bogus", 1)
