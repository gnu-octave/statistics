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
## An object of this class has the following properties:
##
## @itemize
## @item @qcode{TransformWeights} — the @math{P * Q} matrix of learned
## transformation weights (unit-length columns).
## @item @qcode{Mu}, @qcode{Sigma} — the per-predictor mean and standard
## deviation used when @qcode{'Standardize'} is @qcode{true} (empty otherwise).
## @item @qcode{FitInfo} — a structure with the number of @qcode{Iteration}s and
## the final @qcode{Objective} value of the fit.
## @item @qcode{ModelParameters} — a structure of the options used for the fit.
## @item @qcode{NumPredictors}, @qcode{NumLearnedFeatures} — the input dimension
## @math{P} and the number of learned features @math{Q}.
## @item @qcode{InitialTransformWeights} — the starting weights used by the fit.
## @end itemize
##
## @seealso{rica, sparsefilt}
## @end deftp

classdef ReconstructionICA

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
                     "GradientTolerance", 1e-6, "StepTolerance", 1e-6);
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
          otherwise
            error ("rica: unknown parameter name '%s'.", name);
        endswitch
      endfor

      if (! strcmp (opts.ContrastFcn, "logcosh"))
        error ("rica: only the 'logcosh' ContrastFcn is supported.");
      endif

      ## Standardize (center and scale) the data if requested
      if (opts.Standardize)
        this.Mu = mean (X);
        this.Sigma = std (X);
        sig = this.Sigma;
        sig(sig == 0) = 1;
        X = (X - this.Mu) ./ sig;
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
      ofun = @(wv) __rica_objective__ (wv, X, p, Q, lambda);
      fmopts = optimset ("GradObj", "on", "MaxIter", opts.IterationLimit, ...
                         "TolFun", 1e-10, "TolX", 1e-10, "Display", "off");
      [wv, fval, ~, output] = fminunc (ofun, W0(:), fmopts);

      W = reshape (wv, p, Q);
      W = W ./ sqrt (sum (W .^ 2, 1));   ## unit-length columns
      this.TransformWeights = W;

      this.NumPredictors = p;
      this.NumLearnedFeatures = Q;
      niter = output.iterations;
      this.FitInfo = struct ("Iteration", niter, "Objective", fval);
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
        sig = this.Sigma;
        sig(sig == 0) = 1;
        X = (X - this.Mu) ./ sig;
      endif
      Z = X * this.TransformWeights;
    endfunction

  endmethods

endclassdef

## RICA objective and gradient with respect to the raw (un-normalized) weights.
function [f, g] = __rica_objective__ (wv, X, p, Q, lambda)

  W = reshape (wv, p, Q);
  nrm = sqrt (sum (W .^ 2, 1));
  Wn = W ./ nrm;
  Z = X * Wn;
  E = X * (Wn * Wn') - X;
  f = lambda * sum (E(:) .^ 2) + sum (sum (0.5 * log (cosh (2 * Z))));

  if (nargout > 1)
    gWn = lambda * 2 * (X' * E * Wn + E' * X * Wn) + X' * tanh (2 * Z);
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
%! Xs = (X - Mdl.Mu) ./ Mdl.Sigma;
%! assert_equal (Z, Xs * Mdl.TransformWeights, 1e-12);

%!test
%! ## The default constructor returns an empty model.
%! Mdl = ReconstructionICA ();
%! assert_equal (isempty (Mdl.TransformWeights), true);

%!error<rica: only the 'logcosh' ContrastFcn is supported.> ...
%! ReconstructionICA (ones (6, 5), 2, "ContrastFcn", "exp")
%!error<rica: unknown parameter name 'bogus'.> ...
%! ReconstructionICA (ones (6, 5), 2, "bogus", 1)
