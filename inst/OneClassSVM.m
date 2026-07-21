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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {statistics} {@var{Mdl} =} OneClassSVM (@var{X})
## @deftypefnx {statistics} {@var{Mdl} =} OneClassSVM (@var{X}, @var{name}, @var{value})
##
## One-class support vector machine model for anomaly detection.
##
## A @code{OneClassSVM} object stores a one-class support vector machine fitted
## to a set of observations in an expanded feature space, and detects anomalies
## through the @code{isanomaly} method.  Create a model with the @code{ocsvm}
## function rather than by calling this constructor directly.
##
## The model maps the data to a randomized feature space that approximates a
## Gaussian kernel and fits a linear boundary that encloses the bulk of the
## observations; points outside the boundary receive higher anomaly scores.
##
## @seealso{ocsvm, isanomaly}
## @end deftypefn

classdef OneClassSVM

  properties (GetAccess = public, SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {OneClassSVM} {property} KernelScale
    ## The scale of the Gaussian kernel approximated by the feature expansion.
    ## @end deftp
    KernelScale = [];

    ## -*- texinfo -*-
    ## @deftp {OneClassSVM} {property} Lambda
    ## The strength of the ridge (L2) regularization term.
    ## @end deftp
    Lambda = [];

    ## -*- texinfo -*-
    ## @deftp {OneClassSVM} {property} NumExpansionDimensions
    ## The number of dimensions of the expanded feature space.
    ## @end deftp
    NumExpansionDimensions = [];

    ## -*- texinfo -*-
    ## @deftp {OneClassSVM} {property} Mu
    ## The predictor means used for standardization, or empty.
    ## @end deftp
    Mu = [];

    ## -*- texinfo -*-
    ## @deftp {OneClassSVM} {property} Sigma
    ## The predictor standard deviations used for standardization, or empty.
    ## @end deftp
    Sigma = [];

    ## -*- texinfo -*-
    ## @deftp {OneClassSVM} {property} ContaminationFraction
    ## The assumed fraction of anomalies in the training data, in @math{[0, 1]}.
    ## @end deftp
    ContaminationFraction = [];

    ## -*- texinfo -*-
    ## @deftp {OneClassSVM} {property} ScoreThreshold
    ## The score above which an observation is flagged as an anomaly.
    ## @end deftp
    ScoreThreshold = [];

  endproperties

  properties (GetAccess = public, SetAccess = private, Hidden)
    W_      = [];    # random projection weights (D-by-P)
    b_      = [];    # random phases (D-by-1)
    beta_   = [];    # fitted linear coefficients (D-by-1)
    scores_ = [];    # anomaly scores of the training observations
    tf_     = [];    # anomaly flags of the training observations
  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {OneClassSVM} {@var{Mdl} =} OneClassSVM (@var{X})
    ## @deftypefnx {OneClassSVM} {@var{Mdl} =} OneClassSVM (@var{X}, @var{name}, @var{value})
    ##
    ## Fit a one-class support vector machine to the @math{N}-by-@math{P} matrix
    ## @var{X}.  Prefer the @code{ocsvm} function to this constructor.
    ##
    ## @end deftypefn
    function obj = OneClassSVM (X, varargin)

      if (nargin < 1)
        error ("ocsvm: too few input arguments.");
      endif
      if (! isnumeric (X) || ! isreal (X) || ndims (X) != 2 || isempty (X))
        error ("ocsvm: X must be a nonempty real numeric matrix.");
      endif
      [n, p] = size (X);

      ## Defaults ("auto" values resolved below)
      kernelscale = "auto";
      lambda      = "auto";
      numdims     = "auto";
      standardize = false;
      contam      = 0;

      ## Parse Name-Value pairs
      if (mod (numel (varargin), 2) != 0)
        error ("ocsvm: each NAME must be followed by a VALUE.");
      endif
      while (numel (varargin) > 0)
        name = varargin{1};
        val  = varargin{2};
        if (! ischar (name))
          error ("ocsvm: optional argument names must be strings.");
        endif
        switch (tolower (name))
          case "kernelscale"
            kernelscale = val;
          case "lambda"
            lambda = val;
          case "numexpansiondimensions"
            numdims = val;
          case "standardizedata"
            standardize = logical (val);
          case "contaminationfraction"
            contam = val;
          otherwise
            error ("ocsvm: unknown parameter name '%s'.", name);
        endswitch
        varargin(1:2) = [];
      endwhile

      ## Resolve and validate the "auto" defaults.
      if (ischar (kernelscale) && strcmpi (kernelscale, "auto"))
        kernelscale = 1;
      elseif (! isscalar (kernelscale) || ! isnumeric (kernelscale)
              || kernelscale <= 0)
        error ("ocsvm: KERNELSCALE must be a positive scalar or 'auto'.");
      endif
      if (ischar (lambda) && strcmpi (lambda, "auto"))
        lambda = 1 / n;
      elseif (! isscalar (lambda) || ! isnumeric (lambda) || lambda < 0)
        error ("ocsvm: LAMBDA must be a nonnegative scalar or 'auto'.");
      endif
      if (ischar (numdims) && strcmpi (numdims, "auto"))
        numdims = 2 ^ ceil (log2 (max (n, 2)));
      elseif (! isscalar (numdims) || ! isnumeric (numdims) || numdims < 1
              || fix (numdims) != numdims)
        error (strcat ("ocsvm: NUMEXPANSIONDIMENSIONS must be a positive", ...
                       " integer or 'auto'."));
      endif
      if (! isscalar (contam) || ! isnumeric (contam) || contam < 0
                              || contam > 1)
        error ("ocsvm: CONTAMINATIONFRACTION must be a scalar in [0, 1].");
      endif

      ## Optionally standardize the predictors.
      if (standardize)
        mu = mean (X, 1);
        sigma = std (X, 0, 1);
        sigma(sigma == 0) = 1;
        Xs = (X - mu) ./ sigma;
      else
        mu = [];
        sigma = [];
        Xs = X;
      endif

      ## Random Fourier feature map approximating a Gaussian kernel, then a
      ## squared-hinge one-class SVM with ridge regularization.
      W = randn (numdims, p) / kernelscale;
      b = 2 * pi * rand (numdims, 1);
      Phi = OneClassSVM.features_ (Xs, W, b);
      beta = OneClassSVM.fitOneClass_ (Phi, lambda);
      scores = - (Phi * beta);

      if (contam == 0)
        thr = max (scores);
      else
        thr = quantile (scores, 1 - contam);
      endif

      obj.KernelScale            = kernelscale;
      obj.Lambda                 = lambda;
      obj.NumExpansionDimensions = numdims;
      obj.Mu                     = mu;
      obj.Sigma                  = sigma;
      obj.ContaminationFraction  = contam;
      obj.ScoreThreshold         = thr;
      obj.W_                     = W;
      obj.b_                     = b;
      obj.beta_                  = beta;
      obj.scores_                = scores;
      obj.tf_                    = scores > thr;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {OneClassSVM} {@var{tf} =} isanomaly (@var{Mdl}, @var{Xnew})
    ## @deftypefnx {OneClassSVM} {[@var{tf}, @var{scores}] =} isanomaly (@var{Mdl}, @var{Xnew})
    ## @deftypefnx {OneClassSVM} {[@dots{}] =} isanomaly (@dots{}, @qcode{'ScoreThreshold'}, @var{t})
    ##
    ## Detect anomalies in the new observations @var{Xnew} using the fitted
    ## model @var{Mdl}.  Returns the logical vector @var{tf} flagging anomalies
    ## and the anomaly @var{scores}.  The threshold defaults to
    ## @code{@var{Mdl}.ScoreThreshold} and may be overridden with the
    ## @qcode{'ScoreThreshold'} name-value argument.
    ##
    ## @end deftypefn
    function [tf, scores] = isanomaly (obj, Xnew, varargin)

      if (nargin < 2)
        error ("isanomaly: too few input arguments.");
      endif
      if (! isnumeric (Xnew) || ! isreal (Xnew) || ndims (Xnew) != 2
                             || isempty (Xnew))
        error ("isanomaly: XNEW must be a nonempty real numeric matrix.");
      endif
      if (columns (Xnew) != columns (obj.W_))
        error ("isanomaly: XNEW must have the same number of columns as X.");
      endif

      thr = obj.ScoreThreshold;
      if (mod (numel (varargin), 2) != 0)
        error ("isanomaly: each NAME must be followed by a VALUE.");
      endif
      while (numel (varargin) > 0)
        if (! ischar (varargin{1}) || ! strcmpi (varargin{1}, "ScoreThreshold"))
          error ("isanomaly: unknown parameter name.");
        endif
        thr = varargin{2};
        varargin(1:2) = [];
      endwhile

      if (! isempty (obj.Mu))
        Xnew = (Xnew - obj.Mu) ./ obj.Sigma;
      endif
      Phi = OneClassSVM.features_ (Xnew, obj.W_, obj.b_);
      scores = - (Phi * obj.beta_);
      tf = scores > thr;

    endfunction

    function disp (obj)
      printf ("  OneClassSVM model\n");
      printf ("    KernelScale:            %g\n", obj.KernelScale);
      printf ("    Lambda:                 %g\n", obj.Lambda);
      printf ("    NumExpansionDimensions: %d\n", obj.NumExpansionDimensions);
      printf ("    ContaminationFraction:  %g\n", obj.ContaminationFraction);
      printf ("    ScoreThreshold:         %g\n", obj.ScoreThreshold);
    endfunction

  endmethods

  methods (Static, Access = private)

    ## Random Fourier features approximating a Gaussian kernel.
    function Phi = features_ (X, W, b)
      D = rows (W);
      Phi = sqrt (2 / D) * cos (X * W' + b');
    endfunction

    ## Squared-hinge one-class SVM (unit margin from the origin) with ridge
    ## regularization, minimized by Newton iterations over the active set.
    function beta = fitOneClass_ (Phi, lambda)
      [n, D] = size (Phi);
      beta = zeros (D, 1);
      I = eye (D);
      for it = 1:100
        d = Phi * beta;
        active = d < 1;                 # observations inside the margin
        Pa = Phi(active, :);
        r = 1 - d(active);
        g = lambda * beta - (Pa' * r) / n;
        H = lambda * I + (Pa' * Pa) / n;
        step = H \ g;
        beta = beta - step;
        if (norm (step) < 1e-10)
          break;
        endif
      endfor
    endfunction

  endmethods

endclassdef

## Direct construction and isanomaly dispatch on the class
%!test
%! rand ("state", 8);
%! randn ("state", 8);
%! X = [randn(60,2)*0.3; 9 9; -8 7];
%! Mdl = OneClassSVM (X, "KernelScale", 2, "NumExpansionDimensions", 64);
%! assert (isa (Mdl, "OneClassSVM"));
%! assert_equal (Mdl.NumExpansionDimensions, 64);
%! [tf, scores] = isanomaly (Mdl, [0 0; 12 12]);
%! assert (scores(2) > scores(1));                 # far point is more anomalous
%! assert (islogical (tf));
