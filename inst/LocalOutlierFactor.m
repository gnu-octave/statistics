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
## @deftypefn {statistics} {@var{Mdl} =} LocalOutlierFactor (@var{X})
## @deftypefnx {statistics} {@var{Mdl} =} LocalOutlierFactor (@var{X}, @var{name}, @var{value})
##
## Local Outlier Factor model for anomaly detection.
##
## A @code{LocalOutlierFactor} object stores a Local Outlier Factor (LOF) model
## fitted to a set of observations, and detects anomalies among those or new
## observations through the @code{isanomaly} method.  Create a model with the
## @code{lof} function rather than by calling this constructor directly.
##
## The LOF of an observation compares its local density with the local density
## of its neighbors; a value near 1 indicates an inlier, whereas a value well
## above 1 indicates an outlier that lies in a sparser region than its
## neighbors.
##
## @seealso{lof, isanomaly}
## @end deftypefn

classdef LocalOutlierFactor

  properties (GetAccess = public, SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {LocalOutlierFactor} {property} NumNeighbors
    ## The number of nearest neighbors used to compute the Local Outlier Factor.
    ## @end deftp
    NumNeighbors = [];

    ## -*- texinfo -*-
    ## @deftp {LocalOutlierFactor} {property} Distance
    ## The distance metric used to find neighbors, as a character vector.
    ## @end deftp
    Distance = [];

    ## -*- texinfo -*-
    ## @deftp {LocalOutlierFactor} {property} ContaminationFraction
    ## The assumed fraction of anomalies in the training data, in @math{[0, 1]}.
    ## @end deftp
    ContaminationFraction = [];

    ## -*- texinfo -*-
    ## @deftp {LocalOutlierFactor} {property} ScoreThreshold
    ## The score above which an observation is flagged as an anomaly.
    ## @end deftp
    ScoreThreshold = [];

  endproperties

  properties (GetAccess = public, SetAccess = private, Hidden)
    X_             = [];    # training observations
    DistParameter_ = [];    # metric parameter (exponent or covariance)
    kdist_         = [];    # k-distance of each training observation
    lrd_           = [];    # local reachability density of each observation
    scores_        = [];    # LOF scores of the training observations
    tf_            = [];    # anomaly flags of the training observations
  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {LocalOutlierFactor} {@var{Mdl} =} LocalOutlierFactor (@var{X})
    ## @deftypefnx {LocalOutlierFactor} {@var{Mdl} =} LocalOutlierFactor (@var{X}, @var{name}, @var{value})
    ##
    ## Fit a Local Outlier Factor model to the @math{N}-by-@math{P} matrix
    ## @var{X}.  Prefer the @code{lof} function to this constructor.
    ##
    ## @end deftypefn
    function obj = LocalOutlierFactor (X, varargin)

      if (nargin < 1)
        error ("lof: too few input arguments.");
      endif
      if (! isnumeric (X) || ! isreal (X) || ndims (X) != 2 || isempty (X))
        error ("lof: X must be a nonempty real numeric matrix.");
      endif
      [n, p] = size (X);

      ## Defaults
      k        = min (20, rows (unique (X, "rows")) - 1);
      distance = "euclidean";
      contam   = 0;
      exponent = 2;
      covmat   = [];

      ## Parse Name-Value pairs
      if (mod (numel (varargin), 2) != 0)
        error ("lof: each NAME must be followed by a VALUE.");
      endif
      while (numel (varargin) > 0)
        name = varargin{1};
        val  = varargin{2};
        if (! ischar (name))
          error ("lof: optional argument names must be strings.");
        endif
        switch (tolower (name))
          case "numneighbors"
            k = val;
          case "distance"
            distance = tolower (val);
          case "contaminationfraction"
            contam = val;
          case "exponent"
            exponent = val;
          case "cov"
            covmat = val;
          otherwise
            error ("lof: unknown parameter name '%s'.", name);
        endswitch
        varargin(1:2) = [];
      endwhile

      ## Validate options
      metrics = {"euclidean", "seuclidean", "mahalanobis", "cityblock", ...
                 "minkowski", "chebychev", "cosine", "correlation", ...
                 "hamming", "jaccard", "spearman"};
      if (! any (strcmp (distance, metrics)))
        error ("lof: unsupported distance metric '%s'.", distance);
      endif
      if (! isscalar (k) || ! isnumeric (k) || k < 1 || fix (k) != k || k >= n)
        error ("lof: NUMNEIGHBORS must be a positive integer less than N.");
      endif
      if (! isscalar (contam) || ! isnumeric (contam) || contam < 0
                              || contam > 1)
        error ("lof: CONTAMINATIONFRACTION must be a scalar in [0, 1].");
      endif

      ## Select the metric parameter for pdist2.
      if (strcmp (distance, "minkowski"))
        distparam = exponent;
      elseif (strcmp (distance, "mahalanobis"))
        distparam = covmat;
      else
        distparam = [];
      endif

      ## Pairwise distances, excluding each point from its own neighborhood.
      D = LocalOutlierFactor.pdists_ (X, X, distance, distparam);
      D(1:(n + 1):end) = Inf;
      [kdist, lrd, scores] = LocalOutlierFactor.lofTrain_ (D, k);

      if (contam == 0)
        thr = max (scores);
      else
        thr = quantile (scores, 1 - contam);
      endif

      obj.NumNeighbors          = k;
      obj.Distance              = distance;
      obj.ContaminationFraction = contam;
      obj.ScoreThreshold        = thr;
      obj.X_                    = X;
      obj.DistParameter_        = distparam;
      obj.kdist_                = kdist;
      obj.lrd_                  = lrd;
      obj.scores_               = scores;
      obj.tf_                   = scores > thr;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LocalOutlierFactor} {@var{tf} =} isanomaly (@var{Mdl}, @var{Xnew})
    ## @deftypefnx {LocalOutlierFactor} {[@var{tf}, @var{scores}] =} isanomaly (@var{Mdl}, @var{Xnew})
    ## @deftypefnx {LocalOutlierFactor} {[@dots{}] =} isanomaly (@dots{}, @qcode{'ScoreThreshold'}, @var{t})
    ##
    ## Detect anomalies in the new observations @var{Xnew} using the fitted
    ## model @var{Mdl}.  Returns the logical vector @var{tf} flagging anomalies
    ## and the Local Outlier Factor @var{scores}, each computed from the nearest
    ## neighbors of @var{Xnew} in the training data.  The threshold defaults to
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
      if (columns (Xnew) != columns (obj.X_))
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

      Dq = LocalOutlierFactor.pdists_ (Xnew, obj.X_, obj.Distance, ...
                                       obj.DistParameter_);
      scores = LocalOutlierFactor.lofQuery_ (Dq, obj.NumNeighbors, ...
                                             obj.kdist_, obj.lrd_);
      tf = scores > thr;

    endfunction

    function disp (obj)
      printf ("  LocalOutlierFactor model\n");
      printf ("    NumNeighbors:          %d\n", obj.NumNeighbors);
      printf ("    Distance:              %s\n", obj.Distance);
      printf ("    ContaminationFraction: %g\n", obj.ContaminationFraction);
      printf ("    ScoreThreshold:        %g\n", obj.ScoreThreshold);
    endfunction

  endmethods

  methods (Static, Access = private)

    ## Pairwise distances between the rows of A and B in the given metric.
    function D = pdists_ (A, B, distance, distparam)
      if (isempty (distparam))
        D = pdist2 (A, B, distance);
      else
        D = pdist2 (A, B, distance, distparam);
      endif
    endfunction

    ## LOF of the training observations from the self-excluded distance matrix.
    function [kdist, lrd, scores] = lofTrain_ (D, k)
      n = rows (D);
      [Ds, ord] = sort (D, 2);
      NB = ord(:, 1:k);
      kdist = Ds(:, k);
      lrd = zeros (n, 1);
      for i = 1:n
        rd = max (kdist(NB(i,:)), D(i, NB(i,:))');
        lrd(i) = 1 / mean (rd);
      endfor
      scores = zeros (n, 1);
      for i = 1:n
        scores(i) = mean (lrd(NB(i,:))) / lrd(i);
      endfor
    endfunction

    ## LOF of query rows against the training neighbors, k-distances, and lrd.
    function scores = lofQuery_ (Dq, k, kdist, lrd)
      m = rows (Dq);
      [~, ord] = sort (Dq, 2);
      NB = ord(:, 1:k);
      scores = zeros (m, 1);
      for i = 1:m
        rd = max (kdist(NB(i,:)), Dq(i, NB(i,:))');
        lrdq = 1 / mean (rd);
        scores(i) = mean (lrd(NB(i,:))) / lrdq;
      endfor
    endfunction

  endmethods

endclassdef
