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
## @deftypefn {statistics} {@var{Mdl} =} IsolationForest (@var{X})
## @deftypefnx {statistics} {@var{Mdl} =} IsolationForest (@var{X}, @var{name}, @var{value})
##
## Isolation Forest model for anomaly detection.
##
## An @code{IsolationForest} object stores an ensemble of isolation trees fitted
## to a set of observations and detects anomalies through the @code{isanomaly}
## method.  Create a model with the @code{iforest} function rather than by
## calling this constructor directly.
##
## Anomalies are easier to isolate, so they sit closer to the root of a random
## isolation tree; the shorter its average path length across the ensemble, the
## higher an observation's anomaly score.
##
## @seealso{iforest, isanomaly}
## @end deftypefn

classdef IsolationForest

  properties (GetAccess = public, SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {IsolationForest} {property} NumLearners
    ## The number of isolation trees in the ensemble.
    ## @end deftp
    NumLearners = [];

    ## -*- texinfo -*-
    ## @deftp {IsolationForest} {property} NumObservationsPerLearner
    ## The number of observations subsampled to grow each isolation tree.
    ## @end deftp
    NumObservationsPerLearner = [];

    ## -*- texinfo -*-
    ## @deftp {IsolationForest} {property} ContaminationFraction
    ## The assumed fraction of anomalies in the training data, in @math{[0, 1]}.
    ## @end deftp
    ContaminationFraction = [];

    ## -*- texinfo -*-
    ## @deftp {IsolationForest} {property} ScoreThreshold
    ## The score above which an observation is flagged as an anomaly.
    ## @end deftp
    ScoreThreshold = [];

  endproperties

  properties (GetAccess = public, SetAccess = private, Hidden)
    Trees_   = {};    # ensemble of isolation trees (nested structs)
    scores_  = [];    # anomaly scores of the training observations
    tf_      = [];    # anomaly flags of the training observations
  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {IsolationForest} {@var{Mdl} =} IsolationForest (@var{X})
    ## @deftypefnx {IsolationForest} {@var{Mdl} =} IsolationForest (@var{X}, @var{name}, @var{value})
    ##
    ## Fit an isolation forest to the @math{N}-by-@math{P} matrix @var{X}.
    ## Prefer the @code{iforest} function to this constructor.
    ##
    ## @end deftypefn
    function obj = IsolationForest (X, varargin)

      if (nargin < 1)
        error ("iforest: too few input arguments.");
      endif
      if (! isnumeric (X) || ! isreal (X) || ndims (X) != 2 || isempty (X))
        error ("iforest: X must be a nonempty real numeric matrix.");
      endif
      [n, p] = size (X);

      ## Defaults
      numlearners = 100;
      numobs      = min (n, 256);
      contam      = 0;

      ## Parse Name-Value pairs
      if (mod (numel (varargin), 2) != 0)
        error ("iforest: each NAME must be followed by a VALUE.");
      endif
      while (numel (varargin) > 0)
        name = varargin{1};
        val  = varargin{2};
        if (! ischar (name))
          error ("iforest: optional argument names must be strings.");
        endif
        switch (tolower (name))
          case "numlearners"
            numlearners = val;
          case "numobservationsperlearner"
            numobs = val;
          case "contaminationfraction"
            contam = val;
          otherwise
            error ("iforest: unknown parameter name '%s'.", name);
        endswitch
        varargin(1:2) = [];
      endwhile

      ## Validate options
      if (! isscalar (numlearners) || ! isnumeric (numlearners)
          || numlearners < 1 || fix (numlearners) != numlearners)
        error ("iforest: NUMLEARNERS must be a positive integer scalar.");
      endif
      if (! isscalar (numobs) || ! isnumeric (numobs) || numobs < 3
          || fix (numobs) != numobs || numobs > n)
        error (strcat ("iforest: NUMOBSERVATIONSPERLEARNER must be an", ...
                       " integer in [3, N]."));
      endif
      if (! isscalar (contam) || ! isnumeric (contam) || contam < 0
                              || contam > 1)
        error ("iforest: CONTAMINATIONFRACTION must be a scalar in [0, 1].");
      endif

      ## Grow the ensemble, each tree from a random subsample of NUMOBS rows.
      maxdepth = ceil (log2 (numobs));
      trees = cell (numlearners, 1);
      for t = 1:numlearners
        idx = randperm (n, numobs);
        trees{t} = IsolationForest.buildTree_ (X(idx, :), 0, maxdepth);
      endfor

      scores = IsolationForest.scoreAll_ (X, trees, numobs);
      if (contam == 0)
        thr = max (scores);
      else
        thr = quantile (scores, 1 - contam);
      endif

      obj.NumLearners               = numlearners;
      obj.NumObservationsPerLearner = numobs;
      obj.ContaminationFraction     = contam;
      obj.ScoreThreshold            = thr;
      obj.Trees_                    = trees;
      obj.scores_                   = scores;
      obj.tf_                       = scores > thr;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {IsolationForest} {@var{tf} =} isanomaly (@var{Mdl}, @var{Xnew})
    ## @deftypefnx {IsolationForest} {[@var{tf}, @var{scores}] =} isanomaly (@var{Mdl}, @var{Xnew})
    ## @deftypefnx {IsolationForest} {[@dots{}] =} isanomaly (@dots{}, @qcode{'ScoreThreshold'}, @var{t})
    ##
    ## Detect anomalies in the new observations @var{Xnew} using the fitted
    ## model @var{Mdl}.  Returns the logical vector @var{tf} flagging anomalies
    ## and the anomaly @var{scores}, each obtained by dropping @var{Xnew}
    ## through the isolation trees.  The threshold defaults to
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

      scores = IsolationForest.scoreAll_ (Xnew, obj.Trees_, ...
                                          obj.NumObservationsPerLearner);
      tf = scores > thr;

    endfunction

    function disp (obj)
      printf ("  IsolationForest model\n");
      printf ("    NumLearners:               %d\n", obj.NumLearners);
      printf ("    NumObservationsPerLearner: %d\n", ...
              obj.NumObservationsPerLearner);
      printf ("    ContaminationFraction:     %g\n", obj.ContaminationFraction);
      printf ("    ScoreThreshold:            %g\n", obj.ScoreThreshold);
    endfunction

  endmethods

  methods (Static, Access = private)

    ## Grow one isolation tree by recursive random splits.
    function node = buildTree_ (X, depth, maxdepth)
      n = rows (X);
      if (depth >= maxdepth || n <= 1)
        node = struct ("leaf", true, "size", n);
        return;
      endif
      mn = min (X, [], 1);
      mx = max (X, [], 1);
      valid = find (mx > mn);
      if (isempty (valid))               # all rows identical
        node = struct ("leaf", true, "size", n);
        return;
      endif
      q = valid(randi (numel (valid)));
      sp = mn(q) + (mx(q) - mn(q)) * rand;
      isleft = X(:, q) < sp;
      if (! any (isleft) || all (isleft)) # degenerate split
        node = struct ("leaf", true, "size", n);
        return;
      endif
      L = IsolationForest.buildTree_ (X(isleft, :), depth + 1, maxdepth);
      R = IsolationForest.buildTree_ (X(! isleft, :), depth + 1, maxdepth);
      node = struct ("leaf", false, "feature", q, "split", sp, ...
                     "left", L, "right", R);
    endfunction

    ## Path length of a single observation down one tree.
    function h = pathLength_ (x, node, depth)
      if (node.leaf)
        h = depth + IsolationForest.cFactor_ (node.size);
        return;
      endif
      if (x(node.feature) < node.split)
        h = IsolationForest.pathLength_ (x, node.left, depth + 1);
      else
        h = IsolationForest.pathLength_ (x, node.right, depth + 1);
      endif
    endfunction

    ## Average path length of an unsuccessful search in a binary search tree.
    function c = cFactor_ (n)
      if (n <= 1)
        c = 0;
      else
        H = log (n - 1) + 0.5772156649015329;
        c = 2 * H - 2 * (n - 1) / n;
      endif
    endfunction

    ## Anomaly scores of the rows of P over the whole ensemble.
    function scores = scoreAll_ (P, trees, psi)
      m = rows (P);
      T = numel (trees);
      H = zeros (m, 1);
      for t = 1:T
        node = trees{t};
        for i = 1:m
          H(i) += IsolationForest.pathLength_ (P(i, :), node, 0);
        endfor
      endfor
      Eh = H / T;
      scores = 2 .^ (- Eh / IsolationForest.cFactor_ (psi));
    endfunction

  endmethods

endclassdef

## Direct construction and isanomaly dispatch on the class
%!test
%! rand ("state", 5);
%! X = [randn(40,2)*0.3; 10 10; -9 8];
%! Mdl = IsolationForest (X, "NumLearners", 50);
%! assert (isa (Mdl, "IsolationForest"));
%! assert_equal (Mdl.NumLearners, 50);
%! assert_equal (Mdl.NumObservationsPerLearner, 42);
%! [tf, scores] = isanomaly (Mdl, [0 0; 12 12]);
%! assert (all (scores >= 0 & scores <= 1));
%! assert (scores(2) > scores(1));                 # far point is more anomalous
%! assert (islogical (tf));
