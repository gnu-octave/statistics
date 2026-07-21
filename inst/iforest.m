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
## @deftypefn  {statistics} {@var{Mdl} =} iforest (@var{X})
## @deftypefnx {statistics} {[@var{Mdl}, @var{tf}] =} iforest (@var{X})
## @deftypefnx {statistics} {[@var{Mdl}, @var{tf}, @var{scores}] =} iforest (@var{X})
## @deftypefnx {statistics} {[@dots{}] =} iforest (@dots{}, @var{name}, @var{value})
##
## Detect anomalies with an isolation forest.
##
## @code{@var{Mdl} = iforest (@var{X})} fits an isolation forest to the
## @math{N}-by-@math{P} matrix @var{X}, whose rows are observations and columns
## are variables, and returns an @code{IsolationForest} object @var{Mdl}.
##
## @code{[@var{Mdl}, @var{tf}, @var{scores}] = iforest (@var{X})} also returns
## the @math{N}-by-1 logical vector @var{tf} flagging the anomalous observations
## and the @math{N}-by-1 vector @var{scores} of anomaly scores in the range
## @math{[0, 1]}.  A higher score indicates an observation that is more easily
## isolated, and therefore more likely to be an anomaly.
##
## The score of an observation is @code{2^(-E[h] / c)}, where @math{E[h]} is its
## average path length over the isolation trees and @math{c} is the expected
## path length of an unsuccessful search in a binary tree of
## @var{NumObservationsPerLearner} nodes.  Each tree is grown from a random
## subsample of the data by recursively splitting on a random variable at a
## random value, so anomalies, being easier to isolate, obtain shorter paths.
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.34 0.66
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'NumLearners'} @tab the number of isolation trees, a positive
## integer (default 100).
##
## @item @qcode{'NumObservationsPerLearner'} @tab the subsample size used to
## grow each tree, an integer in @math{[3, N]} (default @code{min (@var{N},
## 256)}).
##
## @item @qcode{'ContaminationFraction'} @tab the assumed fraction of anomalies
## in @var{X}, a scalar in @math{[0, 1]} (default 0).  It sets
## @code{@var{Mdl}.ScoreThreshold} to @code{quantile (@var{scores}, 1 -
## @var{ContaminationFraction})}; when it is 0 the threshold is the maximum
## score and no training observation is flagged.
## @end multitable
##
## Because the trees are grown from random subsamples and random splits, the
## scores depend on the state of the random number generator and are not
## reproducible across runs unless the generator is seeded.
##
## Use the @code{isanomaly} method of @var{Mdl} to detect anomalies in new data.
##
## @seealso{IsolationForest, isanomaly, lof, robustcov}
## @end deftypefn

function [Mdl, tf, scores] = iforest (X, varargin)

  if (nargin < 1)
    error ("iforest: too few input arguments.");
  endif

  Mdl = IsolationForest (X, varargin{:});
  tf = Mdl.tf_;
  scores = Mdl.scores_;

endfunction

%!demo
%! ## Flag a handful of outliers around a Gaussian cluster.
%! X = [randn(200,2); 6 + randn(10,2)];
%! [Mdl, tf, scores] = iforest (X, "ContaminationFraction", 0.05);
%! gscatter (X(:,1), X(:,2), tf);
%! title ("iforest: inliers vs. flagged anomalies");

## Well-separated outliers score higher and are flagged (RNG seeded)
%!test
%! rand ("state", 42);
%! X = [randn(60,2)*0.3; 12 12; -11 10; 10 -12];
%! [Mdl, tf, scores] = iforest (X, "ContaminationFraction", 3/63);
%! assert_equal (Mdl.NumLearners, 100);
%! assert_equal (Mdl.NumObservationsPerLearner, 63);
%! assert (all (scores >= 0 & scores <= 1));
%! assert (all (scores(61:63) > max (scores(1:60))));   # outliers rank highest
%! assert_equal (tf, logical ([false(60,1); true(3,1)]));

## ScoreThreshold follows the quantile rule; contamination 0 flags none
%!test
%! rand ("state", 7);
%! X = [randn(50,2)*0.3; 9 9; -8 7];
%! [Mdl, tf, scores] = iforest (X);
%! assert_equal (Mdl.ContaminationFraction, 0);
%! assert_equal (Mdl.ScoreThreshold, max (scores), 1e-12);
%! assert_equal (tf, false (52, 1));
%! [Mdl2, tf2, s2] = iforest (X, "ContaminationFraction", 0.1);
%! assert_equal (Mdl2.ScoreThreshold, quantile (s2, 0.9), 1e-12);
%! assert_equal (tf2, s2 > Mdl2.ScoreThreshold);

## NumObservationsPerLearner caps at 256 for large samples
%!test
%! rand ("state", 1);
%! X = randn (400, 2);
%! Mdl = iforest (X);
%! assert_equal (Mdl.NumObservationsPerLearner, 256);

## isanomaly scores new observations against the trained forest
%!test
%! rand ("state", 3);
%! X = [randn(60,2)*0.3; 12 12; -11 10; 10 -12];
%! Mdl = iforest (X, "ContaminationFraction", 3/63);
%! [tf, scores] = isanomaly (Mdl, [0 0; 15 15]);
%! assert (all (scores >= 0 & scores <= 1));
%! assert (scores(2) > scores(1));            # far point scores higher
%! assert_equal (tf, logical ([false; true]));

## Test input validation
%!error <iforest: too few input arguments.> iforest ()
%!error <iforest: X must be a nonempty real numeric matrix.> iforest ([])
%!error <iforest: X must be a nonempty real numeric matrix.> iforest ("a")
%!error <iforest: each NAME must be followed by a VALUE.> ...
%! iforest (randn (10,2), "NumLearners")
%!error <iforest: unknown parameter name 'foo'.> ...
%! iforest (randn (10,2), "foo", "bar")
%!error <iforest: NUMLEARNERS must be a positive integer scalar.> ...
%! iforest (randn (10,2), "NumLearners", 0)
%!error <iforest: NUMOBSERVATIONSPERLEARNER must be an integer in .3, N..> ...
%! iforest (randn (10,2), "NumObservationsPerLearner", 2)
%!error <iforest: NUMOBSERVATIONSPERLEARNER must be an integer in .3, N..> ...
%! iforest (randn (10,2), "NumObservationsPerLearner", 20)
%!error <iforest: CONTAMINATIONFRACTION must be a scalar in .0, 1..> ...
%! iforest (randn (10,2), "ContaminationFraction", 2)
