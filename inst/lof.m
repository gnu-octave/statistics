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
## @deftypefn  {statistics} {@var{Mdl} =} lof (@var{X})
## @deftypefnx {statistics} {[@var{Mdl}, @var{tf}] =} lof (@var{X})
## @deftypefnx {statistics} {[@var{Mdl}, @var{tf}, @var{scores}] =} lof (@var{X})
## @deftypefnx {statistics} {[@dots{}] =} lof (@dots{}, @var{name}, @var{value})
##
## Detect anomalies with the Local Outlier Factor (LOF) method.
##
## @code{@var{Mdl} = lof (@var{X})} fits a Local Outlier Factor model to the
## @math{N}-by-@math{P} matrix @var{X}, whose rows are observations and columns
## are variables, and returns a @code{LocalOutlierFactor} object @var{Mdl}.
##
## @code{[@var{Mdl}, @var{tf}, @var{scores}] = lof (@var{X})} also returns the
## @math{N}-by-1 logical vector @var{tf} flagging the anomalous observations and
## the @math{N}-by-1 vector @var{scores} of LOF values.  A score near 1
## indicates an inlier, whereas a score well above 1 indicates an outlier lying
## in a region sparser than its neighbors.
##
## The Local Outlier Factor of an observation is the average ratio of the local
## reachability density of its @var{NumNeighbors} nearest neighbors to its own
## local reachability density, where the local reachability density is the
## inverse mean reachability distance to those neighbors and the reachability
## distance from @math{p} to @math{o} is @code{max (k-distance (o), d (p, o))}.
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.28 0.72
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'NumNeighbors'} @tab the number of nearest neighbors, a positive
## integer less than @math{N}.  The default is @code{min (20, @var{u} - 1)},
## where @var{u} is the number of unique observations.
##
## @item @qcode{'Distance'} @tab the distance metric used to find neighbors, one
## of the metrics accepted by @code{pdist2} (@qcode{'euclidean'} by default).
##
## @item @qcode{'ContaminationFraction'} @tab the assumed fraction of anomalies
## in @var{X}, a scalar in @math{[0, 1]} (default 0).  It sets
## @code{@var{Mdl}.ScoreThreshold} to @code{quantile (@var{scores}, 1 -
## @var{ContaminationFraction})}; when it is 0 the threshold is the maximum
## score and no training observation is flagged.
##
## @item @qcode{'Exponent'} @tab the Minkowski distance exponent (default 2),
## used only with the @qcode{'minkowski'} distance.
##
## @item @qcode{'Cov'} @tab the covariance matrix used only with the
## @qcode{'mahalanobis'} distance.
## @end multitable
##
## Use the @code{isanomaly} method of @var{Mdl} to detect anomalies in new data.
##
## @seealso{LocalOutlierFactor, isanomaly, dbscan, robustcov}
## @end deftypefn

function [Mdl, tf, scores] = lof (X, varargin)

  if (nargin < 1)
    error ("lof: too few input arguments.");
  endif

  Mdl = LocalOutlierFactor (X, varargin{:});
  tf = Mdl.tf_;
  scores = Mdl.scores_;

endfunction

%!demo
%! ## Flag a handful of outliers around a Gaussian cluster.
%! X = [randn(100,2); 4 + randn(6,2)];
%! [Mdl, tf, scores] = lof (X, "ContaminationFraction", 0.05);
%! gscatter (X(:,1), X(:,2), tf);
%! title ("lof: inliers vs. flagged anomalies");

## MATLAB parity: LOF scores are exact to machine precision (k = 5)
%!test
%! X = [0 0; 0.1 0.1; 0.2 -0.1; -0.1 0.2; 0.1 -0.2; -0.2 0.1; 0.15 0.05; ...
%!      -0.05 -0.15; 0.05 0.12; -0.12 -0.05; 5 5; -4 3];
%! [Mdl, tf, scores] = lof (X, "NumNeighbors", 5);
%! exp_scores = [1.066572826392391; 0.965333187918473; 1.035972400143920; ...
%!               1.019798243339781; 1.055271577857807; 1.047371774849485; ...
%!               0.941746286231341; 0.987111513710991; 0.948358372496933; ...
%!               1.004238203765983; 29.441894930160952; 20.058918161335612];
%! assert_equal (scores, exp_scores, 1e-12);
%! assert_equal (tf, false (12, 1));               # contamination 0 flags none
%! assert_equal (Mdl.ScoreThreshold, 29.441894930160952, 1e-12);

## MATLAB parity: ContaminationFraction sets the quantile threshold and flags
%!test
%! X = [0 0; 0.1 0.1; 0.2 -0.1; -0.1 0.2; 0.1 -0.2; -0.2 0.1; 0.15 0.05; ...
%!      -0.05 -0.15; 0.05 0.12; -0.12 -0.05; 5 5; -4 3];
%! [Mdl, tf, scores] = lof (X, "NumNeighbors", 5, "ContaminationFraction", 0.2);
%! assert_equal (Mdl.ScoreThreshold, 2.965807359886740, 1e-12);
%! assert_equal (tf, logical ([0;0;0;0;0;0;0;0;0;0;1;1]));

## MATLAB parity: default NumNeighbors is min (20, unique-1)
%!test
%! X = [0 0; 0.1 0.1; 0.2 -0.1; -0.1 0.2; 0.1 -0.2; -0.2 0.1; 0.15 0.05; ...
%!      -0.05 -0.15; 0.05 0.12; -0.12 -0.05; 5 5; -4 3];
%! Mdl = lof (X);
%! assert_equal (Mdl.NumNeighbors, 11);

## MATLAB parity: isanomaly reproduces LOF for new observations
%!test
%! X = [0 0; 0.1 0.1; 0.2 -0.1; -0.1 0.2; 0.1 -0.2; -0.2 0.1; 0.15 0.05; ...
%!      -0.05 -0.15; 0.05 0.12; -0.12 -0.05; 5 5; -4 3];
%! Mdl = lof (X, "NumNeighbors", 5);
%! [tf, scores] = isanomaly (Mdl, [0 0.05; 6 6; -0.1 -0.1]);
%! assert_equal (scores, [0.954484172585537; 27.876003442184082; ...
%!                        1.001784195086020], 1e-12);
%! assert_equal (tf, logical ([0; 0; 0]));         # cutoff 29.44 from training
%! [tf2, ~] = isanomaly (Mdl, [6 6], "ScoreThreshold", 5);
%! assert_equal (tf2, true);

## A non-euclidean metric is passed through to pdist2
%!test
%! X = [0 0; 0.1 0.1; 0.2 -0.1; -0.1 0.2; 0.1 -0.2; -0.2 0.1; 0.15 0.05; ...
%!      -0.05 -0.15; 0.05 0.12; -0.12 -0.05; 5 5; -4 3];
%! [~, ~, scores] = lof (X, "NumNeighbors", 5, "Distance", "cityblock");
%! assert_equal (scores(11), 34.198447151536712, 1e-10);
%! assert (scores(1) > 1 && scores(7) < 1);

## Test input validation
%!error <lof: too few input arguments.> lof ()
%!error <lof: X must be a nonempty real numeric matrix.> lof ([])
%!error <lof: X must be a nonempty real numeric matrix.> lof ("a")
%!error <lof: each NAME must be followed by a VALUE.> lof (ones (5,2), "Distance")
%!error <lof: unknown parameter name 'foo'.> lof (ones (5,2), "foo", "bar")
%!error <lof: unsupported distance metric 'taxicab'.> ...
%! lof (ones (5,2), "Distance", "taxicab")
%!error <lof: NUMNEIGHBORS must be a positive integer less than N.> ...
%! lof (ones (5,2), "NumNeighbors", 0)
%!error <lof: NUMNEIGHBORS must be a positive integer less than N.> ...
%! lof (ones (5,2), "NumNeighbors", 5)
%!error <lof: CONTAMINATIONFRACTION must be a scalar in .0, 1..> ...
%! lof (magic (5), "NumNeighbors", 2, "ContaminationFraction", 1.5)
%!error <isanomaly: XNEW must have the same number of columns as X.> ...
%! isanomaly (lof (magic (6)), ones (3,3))
