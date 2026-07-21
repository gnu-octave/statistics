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
## @deftypefn  {statistics} {@var{Mdl} =} ocsvm (@var{X})
## @deftypefnx {statistics} {[@var{Mdl}, @var{tf}] =} ocsvm (@var{X})
## @deftypefnx {statistics} {[@var{Mdl}, @var{tf}, @var{scores}] =} ocsvm (@var{X})
## @deftypefnx {statistics} {[@dots{}] =} ocsvm (@dots{}, @var{name}, @var{value})
##
## Detect anomalies with a one-class support vector machine.
##
## @code{@var{Mdl} = ocsvm (@var{X})} fits a one-class support vector machine to
## the @math{N}-by-@math{P} matrix @var{X}, whose rows are observations and
## columns are variables, and returns a @code{OneClassSVM} object @var{Mdl}.
##
## @code{[@var{Mdl}, @var{tf}, @var{scores}] = ocsvm (@var{X})} also returns the
## @math{N}-by-1 logical vector @var{tf} flagging the anomalous observations and
## the @math{N}-by-1 vector @var{scores} of anomaly scores.  A higher score
## indicates an observation that lies further outside the boundary enclosing the
## data, and is therefore more likely to be an anomaly.
##
## The observations are mapped to a randomized feature space that approximates a
## Gaussian kernel of scale @var{KernelScale} using @var{NumExpansionDimensions}
## features, and a linear one-class boundary is fitted there with ridge
## regularization of strength @var{Lambda}.
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.34 0.66
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'KernelScale'} @tab the scale of the approximated Gaussian
## kernel, a positive scalar or @qcode{'auto'} (default).
##
## @item @qcode{'Lambda'} @tab the ridge regularization strength, a nonnegative
## scalar or @qcode{'auto'} (default).
##
## @item @qcode{'NumExpansionDimensions'} @tab the number of expanded feature
## dimensions, a positive integer or @qcode{'auto'} (default).
##
## @item @qcode{'StandardizeData'} @tab a logical scalar (default
## @code{false}); when @code{true} each predictor is centered and scaled and the
## means and standard deviations are stored in @code{@var{Mdl}.Mu} and
## @code{@var{Mdl}.Sigma}.
##
## @item @qcode{'ContaminationFraction'} @tab the assumed fraction of anomalies
## in @var{X}, a scalar in @math{[0, 1]} (default 0).  It sets
## @code{@var{Mdl}.ScoreThreshold} to @code{quantile (@var{scores}, 1 -
## @var{ContaminationFraction})}; when it is 0 the threshold is the maximum
## score and no training observation is flagged.
## @end multitable
##
## The feature expansion uses random projections, so the scores depend on the
## state of the random number generator and are not reproducible across runs
## unless the generator is seeded.  The @qcode{'auto'} selections and the fitted
## model differ from MATLAB's implementation, which uses a different feature
## expansion and solver.
##
## For a deterministic, classic one-class support vector machine (a
## nu-SVM with an exact kernel, computed through @code{libsvm}), use
## @code{fitcsvm} with a single class in the response or with the @qcode{'Nu'}
## name-value argument; that path returns a @code{ClassificationSVM} object
## whose @code{predict} method labels observations, rather than the
## anomaly-scoring interface provided here.
##
## Use the @code{isanomaly} method of @var{Mdl} to detect anomalies in new data.
##
## @seealso{OneClassSVM, isanomaly, iforest, lof, fitcsvm}
## @end deftypefn

function [Mdl, tf, scores] = ocsvm (X, varargin)

  if (nargin < 1)
    error ("ocsvm: too few input arguments.");
  endif

  Mdl = OneClassSVM (X, varargin{:});
  tf = Mdl.tf_;
  scores = Mdl.scores_;

endfunction

%!demo
%! ## Flag a handful of outliers around a Gaussian cluster.
%! X = [randn(200,2); 5 + randn(10,2)];
%! [Mdl, tf, scores] = ocsvm (X, "KernelScale", 2, ...
%!                            "ContaminationFraction", 0.05);
%! gscatter (X(:,1), X(:,2), tf);
%! title ("ocsvm: inliers vs. flagged anomalies");

## Well-separated outliers score higher and are flagged (RNG seeded)
%!test
%! rand ("state", 42);
%! randn ("state", 42);
%! X = [randn(60,2)*0.3; 10 10; -9 8; 8 -10];
%! [Mdl, tf, scores] = ocsvm (X, "KernelScale", 2, ...
%!                            "NumExpansionDimensions", 128, ...
%!                            "ContaminationFraction", 3/63);
%! assert (all (scores(61:63) > max (scores(1:60))));   # outliers rank highest
%! assert_equal (tf, logical ([false(60,1); true(3,1)]));

## ScoreThreshold follows the quantile rule; contamination 0 flags none
%!test
%! rand ("state", 7);
%! randn ("state", 7);
%! X = [randn(50,2)*0.3; 9 9; -8 7];
%! [Mdl, tf, scores] = ocsvm (X, "KernelScale", 2);
%! assert_equal (Mdl.ContaminationFraction, 0);
%! assert_equal (Mdl.ScoreThreshold, max (scores), 1e-12);
%! assert_equal (tf, false (52, 1));
%! [Mdl2, tf2, s2] = ocsvm (X, "KernelScale", 2, "ContaminationFraction", 0.1);
%! assert_equal (Mdl2.ScoreThreshold, quantile (s2, 0.9), 1e-12);
%! assert_equal (tf2, s2 > Mdl2.ScoreThreshold);

## 'auto' defaults are resolved to concrete numeric values
%!test
%! rand ("state", 1);
%! randn ("state", 1);
%! X = randn (64, 3);
%! Mdl = ocsvm (X);
%! assert_equal (Mdl.KernelScale, 1);
%! assert_equal (Mdl.Lambda, 1/64, 1e-12);
%! assert_equal (Mdl.NumExpansionDimensions, 64);   # 2^ceil(log2(64))

## StandardizeData stores Mu and Sigma and isanomaly reuses them
%!test
%! rand ("state", 3);
%! randn ("state", 3);
%! X = [randn(60,2)*[3 0; 0 0.2] + [5 -2]; 20 5];
%! Mdl = ocsvm (X, "KernelScale", 2, "StandardizeData", true);
%! assert_equal (numel (Mdl.Mu), 2);
%! assert_equal (numel (Mdl.Sigma), 2);
%! [tf, scores] = isanomaly (Mdl, [5 -2; 20 5]);
%! assert (scores(2) > scores(1));                 # the outlier scores higher

## isanomaly scores new observations against the trained model
%!test
%! rand ("state", 5);
%! randn ("state", 5);
%! X = [randn(60,2)*0.3; 10 10; -9 8];
%! Mdl = ocsvm (X, "KernelScale", 2, "ContaminationFraction", 2/62);
%! [tf, scores] = isanomaly (Mdl, [0 0; 15 15]);
%! assert (scores(2) > scores(1));
%! assert_equal (tf, logical ([false; true]));

## Test input validation
%!error <ocsvm: too few input arguments.> ocsvm ()
%!error <ocsvm: X must be a nonempty real numeric matrix.> ocsvm ([])
%!error <ocsvm: X must be a nonempty real numeric matrix.> ocsvm ("a")
%!error <ocsvm: each NAME must be followed by a VALUE.> ...
%! ocsvm (randn (10,2), "KernelScale")
%!error <ocsvm: unknown parameter name 'foo'.> ...
%! ocsvm (randn (10,2), "foo", "bar")
%!error <ocsvm: KERNELSCALE must be a positive scalar or 'auto'.> ...
%! ocsvm (randn (10,2), "KernelScale", -1)
%!error <ocsvm: LAMBDA must be a nonnegative scalar or 'auto'.> ...
%! ocsvm (randn (10,2), "Lambda", -1)
%!error <ocsvm: NUMEXPANSIONDIMENSIONS must be a positive integer or 'auto'.> ...
%! ocsvm (randn (10,2), "NumExpansionDimensions", 0)
%!error <ocsvm: CONTAMINATIONFRACTION must be a scalar in .0, 1..> ...
%! ocsvm (randn (10,2), "ContaminationFraction", 2)
