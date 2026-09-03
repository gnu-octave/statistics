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
## @deftypefn  {statistics} {@var{Mdl} =} rica (@var{X}, @var{Q})
## @deftypefnx {statistics} {@var{Mdl} =} rica (@var{X}, @var{Q}, @var{Name}, @var{Value})
##
## Reconstruction independent component analysis (RICA) for feature extraction.
##
## @code{@var{Mdl} = rica (@var{X}, @var{Q})} learns @var{Q} features from the
## @math{N * P} data matrix @var{X} (rows are observations, columns are
## predictors) and returns a @qcode{ReconstructionICA} object @var{Mdl}.  Apply
## the learned transformation to data with @code{transform (@var{Mdl}, @var{X})},
## which returns @code{@var{X} * @var{Mdl}.TransformWeights}.
##
## The @math{P * @var{Q}} weight matrix (with unit-length columns) minimizes the
## objective
##
## @example
## @var{Lambda} * ||@var{X} * @var{W} * @var{W}' - @var{X}||_F^2
##       + sum (sum (@var{g} (@var{X} * @var{W})))
## @end example
##
## @noindent
## over the transformation weights @var{W}, combining a reconstruction cost with
## a sparsity contrast @var{g} applied elementwise and selected by
## @qcode{'ContrastFcn'}.
##
## Name/Value pairs:
##
## @table @asis
## @item @qcode{'IterationLimit'}
## Maximum number of iterations (default 1000).
##
## @item @qcode{'Lambda'}
## Weight of the reconstruction term (default 1).
##
## @item @qcode{'Standardize'}
## Logical; center and scale each predictor before fitting (default
## @qcode{false}).
##
## @item @qcode{'ContrastFcn'}
## The sparsity contrast @var{g} applied to each element @math{z} of
## @math{@var{X} * @var{W}}: @qcode{'logcosh'} (default), which is
## @math{0.5 * log (cosh (2 * z))}; @qcode{'exp'}, which is
## @math{-exp (-z^2 / 2)}; or @qcode{'sqrt'}, which is
## @math{sqrt (z^2 + 1e-8)}, a smooth stand-in for @math{abs (z)}.
##
## @item @qcode{'InitialTransformWeights'}
## A @math{P * @var{Q}} initial value for the weights.  The default is random.
##
## @item @qcode{'NonGaussianityIndicator'}
## A @var{Q}-element vector of @math{+1} and @math{-1}, one per learned
## feature, setting the sign its contrast term carries in the objective:
## @math{+1} seeks a super-Gaussian (sparse) feature and @math{-1} a
## sub-Gaussian (spread) one.  The default is all @math{+1}.
##
## @item @qcode{'GradientTolerance'}, @qcode{'StepTolerance'}
## Stop once the gradient's or the step's infinity norm falls to or below the
## given value (default @qcode{1e-6} each).  They govern the fit only under
## @qcode{'Solver', 'lbfgs'}; the @qcode{'quasinewton'} solver runs to its own
## tighter internal tolerances and records these without acting on them.
##
## @item @qcode{'Solver'}
## @qcode{'quasinewton'} (default) minimizes through Octave's @code{fminunc},
## which carries a full inverse Hessian.  @qcode{'lbfgs'} selects the
## limited-memory BFGS solver MATLAB uses, holding as many curvature pairs as
## the transform has parameters, and is several times faster here.  It stops
## where @qcode{'GradientTolerance'} and @qcode{'StepTolerance'} say to, so a
## value tighter than the default carries it further.
## @end table
##
## @subheading Note on reproducibility
##
## The RICA objective is not convex and is minimized by a quasi-Newton solver, so
## the learned weights depend on the starting point and the solver, and are only
## defined up to a permutation and sign of the feature columns.  Different runs
## (or different software, including MATLAB) may return different weights that
## nonetheless describe an equally valid feature transformation.  Fix
## @qcode{'InitialTransformWeights'} for a reproducible result.
##
## @seealso{ReconstructionICA, sparsefilt, pca}
## @end deftypefn

function Mdl = rica (X, Q, varargin)

  ## Input validation
  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (X) && ismatrix (X) && ndims (X) == 2))
    error ("rica: X must be a numeric matrix.");
  endif
  if (! (isreal (X) && all (isfinite (X(:)))))
    error ("rica: X must be real and finite.");
  endif

  [n, p] = size (X);

  if (! (isscalar (Q) && isnumeric (Q) && Q >= 1 && Q == fix (Q)))
    error ("rica: Q must be a positive integer.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("rica: Name/Value arguments must come in pairs.");
  endif

  ## Validate a supplied initial-weights matrix up front for a clear message
  for k = 1:2:numel (varargin)
    if (ischar (varargin{k}) && strcmpi (varargin{k}, "InitialTransformWeights"))
      W0 = varargin{k+1};
      if (! (isnumeric (W0) && isequal (size (W0), [p, Q])))
        error ("rica: 'InitialTransformWeights' must be a %d-by-%d matrix.", p, Q);
      endif
    endif
  endfor

  Mdl = ReconstructionICA (X, Q, varargin{:});

endfunction

%!demo
%! ## Learn two features from data with the default (random) start.
%! X = [1, 2, 3, 4; 2, 3, 4, 5; -1, 0, 1, 2; 3, 1, 4, 1; 0, 2, 1, 3];
%! Mdl = rica (X, 2, "IterationLimit", 200);
%! Z = transform (Mdl, X)

## The RICA objective is verified against MATLAB R2023b for
## X = reshape (mod ((1:60)*7, 13), 12, 5) - 6 with fixed initial weights.

%!shared X, W0
%! X = reshape (mod ((1:60)*7, 13), 12, 5) - 6;
%! W0 = reshape (mod ((1:15)*3, 7), 5, 3) - 3;

%!test
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!             "Standardize", false, "IterationLimit", 1000);
%! ## transform is X * TransformWeights, with unit-length weight columns
%! assert_equal (transform (Mdl, X), X * Mdl.TransformWeights, 1e-12);
%! assert_equal (sqrt (sum (Mdl.TransformWeights .^ 2, 1)), [1, 1, 1], 1e-10);
%! assert_equal (size (Mdl.TransformWeights), [5, 3]);

%!test
%! ## The solver reaches a low objective value (MATLAB's is 194.5156).
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!             "Standardize", false, "IterationLimit", 1000);
%! assert_equal (Mdl.FitInfo.Objective(end) < 195, true);

%!test
%! ## The objective value stored in FitInfo matches a direct evaluation.
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!             "Standardize", false, "IterationLimit", 1000);
%! W = Mdl.TransformWeights;
%! Z = X * W;
%! f = sum (sumsq (X * (W * W') - X)) + sum (sum (0.5 * log (cosh (2 * Z))));
%! assert_equal (Mdl.FitInfo.Objective(end), f, 1e-6);

%!test
%! ## The 'exp' contrast objective matches a direct evaluation.
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!             "Standardize", false, "IterationLimit", 200, ...
%!             "ContrastFcn", "exp");
%! W = Mdl.TransformWeights;
%! Z = X * W;
%! f = sum (sumsq (X * (W * W') - X)) + sum (sum (-exp (-Z .^ 2 / 2)));
%! assert_equal (Mdl.FitInfo.Objective(end), f, 1e-6);

%!test
%! ## The 'sqrt' contrast objective matches a direct evaluation, the smoothing
%! ## constant included.
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!             "Standardize", false, "IterationLimit", 200, ...
%!             "ContrastFcn", "sqrt");
%! W = Mdl.TransformWeights;
%! Z = X * W;
%! f = sum (sumsq (X * (W * W') - X)) + sum (sum (sqrt (Z .^ 2 + 1e-8)));
%! assert_equal (Mdl.FitInfo.Objective(end), f, 1e-6);

%!test
%! ## The smoothing constant is visible only where a projection is zero, so it
%! ## is pinned there: one column of X projected onto a weight orthogonal to it.
%! Xz = [1 0; -2 0; 3 0; 0 0];
%! Wz = [0; 1];
%! M = rica (Xz, 1, "InitialTransformWeights", Wz, "IterationLimit", 0, ...
%!           "Standardize", false, "ContrastFcn", "sqrt");
%! assert_equal (M.FitInfo.Objective(1), 14 + 4 * sqrt (1e-8), 1e-12);

%!test
%! ## Each contrast gives its own objective; they are not interchangeable.
%! args = {"InitialTransformWeights", W0, "Standardize", false, ...
%!         "IterationLimit", 0};
%! f1 = rica (X, 3, args{:}, "ContrastFcn", "logcosh").FitInfo.Objective(1);
%! f2 = rica (X, 3, args{:}, "ContrastFcn", "exp").FitInfo.Objective(1);
%! f3 = rica (X, 3, args{:}, "ContrastFcn", "sqrt").FitInfo.Objective(1);
%! assert_equal (f1 != f2 && f2 != f3 && f1 != f3, true);

%!test
%! ## Standardize centers and scales the data before transforming.
%! Mdl = rica (X, 2, "Standardize", true, "IterationLimit", 100, ...
%!             "InitialTransformWeights", W0(:,1:2));
%! assert_equal (size (Mdl.Mu), [5, 1]);
%! assert_equal (size (Mdl.Sigma), [5, 1]);

%!test
%! ## FitInfo carries the whole trajectory, not just its final value: two
%! ## columns of equal length, Iteration the 0-based index, Objective the
%! ## objective at the starting weights first and the solution last
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!             "Standardize", false, "IterationLimit", 1000);
%! it = Mdl.FitInfo.Iteration;
%! ob = Mdl.FitInfo.Objective;
%! assert_equal (columns (it), 1);
%! assert_equal (size (ob), size (it));
%! assert_equal (it, (0:numel (ob) - 1)');
%! assert_equal (all (diff (ob) <= 1e-10), true);
%! ## the first entry is the objective at the starting weights, which is what
%! ## a fit allowed no iterations at all reports
%! M0 = rica (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!            "Standardize", false, "IterationLimit", 0);
%! assert_equal (numel (M0.FitInfo.Objective), 1);
%! assert_equal (ob(1), M0.FitInfo.Objective(1), 1e-12);

%!test
%! ## a capped run stops with one entry per iteration plus the starting point
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!             "Standardize", false, "IterationLimit", 5);
%! assert_equal (numel (Mdl.FitInfo.Objective) <= 6, true);
%! assert_equal (Mdl.FitInfo.Iteration(1), 0);

## Test input validation
%!test
%! ## 'lbfgs' reaches the optimum the default solver reaches
%! args = {"InitialTransformWeights", W0, "Lambda", 1, ...
%!         "Standardize", false, "IterationLimit", 1000};
%! Mq = rica (X, 3, args{:});
%! Ml = rica (X, 3, args{:}, "Solver", "lbfgs");
%! assert_equal (Ml.FitInfo.Objective(end), Mq.FitInfo.Objective(end), 1e-6);

%!test
%! ## The trajectory starts at the initial weights whichever solver ran
%! args = {"InitialTransformWeights", W0, "IterationLimit", 20};
%! Mq = rica (X, 3, args{:});
%! Ml = rica (X, 3, args{:}, "Solver", "lbfgs");
%! assert_equal (Ml.FitInfo.Iteration(1), 0);
%! assert_equal (Ml.FitInfo.Objective(1), Mq.FitInfo.Objective(1), 1e-10);

%!test
%! ## The solver that ran is recorded
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "Solver", "lbfgs");
%! assert_equal (Mdl.ModelParameters.Solver, "lbfgs");

%!error<rica: 'Solver' must be 'quasinewton' or 'lbfgs'.> ...
%! rica (X, 3, "Solver", "bogus")
%!error<rica: 'Solver' must be 'quasinewton' or 'lbfgs'.> ...
%! rica (X, 3, "Solver", 5)

%!error<Invalid call to rica> rica (ones (5, 3))
%!error<rica: X must be a numeric matrix.> rica ({1, 2}, 1)
%!error<rica: X must be real and finite.> rica ([1, Inf; 2, 3], 1)
%!error<rica: Q must be a positive integer.> rica (ones (5, 3), 0)
%!error<rica: Q must be a positive integer.> rica (ones (5, 3), 1.5)
%!error<rica: 'InitialTransformWeights' must be a 5-by-2 matrix.> ...
%! rica (ones (4, 5), 2, "InitialTransformWeights", ones (3, 3))
%!test
%! ## The NonGaussianityIndicator pair reaches the model through rica.
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "IterationLimit", 100, ...
%!             "NonGaussianityIndicator", [1; -1; 1]);
%! assert_equal (Mdl.NonGaussianityIndicator, [1; -1; 1]);

%!error<rica: 'ContrastFcn' must be 'logcosh', 'exp', or 'sqrt'.> ...
%! rica (ones (4, 5), 2, "ContrastFcn", "bogus", "InitialTransformWeights", ones (5, 2))
%!error<rica: unknown parameter name 'bogus'.> ...
%! rica (ones (4, 5), 2, "bogus", 1)
