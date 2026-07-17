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
## @deftypefn  {statistics} {@var{Mdl} =} sparsefilt (@var{X}, @var{Q})
## @deftypefnx {statistics} {@var{Mdl} =} sparsefilt (@var{X}, @var{Q}, @var{Name}, @var{Value})
##
## Sparse filtering for feature extraction.
##
## @code{@var{Mdl} = sparsefilt (@var{X}, @var{Q})} learns @var{Q} features from
## the @math{N * P} data matrix @var{X} (rows are observations, columns are
## predictors) and returns a @qcode{SparseFiltering} object @var{Mdl}.  Apply the
## learned transformation to data with @code{transform (@var{Mdl}, @var{X})}.
##
## The @math{N * @var{Q}} features returned by @code{transform} are the
## soft-absolute activations @code{sqrt ((@var{X} * @var{W}) .^ 2 + 1e-8)},
## normalized first across observations (each feature) and then across features
## (each observation).  The @math{P * @var{Q}} weight matrix @var{W} minimizes
## the sum of those features plus an L2 penalty @code{@var{Lambda} * ||@var{W}
## ||_F^2}, driving the features to be sparse.
##
## Name/Value pairs:
##
## @table @asis
## @item @qcode{'IterationLimit'}
## Maximum number of iterations (default 1000).
##
## @item @qcode{'Lambda'}
## Weight of the L2 penalty on the transform weights (default 1).
##
## @item @qcode{'Standardize'}
## Logical; center and scale each predictor before fitting (default
## @qcode{false}).
##
## @item @qcode{'InitialTransformWeights'}
## A @math{P * @var{Q}} initial value for the weights.  The default is random.
##
## @item @qcode{'GradientTolerance'}, @qcode{'StepTolerance'}
## Accepted for MATLAB compatibility and recorded in @qcode{ModelParameters};
## the built-in solver iterates up to @qcode{'IterationLimit'} with tight
## internal tolerances.
## @end table
##
## @subheading Note on reproducibility
##
## The sparse filtering objective is not convex and is minimized by a
## quasi-Newton solver, so the learned weights depend on the starting point and
## the solver.  Different runs (or different software, including MATLAB) may
## return different weights that nonetheless describe an equally valid feature
## transformation.  Fix @qcode{'InitialTransformWeights'} for a reproducible
## result.
##
## @seealso{SparseFiltering, rica, pca}
## @end deftypefn

function Mdl = sparsefilt (X, Q, varargin)

  ## Input validation
  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (X) && ismatrix (X) && ndims (X) == 2))
    error ("sparsefilt: X must be a numeric matrix.");
  endif
  if (! (isreal (X) && all (isfinite (X(:)))))
    error ("sparsefilt: X must be real and finite.");
  endif

  [n, p] = size (X);

  if (! (isscalar (Q) && isnumeric (Q) && Q >= 1 && Q == fix (Q)))
    error ("sparsefilt: Q must be a positive integer.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("sparsefilt: Name/Value arguments must come in pairs.");
  endif

  for k = 1:2:numel (varargin)
    if (ischar (varargin{k}) && strcmpi (varargin{k}, "InitialTransformWeights"))
      W0 = varargin{k+1};
      if (! (isnumeric (W0) && isequal (size (W0), [p, Q])))
        error ("sparsefilt: 'InitialTransformWeights' must be a %d-by-%d matrix.", ...
               p, Q);
      endif
    endif
  endfor

  Mdl = SparseFiltering (X, Q, varargin{:});

endfunction

%!demo
%! ## Learn two sparse features from data with the default (random) start.
%! X = [1, 2, 3, 4; 2, 3, 4, 5; -1, 0, 1, 2; 3, 1, 4, 1; 0, 2, 1, 3];
%! Mdl = sparsefilt (X, 2, "IterationLimit", 200);
%! Z = transform (Mdl, X)

## The sparse filtering objective is verified against MATLAB R2023b for
## X = reshape (mod ((1:60)*7, 13), 12, 5) - 6 with fixed initial weights.

%!shared X, W0
%! X = reshape (mod ((1:60)*7, 13), 12, 5) - 6;
%! W0 = reshape (mod ((1:15)*3, 7), 5, 3) - 3;

%!test
%! Mdl = sparsefilt (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!                   "Standardize", false, "IterationLimit", 1000);
%! Z = transform (Mdl, X);
%! ## features are in [0, 1] and each observation (row) is essentially unit length
%! assert_equal (size (Z), [12, 3]);
%! assert_equal (all (Z(:) >= 0) && all (Z(:) <= 1 + 1e-12), true);
%! ## rows are unit length up to the 1e-8 normalization regularizer
%! assert_equal (sqrt (sum (Z .^ 2, 2)), ones (12, 1), 1e-3);

%!test
%! ## The objective at the initial weights matches MATLAB (73.171833).
%! A = X * W0;
%! F = sqrt (A .^ 2 + 1e-8);
%! F = F ./ sqrt (sum (F .^ 2, 1) + 1e-8);
%! F = F ./ sqrt (sum (F .^ 2, 2) + 1e-8);
%! f0 = sum (F(:)) + sum (W0(:) .^ 2);
%! assert_equal (f0, 73.171833042485730, 1e-9);

%!test
%! ## The solver drives the objective well below its value at the start
%! ## (73.17) -- reaching a good local minimum of the sparse filtering cost.
%! Mdl = sparsefilt (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!                   "Standardize", false, "IterationLimit", 1000);
%! assert_equal (Mdl.FitInfo.Objective < 20, true);

%!test
%! ## Lambda is an L2 penalty on the weights: FitInfo.Objective is
%! ## sum(features) + Lambda * ||W||^2 at the solution.
%! Mdl = sparsefilt (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!                   "Standardize", false, "IterationLimit", 1000);
%! W = Mdl.TransformWeights;
%! Z = transform (Mdl, X);
%! assert_equal (Mdl.FitInfo.Objective, sum (Z(:)) + sum (W(:) .^ 2), 1e-6);

## Test input validation
%!error<Invalid call to sparsefilt> sparsefilt (ones (5, 3))
%!error<sparsefilt: X must be a numeric matrix.> sparsefilt ({1, 2}, 1)
%!error<sparsefilt: X must be real and finite.> sparsefilt ([1, Inf; 2, 3], 1)
%!error<sparsefilt: Q must be a positive integer.> sparsefilt (ones (5, 3), 0)
%!error<sparsefilt: Q must be a positive integer.> sparsefilt (ones (5, 3), 1.5)
%!error<sparsefilt: 'InitialTransformWeights' must be a 5-by-2 matrix.> ...
%! sparsefilt (ones (4, 5), 2, "InitialTransformWeights", ones (3, 3))
%!error<sparsefilt: unknown parameter name 'bogus'.> ...
%! sparsefilt (ones (4, 5), 2, "bogus", 1)
