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
##       + sum (sum (0.5 * log (cosh (2 * @var{X} * @var{W}))))
## @end example
##
## @noindent
## over the transformation weights @var{W}, combining a reconstruction cost with
## a @qcode{'logcosh'} sparsity contrast.
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
## The sparsity contrast.  Only @qcode{'logcosh'} — that is
## @math{0.5 \log (\cosh (2 z))} — is supported.
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
%! assert_equal (Mdl.FitInfo.Objective < 195, true);

%!test
%! ## The objective value stored in FitInfo matches a direct evaluation.
%! Mdl = rica (X, 3, "InitialTransformWeights", W0, "Lambda", 1, ...
%!             "Standardize", false, "IterationLimit", 1000);
%! W = Mdl.TransformWeights;
%! Z = X * W;
%! f = sum (sumsq (X * (W * W') - X)) + sum (sum (0.5 * log (cosh (2 * Z))));
%! assert_equal (Mdl.FitInfo.Objective, f, 1e-6);

%!test
%! ## Standardize centers and scales the data before transforming.
%! Mdl = rica (X, 2, "Standardize", true, "IterationLimit", 100, ...
%!             "InitialTransformWeights", W0(:,1:2));
%! assert_equal (size (Mdl.Mu), [1, 5]);
%! assert_equal (size (Mdl.Sigma), [1, 5]);

## Test input validation
%!error<Invalid call to rica> rica (ones (5, 3))
%!error<rica: X must be a numeric matrix.> rica ({1, 2}, 1)
%!error<rica: X must be real and finite.> rica ([1, Inf; 2, 3], 1)
%!error<rica: Q must be a positive integer.> rica (ones (5, 3), 0)
%!error<rica: Q must be a positive integer.> rica (ones (5, 3), 1.5)
%!error<rica: 'InitialTransformWeights' must be a 5-by-2 matrix.> ...
%! rica (ones (4, 5), 2, "InitialTransformWeights", ones (3, 3))
%!error<rica: only the 'logcosh' ContrastFcn is supported.> ...
%! rica (ones (4, 5), 2, "ContrastFcn", "exp", "InitialTransformWeights", ones (5, 2))
%!error<rica: unknown parameter name 'bogus'.> ...
%! rica (ones (4, 5), 2, "bogus", 1)
