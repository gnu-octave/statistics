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
## @deftypefn  {statistics} {@var{Mdl} =} fitcecoc (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcecoc (@dots{}, @var{name}, @var{value})
##
## Fit a multiclass model from binary learners.
##
## @code{@var{Mdl} = fitcecoc (@var{X}, @var{Y})} turns the multiclass
## problem in @var{X} and @var{Y} into a set of two class problems, fits one
## binary learner to each, and returns them as a @code{ClassificationECOC}
## object.  A coding matrix says which classes each learner calls +1, which
## it calls -1, and which sit it out; a new observation is sent to every
## learner and given the class whose row of that matrix its scores match most
## closely.
##
## @itemize
## @item
## @var{X} must be a @math{NxP} numeric matrix of predictor data.
## @item
## @var{Y} must be a @math{Nx1} vector of class labels, of any type
## @code{ClassNames} accepts, with at least two distinct values.
## @end itemize
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{'Learners'} @tab @tab The binary learner, either a name,
## @qcode{'svm'} (default), @qcode{'tree'}, @qcode{'knn'},
## @qcode{'naivebayes'}, @qcode{'discriminant'}, @qcode{'linear'} or
## @qcode{'kernel'}, or a template from @code{templateSVM} and its siblings,
## which also carries the options that learner is to be fitted with.
##
## @item @qcode{'Coding'} @tab @tab The coding design, either a name
## @code{designecoc} accepts, @qcode{'onevsone'} by default, or a coding
## matrix given outright, which sets @code{CodingName} to @qcode{'custom'}.
##
## @item @qcode{'BinaryLoss'} @tab @tab The loss the binary scores are read
## with.  The default follows the learner: @qcode{'hinge'} for one scoring on
## @math{(-Inf,+Inf)} and @qcode{'quadratic'} for one scoring on @math{[0,1]}.
##
## @item @qcode{'ClassNames'} @tab @tab The classes to fit, and the order
## the rows of the coding matrix, @code{Prior} and @code{Cost} take them in.
##
## @item @qcode{'Cost'} @tab @tab A @math{KxK} matrix of misclassification
## costs.  The default is @code{1 - eye (K)}.
##
## @item @qcode{'Prior'} @tab @tab @qcode{'empirical'} (default),
## @qcode{'uniform'}, or a vector with one element per class.
##
## @item @qcode{'Weights'} @tab @tab A nonnegative numeric vector with one
## element per observation.  The default is uniform.
##
## @item @qcode{'PredictorNames'} @tab @tab A cellstr of predictor names.
##
## @item @qcode{'ResponseName'} @tab @tab The name of the response variable.
##
## @item @qcode{'ScoreTransform'} @tab @tab A transform applied to the
## returned scores.  The default is @qcode{'none'}.
## @end multitable
##
## A linear or kernel learner carries no training data of its own, so
## @code{fitcecoc} returns a @code{CompactClassificationECOC} for those two
## and a @code{ClassificationECOC} for the other five, which is what R2024a
## does.
##
## @code{'FitPosterior'} is refused rather than quietly ignored: it installs
## a fitted score transform on each binary learner, which needs the posterior
## fitting this package does not have yet.
##
## @seealso{ClassificationECOC, CompactClassificationECOC, designecoc,
## templateSVM, templateTree}
## @end deftypefn

function Mdl = fitcecoc (X, Y, varargin)

  ## Input validation
  if (nargin < 2)
    error ("fitcecoc: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("fitcecoc: name-value arguments must be in pairs.");
  endif
  if (rows (X) != rows (Y))
    error ("fitcecoc: number of rows in X and Y must be equal.");
  endif

  ## A cross-validation option makes the fit a cross-validated one: the
  ## model is fitted on all the data first, then partitioned, which is how
  ## every other fitc* here reads these four.
  cvnames = {'kfold', 'holdout', 'leaveout', 'cvpartition'};
  iscv = cellfun (@(a) ischar (a) && isrow (a) ...
                       && any (strcmpi (a, cvnames)), varargin(1:2:end));
  cvargs = {};
  if (any (iscv))
    if (sum (iscv) > 1)
      error (strcat ("fitcecoc: specify only one of 'KFold', 'Holdout',", ...
                     " 'Leaveout' and 'CVPartition'."));
    endif
    at = 2 * find (iscv) - 1;
    cvargs = varargin(at:at+1);
    varargin(at:at+1) = [];
  endif

  Mdl = ClassificationECOC (X, Y, varargin{:});

  if (! isempty (cvargs))
    Mdl = crossval (Mdl, cvargs{:});
    return;
  endif

  ## A learner that keeps no training data leaves nothing for the full class
  ## to hold, so the fit gives back the compact one.  Measured on R2024a.
  if (any (strcmpi (Mdl.ModelParameters.BinaryLearners.Method, ...
                    {'Linear', 'Kernel'})))
    Mdl = compact (Mdl);
  endif

endfunction

## Tests
%!test  # MATLAB parity: the default fit and what it reports
%! load fisheriris
%! Mdl = fitcecoc (meas, species);
%! assert_equal (class (Mdl), 'ClassificationECOC');
%! assert_equal (numel (properties (Mdl)), 22);
%! assert_equal (Mdl.CodingName, 'onevsone');
%! assert_equal (Mdl.BinaryLoss, 'hinge');
%! assert_equal (Mdl.CodingMatrix, [1, 1, 0; -1, 0, 1; 0, -1, -1]);
%! assert_equal (Mdl.LearnerWeights, [2/3, 2/3, 2/3], 1e-12);
%! assert_equal (Mdl.NumObservations, 150);

%!test  # MATLAB parity: the loss of the default fit
%! ## Our binary SVM is LIBSVM where MATLAB's is SMO, so the scores differ in
%! ## the fourth digit on the pair that is not separable; the labels they
%! ## lead to do not.
%! load fisheriris
%! Mdl = fitcecoc (meas, species);
%! assert_equal (resubLoss (Mdl), 0.0066666666666667, 1e-12);

%!test  # MATLAB parity: a tree code reproduces the decoding exactly
%! ## Our trees match MATLAB's on this fixture, so the whole path can be
%! ## compared and not only the labels.  Measured on R2024a.
%! load fisheriris
%! Mdl = fitcecoc (meas, species, 'Learners', 'tree');
%! assert_equal (Mdl.BinaryLoss, 'quadratic');
%! [~, NegLoss] = predict (Mdl, meas([1, 20, 51, 70, 101, 130], :));
%! assert_equal (NegLoss, ...
%!               [0, -1, -2; ...
%!                0, -1, -2; ...
%!                -2, 0, -1; ...
%!                -2, 0, -1; ...
%!                -2, -0.956994328922495, -0.000472589792060491; ...
%!                -2, -0.444444444444445, -0.111111111111111], 1e-12);

%!test  # MATLAB parity: the default binary loss follows the learner
%! load fisheriris
%! for p = {{'svm', 'hinge'}, {'tree', 'quadratic'}, {'knn', 'quadratic'}, ...
%!          {'naivebayes', 'quadratic'}, {'discriminant', 'quadratic'}}
%!   Mdl = fitcecoc (meas, species, 'Learners', p{1}{1});
%!   assert_equal (Mdl.BinaryLoss, p{1}{2});
%! endfor

%!test  # MATLAB parity: a learner keeping no data gives a compact model
%! load fisheriris
%! assert_equal (class (fitcecoc (meas, species, 'Learners', 'linear')), ...
%!               'CompactClassificationECOC');
%! assert_equal (class (fitcecoc (meas, species, 'Learners', 'kernel')), ...
%!               'CompactClassificationECOC');

%!test  # a template carries the options its learner is fitted with
%! load fisheriris
%! Mdl = fitcecoc (meas, species, 'Learners', templateTree ('MaxNumSplits', 2));
%! assert_equal (Mdl.BinaryLearners{1}.ModelParameters.MaxSplits, 2);

%!test  # MATLAB parity: the coding design is taken by name
%! load fisheriris
%! Mdl = fitcecoc (meas, species, 'Coding', 'onevsall');
%! assert_equal (Mdl.CodingName, 'onevsall');
%! assert_equal (Mdl.CodingMatrix, 2 * eye (3) - 1);
%! assert_equal (numel (Mdl.BinaryLearners), 3);

%!test  # MATLAB parity: a coding matrix given outright is 'custom'
%! load fisheriris
%! M = [1, 0, -1; -1, 1, 0; 0, -1, 1];
%! Mdl = fitcecoc (meas, species, 'Coding', M);
%! assert_equal (Mdl.CodingName, 'custom');
%! assert_equal (Mdl.CodingMatrix, M);

%!test  # what an observation was to each learner is its class's row
%! load fisheriris
%! Mdl = fitcecoc (meas, species);
%! assert_equal (Mdl.BinaryY(1,:), [1, 1, 0]);
%! assert_equal (Mdl.BinaryY(51,:), [-1, 0, 1]);
%! assert_equal (Mdl.BinaryY(101,:), [0, -1, -1]);

%!test  # MATLAB parity: a cross-validation option gives a partitioned model
%! load fisheriris
%! CV = fitcecoc (meas, species, 'KFold', 5);
%! assert_equal (class (CV), 'ClassificationPartitionedECOC');
%! assert_equal (CV.KFold, 5);
%! assert_equal (class (fitcecoc (meas, species, 'Holdout', 0.3)), ...
%!               'ClassificationPartitionedECOC');

## Test input validation
%!error<fitcecoc: specify only one of 'KFold', 'Holdout', 'Leaveout' and 'CVPartition'.> ...
%! fitcecoc (ones (8, 2), [1; 2; 1; 2; 1; 2; 1; 2], 'KFold', 2, 'Holdout', 0.3)
%!error<fitcecoc: too few input arguments.> fitcecoc (ones (4, 2))
%!error<fitcecoc: name-value arguments must be in pairs.> ...
%! fitcecoc (ones (4, 2), [1; 2; 1; 2], 'Coding')
%!error<fitcecoc: number of rows in X and Y must be equal.> ...
%! fitcecoc (ones (4, 2), [1; 2; 1])
%!error<ClassificationECOC: 'nosuch' is not a binary learner.> ...
%! fitcecoc (ones (4, 2), [1; 2; 1; 2], 'Learners', 'nosuch')
%!error<ClassificationECOC: 'FitPosterior' is not implemented> ...
%! fitcecoc (ones (4, 2), [1; 2; 1; 2], 'FitPosterior', true)
%!error<ClassificationECOC: invalid parameter name> ...
%! fitcecoc (ones (4, 2), [1; 2; 1; 2], 'NoSuch', 1)
