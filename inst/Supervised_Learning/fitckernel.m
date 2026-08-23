## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{Mdl} =} fitckernel (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitckernel (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{Mdl}, @var{FitInfo}] =} fitckernel (@dots{})
##
## Fit a Gaussian kernel binary classifier.
##
## @code{@var{Mdl} = fitckernel (@var{X}, @var{Y})} returns a
## @qcode{ClassificationKernel} object fitted to the predictor data @var{X}
## and the two class response @var{Y}, where @var{X} is an @math{NxP}
## numeric matrix and @var{Y} has as many rows as @var{X}.
##
## @code{@var{Mdl} = fitckernel (@dots{}, @var{name}, @var{value})} passes
## the given @qcode{Name-Value} pairs to the model.  They are documented
## under @code{ClassificationKernel}, and the ones most often wanted are
## @qcode{'Learner'}, @qcode{'NumExpansionDimensions'},
## @qcode{'KernelScale'}, @qcode{'Lambda'}, @qcode{'BoxConstraint'} and
## @qcode{'Standardize'}.
##
## @code{[@var{Mdl}, @var{FitInfo}] = fitckernel (@dots{})} also returns a
## structure describing the optimization: the objective it reached, the
## gradient it left, and the tolerances it was given.
##
## @code{@var{Mdl} = fitckernel (@dots{}, @var{cvopt}, @var{value})} returns a
## @code{ClassificationPartitionedKernel}
## instead when one of @qcode{'CrossVal'}, @qcode{'KFold'},
## @qcode{'Holdout'}, @qcode{'Leaveout'} and @qcode{'CVPartition'} is
## given.  A cross-validated model describes no single fit, so
## @var{FitInfo} is not available beside it.
##
## @seealso{ClassificationKernel, ClassificationLinear, fitclinear}
## @end deftypefn

function [Mdl, FitInfo] = fitckernel (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitckernel: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("fitckernel: name-value arguments must be in pairs.");
  endif

  ## A cross-validation option asks for a partitioned model rather than a
  ## fitted one, and the two are different classes with different methods.
  ## MATLAB refuses a second output there, having no single fit to describe,
  ## and so does this.
  ## 'CrossVal' is the one that carries a value saying whether to cross
  ## validate at all; the other four ask for it by being present.
  cvNames = {'kfold', 'holdout', 'leaveout', 'cvpartition'};
  crossval = false;
  for k = 1:2:numel (varargin)
    if (! ischar (varargin{k}))
      continue;
    endif
    if (any (strcmpi (varargin{k}, cvNames)))
      crossval = true;
    elseif (strcmpi (varargin{k}, 'crossval'))
      val = varargin{k+1};
      if (! (ischar (val) && any (strcmpi (val, {'on', 'off'}))))
        error ("%s: 'CrossVal' must be either 'on' or 'off'.", 'fitckernel');
      endif
      crossval = crossval || strcmpi (val, 'on');
    endif
  endfor

  if (crossval)
    if (nargout > 1)
      error (strcat ("fitckernel: a cross validated model has no", ...
                     " FitInfo to return; ask for the model alone."));
    endif
    Mdl = ClassificationPartitionedKernel (X, Y, varargin{:});
    return;
  endif

  ## 'CrossVal', 'off' has said its piece and the learner does not take it.
  keep = true (1, numel (varargin));
  for k = 1:2:numel (varargin)
    if (ischar (varargin{k}) && strcmpi (varargin{k}, 'crossval'))
      keep(k:k+1) = false;
    endif
  endfor
  varargin = varargin(keep);

  Mdl = ClassificationKernel (X, Y, varargin{:});

  if (nargout > 1)
    FitInfo = fitInfo_ (Mdl);
  endif

endfunction

%!demo
%! ## Fit a Gaussian kernel classifier to the two overlapping iris species
%! ## and read what the optimization did.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! [Mdl, FitInfo] = fitckernel (X, Y)

%!test
%! ## The driver returns a kernel classifier
%! load fisheriris
%! Mdl = fitckernel (meas(51:end,:), species(51:end));
%! assert_equal (class (Mdl), 'ClassificationKernel');
%! assert_equal (Mdl.NumExpansionDimensions, 128);

%!test
%! ## The options reach the model
%! load fisheriris
%! Mdl = fitckernel (meas(51:end,:), species(51:end), ...
%!                   'Learner', 'logistic', 'KernelScale', 2, ...
%!                   'NumExpansionDimensions', 64);
%! assert_equal (Mdl.Learner, 'logistic');
%! assert_equal (Mdl.KernelScale, 2);
%! assert_equal (Mdl.NumExpansionDimensions, 64);

%!test
%! ## The second output describes the optimization
%! load fisheriris
%! [~, FitInfo] = fitckernel (meas(51:end,:), species(51:end));
%! assert_equal (FitInfo.Solver, 'LBFGS-fast');
%! assert_equal (FitInfo.LossFunction, 'hinge');
%! assert_equal (FitInfo.Lambda, 0.01);

%!test
%! ## A cross-validation option returns a partitioned model instead
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = fitckernel (X, Y, 'KFold', 5);
%! assert_equal (class (CVMdl), 'ClassificationPartitionedKernel');
%! assert_equal (CVMdl.KFold, 5);
%! assert_equal (numel (CVMdl.Trained), 5);

%!test
%! ## 'CrossVal' on gives the ten folds it defaults to, and 'off' the model
%! ## itself
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = fitckernel (X, Y, 'CrossVal', 'on');
%! assert_equal (CVMdl.KFold, 10);
%! assert_equal (class (fitckernel (X, Y, 'CrossVal', 'off')), ...
%!                     'ClassificationKernel');

%!error<fitckernel: a cross validated model has no FitInfo to return; ask for the model alone.> ...
%! [Mdl, FitInfo] = fitckernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'KFold', 3);

## Test input validation
%!error<fitckernel: too few input arguments.> fitckernel (ones (5, 2))
%!error<fitckernel: name-value arguments must be in pairs.> ...
%! fitckernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner')
%!error<ClassificationKernel: 'Learner' must be either 'svm' or 'logistic'.> ...
%! fitckernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner', 'tree')
