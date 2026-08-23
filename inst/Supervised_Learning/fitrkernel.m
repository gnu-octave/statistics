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
## @deftypefn  {statistics} {@var{Mdl} =} fitrkernel (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitrkernel (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{Mdl}, @var{FitInfo}] =} fitrkernel (@dots{})
##
## Fit a Gaussian kernel regression model.
##
## @code{@var{Mdl} = fitrkernel (@var{X}, @var{Y})} returns a
## @qcode{RegressionKernel} object fitted to the predictor data @var{X} and
## the continuous response @var{Y}, where @var{X} is an @math{NxP} numeric
## matrix and @var{Y} an @math{Nx1} numeric vector with as many rows as
## @var{X}.
##
## @code{@var{Mdl} = fitrkernel (@dots{}, @var{name}, @var{value})} passes
## the given @qcode{Name-Value} pairs to the model.  They are documented
## under @code{RegressionKernel}, and the ones most often wanted are
## @qcode{'Learner'}, @qcode{'Epsilon'},
## @qcode{'NumExpansionDimensions'}, @qcode{'KernelScale'},
## @qcode{'Lambda'} and @qcode{'BoxConstraint'}.
##
## @code{[@var{Mdl}, @var{FitInfo}] = fitrkernel (@dots{})} also returns a
## structure describing the optimization: the objective it reached, the
## gradient it left, and the tolerances it was given.
##
## @code{@var{Mdl} = fitrkernel (@dots{}, @var{cvopt}, @var{value})} returns a
## @code{RegressionPartitionedKernel}
## instead when one of @qcode{'CrossVal'}, @qcode{'KFold'},
## @qcode{'Holdout'}, @qcode{'Leaveout'} and @qcode{'CVPartition'} is
## given.  A cross-validated model describes no single fit, so
## @var{FitInfo} is not available beside it.
##
## @seealso{RegressionKernel, RegressionLinear, fitrlinear}
## @end deftypefn

function [Mdl, FitInfo] = fitrkernel (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitrkernel: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("fitrkernel: name-value arguments must be in pairs.");
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
        error ("%s: 'CrossVal' must be either 'on' or 'off'.", 'fitrkernel');
      endif
      crossval = crossval || strcmpi (val, 'on');
    endif
  endfor

  if (crossval)
    if (nargout > 1)
      error (strcat ("fitrkernel: a cross validated model has no", ...
                     " FitInfo to return; ask for the model alone."));
    endif
    Mdl = RegressionPartitionedKernel (X, Y, varargin{:});
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

  Mdl = RegressionKernel (X, Y, varargin{:});

  if (nargout > 1)
    FitInfo = fitInfo_ (Mdl);
  endif

endfunction

%!demo
%! ## Fit a Gaussian kernel regression to fuel consumption and read what the
%! ## optimization did.
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! [Mdl, FitInfo] = fitrkernel (X(ok,:), MPG(ok))

%!test
%! ## The driver returns a kernel regression model
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! Mdl = fitrkernel (X(ok,:), MPG(ok));
%! assert_equal (class (Mdl), 'RegressionKernel');
%! assert_equal (Mdl.Epsilon, 0.926612305411416, 1e-12);

%!test
%! ## The options reach the model
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! Mdl = fitrkernel (X(ok,:), MPG(ok), 'Learner', 'leastsquares', ...
%!                   'Standardize', true, 'NumExpansionDimensions', 64);
%! assert_equal (Mdl.Learner, 'leastsquares');
%! assert_equal (Mdl.NumExpansionDimensions, 64);
%! assert_equal (Mdl.Mu, mean (X(ok,:)), 1e-12);

%!test
%! ## The second output describes the optimization
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! [~, FitInfo] = fitrkernel (X(ok,:), MPG(ok));
%! assert_equal (FitInfo.Solver, 'LBFGS-fast');
%! assert_equal (FitInfo.LossFunction, 'epsiloninsensitive');

%!test
%! ## A cross-validation option returns a partitioned model instead
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! X = X(ok,:);
%! Y = MPG(ok);
%! CVMdl = fitrkernel (X, Y, 'KFold', 5);
%! assert_equal (class (CVMdl), 'RegressionPartitionedKernel');
%! assert_equal (CVMdl.KFold, 5);
%! assert_equal (numel (CVMdl.Trained), 5);

%!test
%! ## 'CrossVal' on gives the ten folds it defaults to, and 'off' the model
%! ## itself
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! X = X(ok,:);
%! Y = MPG(ok);
%! CVMdl = fitrkernel (X, Y, 'CrossVal', 'on');
%! assert_equal (CVMdl.KFold, 10);
%! assert_equal (class (fitrkernel (X, Y, 'CrossVal', 'off')), ...
%!                     'RegressionKernel');

%!error<fitrkernel: a cross validated model has no FitInfo to return; ask for the model alone.> ...
%! [Mdl, FitInfo] = fitrkernel (ones (10, 2), ones (10, 1), 'KFold', 3);

## Test input validation
%!error<fitrkernel: too few input arguments.> fitrkernel (ones (5, 2))
%!error<fitrkernel: name-value arguments must be in pairs.> ...
%! fitrkernel (ones (10, 2), ones (10, 1), 'Learner')
%!error<RegressionKernel: 'Learner' must be either 'svm' or 'leastsquares'.> ...
%! fitrkernel (ones (10, 2), ones (10, 1), 'Learner', 'logistic')
