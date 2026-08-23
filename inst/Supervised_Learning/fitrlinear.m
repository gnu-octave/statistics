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
## @deftypefn  {statistics} {@var{Mdl} =} fitrlinear (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitrlinear (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{Mdl}, @var{FitInfo}] =} fitrlinear (@dots{})
##
## Fit a linear regression model.
##
## @code{@var{Mdl} = fitrlinear (@var{X}, @var{Y})} returns a
## @qcode{RegressionLinear} object fitted to the predictor data @var{X} and
## the continuous response @var{Y}, where @var{X} is an @math{NxP} numeric
## matrix and @var{Y} an @math{Nx1} numeric vector with as many rows as
## @var{X}.
##
## @code{@var{Mdl} = fitrlinear (@dots{}, @var{name}, @var{value})} passes
## the given @qcode{Name-Value} pairs to the model.  They are documented
## under @code{RegressionLinear}, and the ones most often wanted are
## @qcode{'Learner'}, @qcode{'Epsilon'}, @qcode{'Regularization'},
## @qcode{'Lambda'} and @qcode{'Solver'}.
##
## @code{[@var{Mdl}, @var{FitInfo}] = fitrlinear (@dots{})} also returns a
## structure describing the optimization: what it converged to, how far it
## got, and which tolerance stopped it.  Its fields follow the solver, so a
## dual fit reports the dual variables and a mini-batch fit the batch it
## stopped on.
##
## @code{@var{Mdl} = fitrlinear (@dots{}, @var{cvopt}, @var{value})} returns a
## @code{RegressionPartitionedLinear}
## instead when one of @qcode{'CrossVal'}, @qcode{'KFold'},
## @qcode{'Holdout'}, @qcode{'Leaveout'} and @qcode{'CVPartition'} is
## given.  A cross-validated model describes no single fit, so
## @var{FitInfo} is not available beside it.
##
## @seealso{RegressionLinear, RegressionKernel, fitrkernel}
## @end deftypefn

function [Mdl, FitInfo] = fitrlinear (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitrlinear: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("fitrlinear: name-value arguments must be in pairs.");
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
        error ("%s: 'CrossVal' must be either 'on' or 'off'.", 'fitrlinear');
      endif
      crossval = crossval || strcmpi (val, 'on');
    endif
  endfor

  if (crossval)
    if (nargout > 1)
      error (strcat ("fitrlinear: a cross validated model has no", ...
                     " FitInfo to return; ask for the model alone."));
    endif
    Mdl = RegressionPartitionedLinear (X, Y, varargin{:});
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

  Mdl = RegressionLinear (X, Y, varargin{:});

  if (nargout > 1)
    FitInfo = fitInfo_ (Mdl);
  endif

endfunction

%!demo
%! ## Fit a linear regression to fuel consumption and read what the
%! ## optimization did.
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! [Mdl, FitInfo] = fitrlinear (X(ok,:), MPG(ok), 'Learner', 'leastsquares')

%!test
%! ## The driver returns what the class constructor returns
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! M1 = fitrlinear (X(ok,:), MPG(ok));
%! M2 = RegressionLinear (X(ok,:), MPG(ok));
%! assert_equal (class (M1), 'RegressionLinear');
%! assert_equal (M1.Beta, M2.Beta);
%! assert_equal (M1.Epsilon, M2.Epsilon);

%!test
%! ## The options reach the model
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! Mdl = fitrlinear (X(ok,:), MPG(ok), 'Learner', 'leastsquares', ...
%!                   'Lambda', 0.02, 'ResponseName', 'mpg');
%! assert_equal (Mdl.Learner, 'leastsquares');
%! assert_equal (Mdl.Lambda, 0.02);
%! assert_equal (Mdl.ResponseName, 'mpg');

%!test
%! ## The second output describes the optimization
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! [~, FitInfo] = fitrlinear (X(ok,:), MPG(ok));
%! assert_equal (FitInfo.Solver, {'bfgs'});
%! assert_equal (FitInfo.Lambda, 1 / 93, 1e-15);
%! assert_equal (isfield (FitInfo, 'Objective'), true);

%!test
%! ## A cross-validation option returns a partitioned model instead
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! X = X(ok,:);
%! Y = MPG(ok);
%! CVMdl = fitrlinear (X, Y, 'KFold', 5);
%! assert_equal (class (CVMdl), 'RegressionPartitionedLinear');
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
%! CVMdl = fitrlinear (X, Y, 'CrossVal', 'on');
%! assert_equal (CVMdl.KFold, 10);
%! assert_equal (class (fitrlinear (X, Y, 'CrossVal', 'off')), ...
%!                     'RegressionLinear');

%!error<fitrlinear: a cross validated model has no FitInfo to return; ask for the model alone.> ...
%! [Mdl, FitInfo] = fitrlinear (ones (10, 2), ones (10, 1), 'KFold', 3);

## Test input validation
%!error<fitrlinear: too few input arguments.> fitrlinear (ones (5, 2))
%!error<fitrlinear: name-value arguments must be in pairs.> ...
%! fitrlinear (ones (10, 2), ones (10, 1), 'Learner')
%!error<RegressionLinear: 'Learner' must be either 'svm' or 'leastsquares'.> ...
%! fitrlinear (ones (10, 2), ones (10, 1), 'Learner', 'logistic')
