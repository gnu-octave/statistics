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
## @deftypefn {Private Function} {@var{MP} =} linearModelParams (@var{P}, @var{Solver}, @var{LambdaIn}, @var{EpsilonIn}, @var{Beta0}, @var{Bias0}, @var{Verbose}, @var{Type})
##
## The @qcode{ModelParameters} structure of a linear model.
##
## Every option is reported, but only the ones the chosen solver can act on
## carry a value: the rest come back empty.  That is MATLAB's behaviour and
## it is more informative than reporting a default that nothing used.
## Measured on R2024a, one fit per solver:
##
## @multitable @columnfractions 0.30 0.14 0.14 0.14 0.14 0.14
## @headitem field @tab bfgs @tab lbfgs @tab sparsa @tab sgd @tab dual
## @item @qcode{BatchIndex} @tab @tab @tab @tab 0 @tab
## @item @qcode{BatchLimit} @tab @tab @tab @tab 0 @tab
## @item @qcode{BatchSize} @tab @tab @tab @tab 10 @tab
## @item @qcode{DeltaGradientTolerance} @tab @tab @tab @tab @tab yes
## @item @qcode{GradientTolerance} @tab yes @tab yes @tab yes @tab 0 @tab 0
## @item @qcode{HessianHistorySize} @tab @tab 15 @tab @tab @tab
## @item @qcode{IterationLimit} @tab yes @tab yes @tab yes @tab @tab
## @item @qcode{LearnRate} @tab @tab @tab @tab yes @tab
## @item @qcode{LineSearch} @tab yes @tab yes @tab @tab @tab
## @item @qcode{NumCheckConvergence} @tab @tab @tab @tab 0 @tab yes
## @item @qcode{OptimizeLearnRate} @tab @tab @tab @tab yes @tab
## @item @qcode{PassLimit} @tab @tab @tab @tab yes @tab yes
## @item @qcode{TruncationPeriod} @tab @tab @tab @tab @tab
## @end multitable
##
## Two of those are worth noticing.  @qcode{HessianHistorySize} is empty for
## @qcode{'bfgs'} and 15 for @qcode{'lbfgs'}, a full quasi-Newton method
## having no limited memory to size.  And @qcode{GradientTolerance} is
## reported as 0 by the two solvers that do not run a gradient test, which
## by the rule that a tolerance of 0 switches its test off is a true
## statement rather than a missing one.
##
## @qcode{LineSearch} reports @qcode{'strongwolfe'}, which is what
## @code{__lbfgs__} performs.  MATLAB reports @qcode{'weakwolfe'} there.
## The field names the line search that actually ran, so reporting MATLAB's
## word for ours would be false.
##
## @end deftypefn

function MP = linearModelParams (P, Solver, LambdaIn, EpsilonIn, Beta0, ...
                                 Bias0, Verbose, Type)

  quasi = any (ismember (Solver, {'bfgs', 'lbfgs'}));
  batch = any (ismember (Solver, {'sgd', 'asgd'}));
  dual = any (strcmp (Solver, 'dual'));
  iter = quasi || any (strcmp (Solver, 'sparsa'));

  MP = struct ();
  MP.BatchIndex = pick (batch, 0);
  MP.BatchLimit = pick (batch, ifelse (isempty (P.BatchLimit), 0, ...
                                       P.BatchLimit));
  MP.BatchSize = pick (batch, P.BatchSize);
  MP.BetaTolerance = P.BetaTolerance;
  MP.DeltaGradientTolerance = pick (dual, P.DeltaGradientTolerance);
  MP.Epsilon = EpsilonIn;
  MP.FitBias = P.FitBias;
  ## The two solvers that run no gradient test report a tolerance of zero,
  ## which is how a switched-off test is spelled throughout this layer.
  MP.GradientTolerance = P.GradientTolerance;
  if (batch || dual)
    MP.GradientTolerance = 0;
  endif
  MP.HessianHistorySize = pick (any (strcmp (Solver, 'lbfgs')), ...
                                P.HessianHistorySize);
  MP.InitialBeta = Beta0;
  MP.InitialBias = Bias0;
  MP.IterationLimit = pick (iter, P.IterationLimit);
  MP.Learner = P.Learner;
  MP.Lambda = LambdaIn;
  MP.LearnRate = pick (batch, P.LearnRate);
  MP.LineSearch = pick (quasi, 'strongwolfe');
  MP.LossFunction = P.LossFunction;
  MP.NumCheckConvergence = pick (batch || dual, P.NumCheckConvergence);
  if (batch)
    MP.NumCheckConvergence = 0;
  endif
  MP.OptimizeLearnRate = pick (batch, P.OptimizeLearnRate);
  MP.PassLimit = pick (batch || dual, P.PassLimit);
  MP.PostFitBias = P.PostFitBias;
  MP.Regularization = P.Regularization;
  ## Solver arrives as a cell array of names already; wrapping it again
  ## would nest it, which struct () would have unwrapped but a direct
  ## assignment does not.
  MP.Solver = Solver;
  ## We take no RandomStream, so there is never one to report.
  MP.Stream = [];
  MP.TruncationPeriod = pick (batch && strcmp (P.Regularization, 'lasso'), ...
                              P.TruncationPeriod);
  MP.VerbosityLevel = Verbose;
  MP.Method = 'Linear';
  MP.Type = Type;

endfunction

## The value when the solver can act on it, and empty when it cannot.
function v = pick (used, value)
  if (used)
    v = value;
  else
    v = [];
  endif
endfunction

function v = ifelse (c, a, b)
  if (c)
    v = a;
  else
    v = b;
  endif
endfunction
