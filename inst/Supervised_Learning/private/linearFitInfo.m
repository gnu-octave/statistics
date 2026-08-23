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
## @deftypefn {Private Function} {@var{S} =} linearFitInfo (@var{info}, @var{BetaTol}, @var{GradTol}, @var{DeltaTol}, @var{IterLimit}, @var{Solver}, @var{PassLimit}, @var{BatchLimit})
##
## Assemble the @qcode{FitInfo} structure the linear fitting functions
## return beside the model.
##
## @var{info} holds one element per regularization strength, as
## @code{linearSolve} reported it, and every field of @var{S} is a row over
## those strengths, so a single strength gives scalars.
##
## The field set depends on the solver, as it does in MATLAB.  A
## quasi-Newton or proximal fit reports @qcode{IterationLimit}; a
## mini-batch or dual fit reports @qcode{PassLimit}, @qcode{NumPasses} and
## @qcode{BatchLimit} in its place, and adds @qcode{Alpha} for the dual and
## @qcode{BatchIndex} with @qcode{OptimalLearnRate} for the mini-batch
## solvers.  Neither of the latter two runs a gradient test, so both report
## a tolerance of zero and, by the rule that a tolerance of zero means the
## test did not run, a gradient of @qcode{NaN}.
##
## @end deftypefn

function S = linearFitInfo (info, BetaTol, GradTol, DeltaTol, IterLimit, ...
                            Solver, PassLimit, BatchLimit)

  stochastic = any (ismember (Solver, {'sgd', 'asgd', 'dual'}));

  S = struct ();
  S.Lambda = [info.Lambda];
  S.Objective = [info.Objective];
  if (stochastic)
    S.PassLimit = PassLimit;
    S.NumPasses = fieldOrNaN (info, 'NumPasses');
    S.BatchLimit = BatchLimit;
  else
    S.IterationLimit = IterLimit;
  endif
  S.NumIterations = [info.NumIterations];
  if (stochastic)
    S.GradientNorm = nan (1, numel (info));
    S.GradientTolerance = 0;
  else
    S.GradientNorm = [info.GradientNorm];
    S.GradientTolerance = GradTol;
  endif
  S.RelativeChangeInBeta = fieldOrNaN (info, 'RelativeChangeInBeta');
  S.BetaTolerance = BetaTol;
  if (any (strcmp (Solver, 'dual')))
    S.DeltaGradient = fieldOrNaN (info, 'DeltaGradient');
    S.DeltaGradientTolerance = DeltaTol;
  else
    S.DeltaGradient = [];
    S.DeltaGradientTolerance = [];
  endif
  S.TerminationCode = [info.TerminationCode];
  S.TerminationStatus = [info.TerminationStatus];
  if (any (strcmp (Solver, 'dual')))
    S.Alpha = [info.Alpha];
  endif
  if (any (ismember (Solver, {'sgd', 'asgd'})))
    S.BatchIndex = fieldOrNaN (info, 'BatchIndex');
    S.OptimalLearnRate = fieldOrNaN (info, 'OptimalLearnRate');
  endif
  S.History = [];
  S.FitTime = 0;
  S.Solver = Solver;

endfunction

## One row over the regularization strengths, with NaN wherever the solver
## that ran did not report the field at all.
function v = fieldOrNaN (info, name)

  v = nan (1, numel (info));
  for k = 1:numel (info)
    if (isfield (info(k), name) && ! isempty (info(k).(name)))
      v(k) = info(k).(name);
    endif
  endfor

endfunction
