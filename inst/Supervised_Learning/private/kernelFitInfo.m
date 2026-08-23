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
## @deftypefn {Private Function} {@var{S} =} kernelFitInfo (@var{info}, @var{LossFun}, @var{Lambda}, @var{BetaTol}, @var{GradTol})
##
## Assemble the @qcode{FitInfo} structure the kernel fitting functions
## return beside the model.
##
## It is not the structure the linear fitting functions return: a kernel
## model is always fitted by one solver under one regularization strength,
## so the fields are scalars and the ones that describe a choice of solver
## are absent.  The names differ too, @qcode{ObjectiveValue} and
## @qcode{GradientMagnitude} here against @qcode{Objective} and
## @qcode{GradientNorm} there, and both spellings are MATLAB's.
##
## @end deftypefn

function S = kernelFitInfo (info, LossFun, Lambda, BetaTol, GradTol)

  S = struct ();
  S.Solver = 'LBFGS-fast';
  S.LossFunction = LossFun;
  S.Lambda = Lambda;
  S.BetaTolerance = BetaTol;
  S.GradientTolerance = GradTol;
  S.ObjectiveValue = info.Objective;
  S.GradientMagnitude = info.GradientNorm;
  if (isempty (info.RelativeChangeInBeta))
    S.RelativeChangeInBeta = NaN;
  else
    S.RelativeChangeInBeta = info.RelativeChangeInBeta;
  endif
  S.FitTime = 0;
  S.History = [];

endfunction
