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
## @deftypefn  {Private Function} gpPredict (@var{XC}, @var{M})
## @deftypefnx {Private Function} [@var{ypred}, @var{ysd}, @var{yint}] = gpPredict (@var{XC}, @var{M})
##
## Evaluate a fitted Gaussian process at new points.
##
## @var{M} is a structure carrying the fitted model: the active set vectors it
## predicts from, the prediction weights @qcode{Alpha}, the covariance function
## and its parameters, the explicit basis and its coefficients, the noise
## standard deviation, the standardizing location and scale where the model
## standardized its predictors, and the significance level of the intervals.
##
## The mean is the basis term plus the covariance between the new points and
## the active set weighted by @qcode{Alpha}, which is what @qcode{Alpha} is
## for.  The standard deviation is that of the @emph{response}, so it carries
## the noise as well as the uncertainty of the latent function, and the
## interval is the normal quantile of the level times it.  The coefficients of
## the explicit basis are taken as known, which is measurably what MATLAB does:
## carrying their covariance as well moves the answer in the sixth decimal.
##
## @end deftypefn

function [ypred, ysd, yint] = gpPredict (XC, M)

  ## Standardize the new points exactly as the training data was standardized
  if (! isempty (M.Location))
    XC = (XC - M.Location) ./ M.Scale;
  endif

  ## Mean prediction
  Kq = gpCov (XC, M.X, M);
  HC = gpBasis (XC, M.BasisFunction);
  ypred = Kq * M.Alpha;
  if (! isempty (M.Beta))
    ypred += HC * M.Beta(:);
  endif

  if (nargout < 2)
    return;
  endif

  ## The variance of a new response is the prior variance of the latent
  ## function, less what the training data explains of it, plus the noise.
  ## The factorisation is rebuilt here rather than stored: a compact model
  ## keeps its active set and nothing else, and this is the only method that
  ## needs it.
  n = rows (M.X);
  A = gpCov (M.X, M.X, M) + M.Sigma^2 * eye (n);
  V = A \ Kq';
  kqq = gpSelfCov (XC, M);
  vy = kqq - sum (Kq' .* V, 1)' + M.Sigma^2;
  vy(vy < 0) = 0;
  ysd = sqrt (vy);

  if (nargout < 3)
    return;
  endif

  z = norminv (1 - M.CIAlpha / 2);
  yint = [ypred - z * ysd, ypred + z * ysd];

endfunction

## Covariance between two sets of points, built in or user supplied.
function K = gpCov (A, B, M)
  if (is_function_handle (M.KernelFunction))
    K = M.KernelFunction (A, B, M.Theta);
  else
    K = gpKernel (A, B, M.KernelFunction, M.Theta);
  endif
endfunction

## The prior variance of a point with itself.  Every built-in kernel is
## stationary, so this is the signal variance and does not depend on where the
## point is; a supplied kernel need not be, and is asked.
function v = gpSelfCov (XC, M)
  if (is_function_handle (M.KernelFunction))
    v = zeros (rows (XC), 1);
    for i = 1:rows (XC)
      v(i) = M.KernelFunction (XC(i,:), XC(i,:), M.Theta);
    endfor
  else
    v = repmat (M.Theta(end)^2, rows (XC), 1);
  endif
endfunction
