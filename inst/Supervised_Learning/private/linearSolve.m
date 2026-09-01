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
## @deftypefn {Private Function} {[@var{Beta}, @var{Bias}, @var{S}] =} linearSolve (@var{X}, @var{y}, @var{w}, @var{P})
##
## Minimize the regularized linear objective shared by the four linear and
## kernel learners.
##
## The objective is
##
## @example
## sum (@var{w} .* L (@var{X} * Beta + Bias, @var{y})) + Lambda * R (Beta)
## @end example
##
## with @var{w} summing to one, @code{L} the fitted loss of
## @code{linearLoss}, and @code{R} either @math{||Beta||^2 / 2} for ridge or
## @math{||Beta||_1} for lasso.  The bias is never regularized.
##
## @var{P} is a scalar struct holding the resolved fitting options; the
## caller has already validated them and turned every @qcode{'auto'} into a
## number.  @var{S} reports the fit: @qcode{Objective},
## @qcode{NumIterations}, @qcode{GradientNorm}, @qcode{RelativeChangeInBeta},
## @qcode{DeltaGradient}, @qcode{TerminationCode}, @qcode{TerminationStatus},
## @qcode{History} and @qcode{Solver}.
##
## @end deftypefn

function [Beta, Bias, S] = linearSolve (X, y, w, P)

  Beta = P.InitialBeta(:);
  Bias = P.InitialBias;
  if (! P.FitBias)
    Bias = 0;
  endif

  switch (P.Solver)

    case {'lbfgs', 'bfgs'}
      [Beta, Bias, S] = solveQuasiNewton (X, y, w, P, Beta, Bias);

    case 'sparsa'
      [Beta, Bias, S] = solveSpaRSA (X, y, w, P, Beta, Bias);

    case {'sgd', 'asgd'}
      [Beta, Bias, S] = solveSGD (X, y, w, P, Beta, Bias);

    case 'dual'
      [Beta, Bias, S] = solveDual (X, y, w, P, Beta, Bias);

    otherwise
      error ("linearSolve: unknown solver '%s'.", P.Solver);

  endswitch

  ## An explicit refit of the bias against the fitted scores, once the
  ## coefficients are settled.  MATLAB calls this PostFitBias and leaves it
  ## off by default.
  if (P.FitBias && P.PostFitBias)
    Bias = refitBias (X * Beta, y, w, P);
  endif

  S.Objective = linearObjective (X, y, w, P, Beta, Bias);
  S.Solver = P.Solver;

  ## A tolerance of exactly zero does not mean "test against zero", it means
  ## the test does not run, and MATLAB then reports the quantity it governs
  ## as NaN rather than as a number.  Measured on R2024a across a twelve
  ## point sweep of fitclinear.  Note that the engine's LossTolerance is
  ## switched off with -Inf instead: the two conventions differ because
  ## MATLAB's do, each mirroring the quantity it governs, and unifying them
  ## would create a divergence rather than close one.
  if (P.GradientTolerance == 0)
    S.GradientNorm = NaN;
  endif
  if (P.BetaTolerance == 0)
    S.RelativeChangeInBeta = NaN;
  endif

endfunction

## Objective value alone, for reporting.
function f = linearObjective (X, y, w, P, Beta, Bias)

  L = linearLoss (X * Beta + Bias, y, P.LossFunction, P.Epsilon);
  f = sum (w .* L) + P.Lambda * penalty (Beta, P.Regularization);

endfunction

function r = penalty (Beta, reg)

  if (strcmp (reg, 'ridge'))
    r = sum (Beta .^ 2) / 2;
  else
    r = sum (abs (Beta));
  endif

endfunction

## Value and gradient of the smooth part, over the packed parameter vector.
## The lasso penalty is not part of it: SpaRSA handles that one by its prox.
function [f, g] = smoothPart (z, X, y, w, P)

  p = columns (X);
  Beta = z(1:p);
  if (P.FitBias)
    Bias = z(p+1);
  else
    Bias = 0;
  endif

  [L, dL] = linearLoss (X * Beta + Bias, y, P.LossFunction, P.Epsilon);
  f = sum (w .* L);
  v = w .* dL;
  gBeta = X' * v;

  if (strcmp (P.Regularization, 'ridge'))
    f += P.Lambda * sum (Beta .^ 2) / 2;
    gBeta += P.Lambda * Beta;
  endif

  if (P.FitBias)
    gBias = sum (v);
    g = [gBeta; gBias];
  else
    g = gBeta;
  endif

endfunction

## LBFGS and BFGS, both served by the package's __lbfgs__ engine.
function [Beta, Bias, S] = solveQuasiNewton (X, y, w, P, Beta, Bias)

  p = columns (X);
  z0 = Beta;
  if (P.FitBias)
    z0 = [Beta; Bias];
  endif

  opt = struct ();
  opt.IterationLimit = P.IterationLimit;
  opt.GradientTolerance = P.GradientTolerance;
  opt.StepTolerance = 0;
  ## LossTolerance MUST stay -Inf.  The engine's test is on the objective
  ## VALUE rather than on its change, which is MATLAB's convention and not
  ## what the name suggests, so any fit whose objective passes below 1e-6 on
  ## its way down stops there with the coefficients still wrong.  Measured:
  ## an exact linear relation fitted at Lambda 1e-10 stops after 5
  ## iterations with the coefficients out by 1.6e-4, against 11 iterations
  ## and 2.6e-10 with the test switched off.  The BIST named
  ## 'LossTolerance' in RegressionLinear.m fails if this line is removed.
  ##
  ## Note that BetaTolerance is switched off with 0 and this one with -Inf.
  ## The inconsistency is deliberate: each mirrors MATLAB's own convention
  ## for the quantity it governs, and unifying them would create a
  ## divergence rather than close one.
  opt.LossTolerance = -Inf;
  opt.HistorySize = P.HessianHistorySize;
  opt.BetaTolerance = P.BetaTolerance;

  fcn = @(z) smoothPart (z, X, y, w, P);
  [z, info] = __lbfgs__ (fcn, z0, opt);

  Beta = z(1:p);
  if (P.FitBias)
    Bias = z(p+1);
  endif

  S = struct ();
  S.NumIterations = info.Iterations;
  S.GradientNorm = info.Gradient;
  S.RelativeChangeInBeta = info.RelativeChangeInBeta;
  S.DeltaGradient = [];
  [S.TerminationCode, S.TerminationStatus] = terminationOf (info);
  S.History = [];

endfunction

## SpaRSA, the separable-approximation proximal method MATLAB uses for the
## lasso penalty.  Each step minimizes the loss's linear model plus an
## isotropic quadratic of curvature alpha, which is the soft threshold, with
## alpha set by the Barzilai-Borwein ratio and grown until a nonmonotone
## sufficient-decrease test passes.
function [Beta, Bias, S] = solveSpaRSA (X, y, w, P, Beta, Bias)

  p = columns (X);
  memory = 5;
  eta = 2;
  sigma = 0.01;
  alphaMin = 1e-30;
  alphaMax = 1e30;

  z = Beta;
  if (P.FitBias)
    z = [Beta; Bias];
  endif
  [f, g] = smoothPart (z, X, y, w, P);
  obj = f + P.Lambda * penalty (z(1:p), P.Regularization);
  past = obj;
  alpha = 1;

  code = 0;
  relChange = Inf;
  iter = 0;
  for iter = 1:P.IterationLimit

    accepted = false;
    for ls = 1:60
      zNew = z - g / alpha;
      zNew(1:p) = softThreshold (zNew(1:p), P.Lambda / alpha);
      d = zNew - z;
      [fNew, gNew] = smoothPart (zNew, X, y, w, P);
      objNew = fNew + P.Lambda * penalty (zNew(1:p), P.Regularization);
      if (objNew <= max (past) - sigma * alpha / 2 * sum (d .^ 2))
        accepted = true;
        break;
      endif
      alpha *= eta;
    endfor
    if (! accepted)
      code = 0;
      break;
    endif

    ## Relative change in the coefficients, MATLAB's primary stopping test.
    relChange = relativeChange (zNew, z);

    ## Barzilai-Borwein curvature for the next step.
    dg = gNew - g;
    dd = sum (d .^ 2);
    if (dd > 0)
      alpha = (d' * dg) / dd;
      if (! isfinite (alpha) || alpha < alphaMin)
        alpha = alphaMin;
      elseif (alpha > alphaMax)
        alpha = alphaMax;
      endif
    endif

    z = zNew;
    f = fNew;
    g = gNew;
    obj = objNew;
    past = [past, obj];
    if (numel (past) > memory)
      past(1) = [];
    endif

    if (relChange <= P.BetaTolerance)
      code = 1;
      break;
    endif

  endfor

  Beta = z(1:p);
  if (P.FitBias)
    Bias = z(p+1);
  endif

  ## The gradient reported for a lasso fit is that of the smooth part alone,
  ## the penalty having no gradient at the zeros the fit is there to produce.
  S = struct ();
  S.NumIterations = iter;
  S.GradientNorm = max (abs (g));
  S.RelativeChangeInBeta = relChange;
  S.DeltaGradient = [];
  [S.TerminationCode, S.TerminationStatus] = terminationOfCode (code);
  S.History = [];

endfunction

## Stochastic gradient descent, and its averaged variant.  The learning rate
## decays as LearnRate / (1 + Lambda * LearnRate * t) for a ridge penalty and
## stays put for a lasso one, which is soft-thresholded every
## TruncationPeriod mini-batches rather than every one of them.
function [Beta, Bias, S] = solveSGD (X, y, w, P, Beta, Bias)

  [n, p] = size (X);
  averaged = strcmp (P.Solver, 'asgd');
  batch = min (P.BatchSize, n);
  perPass = ceil (n / batch);
  limit = P.PassLimit * perPass;
  if (! isempty (P.BatchLimit))
    limit = min (limit, P.BatchLimit);
  endif

  rate = P.LearnRate;
  ridge = strcmp (P.Regularization, 'ridge');
  betaSum = zeros (p, 1);
  biasSum = 0;
  count = 0;
  t = 0;
  relChange = Inf;
  prev = Beta;
  if (P.FitBias)
    prev = [Beta; Bias];
  endif
  lastObj = Inf;

  for pass = 1:P.PassLimit
    order = randperm (n);
    for b = 1:perPass
      t++;
      if (t > limit)
        break;
      endif
      idx = order((b - 1) * batch + 1 : min (b * batch, n));
      Xb = X(idx,:);
      yb = y(idx);
      wb = w(idx);
      wb = wb / sum (wb);

      [~, dL] = linearLoss (Xb * Beta + Bias, yb, P.LossFunction, P.Epsilon);
      v = wb .* dL;
      gBeta = Xb' * v;
      if (ridge)
        gBeta += P.Lambda * Beta;
        step = rate / (1 + P.Lambda * rate * t);
      else
        step = rate;
      endif

      Beta -= step * gBeta;
      if (P.FitBias)
        Bias -= step * sum (v);
      endif
      if (! ridge && mod (t, P.TruncationPeriod) == 0)
        Beta = softThreshold (Beta, step * P.TruncationPeriod * P.Lambda);
      endif

      if (averaged)
        betaSum += Beta;
        biasSum += Bias;
        count++;
      endif
    endfor

    ## What the pass reached, which drives both the stopping test and, when
    ## OptimizeLearnRate is on, the step size of the next pass: MATLAB
    ## halves the rate whenever a pass leaves the objective higher than it
    ## found it.
    if (P.FitBias)
      relChange = relativeChange ([Beta; Bias], prev);
      prev = [Beta; Bias];
      obj = linearObjective (X, y, w, P, Beta, Bias);
    else
      relChange = relativeChange (Beta, prev);
      prev = Beta;
      obj = linearObjective (X, y, w, P, Beta, 0);
    endif
    if (P.OptimizeLearnRate && mod (pass, P.NumCheckConvergence) == 0)
      if (obj > lastObj)
        rate /= 2;
      endif
      lastObj = obj;
    endif

    if (t > limit || relChange <= P.BetaTolerance)
      break;
    endif
  endfor

  if (averaged && count > 0)
    Beta = betaSum / count;
    Bias = biasSum / count;
  endif

  z = Beta;
  if (P.FitBias)
    z = [Beta; Bias];
  endif
  [~, g] = smoothPart (z, X, y, w, P);

  S = struct ();
  S.NumIterations = t;
  S.NumPasses = pass;
  S.BatchIndex = t;
  S.OptimalLearnRate = rate;
  S.GradientNorm = max (abs (g));
  S.RelativeChangeInBeta = relChange;
  S.DeltaGradient = [];
  [S.TerminationCode, S.TerminationStatus] = terminationOfCode (0);
  S.History = [];

endfunction

## Dual coordinate descent, for a ridge-penalized hinge or
## epsilon-insensitive loss.  One coordinate of the dual is minimized
## exactly at a time and the primal coefficients are carried alongside, so a
## pass costs one sweep over the observations.  Convergence is judged by the
## largest violation of the coordinate optimality conditions, which is what
## MATLAB reports as DeltaGradient.
##
## The intercept enters as a constant predictor, so this solver penalizes it
## where the quasi-Newton solvers leave it free.  That is what a coordinate
## dual can do without the equality constraint an unpenalized intercept
## imposes, it is what LIBLINEAR does, and MATLAB's own dual solver differs
## from its own quasi-Newton one in the same direction and by more: on the
## two overlapping iris species R2024a reaches 0.2760 by 'dual' against
## 0.1576 by 'bfgs'.  Expect the two solvers to disagree, and prefer the
## quasi-Newton answer when the objective as written is what matters.
function [Beta, Bias, S] = solveDual (X, y, w, P, Beta, Bias)

  [n, p] = size (X);
  classification = any (strcmp (P.LossFunction, {'hinge', 'logit'}));
  U = w / P.Lambda;
  sq = sum (X .^ 2, 2);
  ## The intercept enters as a constant predictor, and the value that
  ## constant takes decides how heavily the dual penalizes it: an intercept
  ## fitted through a feature of value B carries a penalty of (Bias / B)^2 / 2
  ## rather than Bias^2 / 2.  LIBLINEAR leaves that value at one and so
  ## regularizes the intercept as hard as a coefficient, which is visibly
  ## wrong when the response is far from zero: on carsmall, whose intercept
  ## is around 24, it costs more than the fit it buys.  Setting it to the
  ## root mean square length of an observation makes the penalty negligible
  ## and the intercept effectively free, which is the objective as written.
  ## It is not free of cost: a constant that large enters every coordinate's
  ## curvature and slows the sweep down, so a fit whose intercept is far from
  ## zero wants a generous PassLimit.  On carsmall the objective falls from
  ## 5.108 at 500 passes to 2.415 at 8000, the last of those below the 2.827
  ## the quasi-Newton solvers stall at.
  biasScale = 1;
  if (P.FitBias)
    biasScale = max (1, sqrt (mean (sum (X .^ 2, 2))));
    sq += biasScale ^ 2;
  endif
  sq = max (sq, eps);

  alpha = zeros (n, 1);
  Beta = zeros (p, 1);
  Bias = 0;
  delta = Inf;
  pass = 0;

  for pass = 1:P.PassLimit
    order = randperm (n);
    delta = 0;
    for k = order
      xk = X(k,:);
      f = xk * Beta + Bias;
      aOld = alpha(k);
      if (classification)
        ## The exact minimizer of the dual along this coordinate, clipped to
        ## its own box.
        g = y(k) * f - 1;
        aNew = min (max (aOld - g / sq(k), 0), U(k));
        step = (aNew - aOld) * y(k);
      else
        ## The same, with the insensitive band entering as an absolute value
        ## that the step is soft thresholded by.
        g = f - y(k);
        u = aOld - g / sq(k);
        t = P.Epsilon / sq(k);
        aNew = sign (u) * max (abs (u) - t, 0);
        aNew = min (max (aNew, -U(k)), U(k));
        step = aNew - aOld;
      endif

      ## How far the coordinate was from its own optimum, in the units of
      ## the gradient.  A coordinate pinned against a bound by a gradient
      ## pushing it outwards does not move and is not in violation.
      delta = max (delta, abs (aNew - aOld) * sq(k));

      if (step != 0)
        Beta += step * xk';
        if (P.FitBias)
          Bias += step * biasScale;
        endif
        alpha(k) = aNew;
      endif
    endfor
    if (mod (pass, P.NumCheckConvergence) == 0
        && delta <= P.DeltaGradientTolerance)
      break;
    endif
  endfor

  S = struct ();
  S.NumIterations = pass * n;
  S.NumPasses = pass;
  S.Alpha = alpha;
  z = Beta;
  if (P.FitBias)
    z = [Beta; Bias];
  endif
  [~, g] = smoothPart (z, X, y, w, P);
  S.GradientNorm = max (abs (g));
  S.RelativeChangeInBeta = [];
  S.DeltaGradient = delta;
  if (delta <= P.DeltaGradientTolerance)
    [S.TerminationCode, S.TerminationStatus] = terminationOfCode (4);
  else
    [S.TerminationCode, S.TerminationStatus] = terminationOfCode (0);
  endif
  S.History = [];

endfunction

function b = softThreshold (b, t)
  b = sign (b) .* max (abs (b) - t, 0);
endfunction

## The relative change MATLAB tests BetaTolerance against.  Measured on
## R2024a: the two-norm of the step over the two-norm of the iterate, the
## bias inside both vectors and no guard on the denominator.  An iterate that
## is identically zero has not moved, so the quotient is 0/0 and NaN, which
## is the value MATLAB reports for this quantity whenever it is not tested.
function r = relativeChange (new, old)
  d = new - old;
  r = norm (d) / norm (new);
endfunction

## Refit the bias alone against the fitted linear part.  Least squares and
## the logit have an interior optimum found by a scalar Newton step; the two
## non-smooth losses are minimized over the finite set of points at which
## their derivative can change sign.
function Bias = refitBias (f, y, w, P)

  switch (P.LossFunction)

    case 'mse'
      Bias = sum (w .* (y - f)) / sum (w);

    case 'epsiloninsensitive'
      knots = unique ([y - f - P.Epsilon; y - f + P.Epsilon]);
      Bias = bestOf (knots, f, y, w, P);

    case 'hinge'
      knots = unique ((1 ./ y) - f);
      Bias = bestOf (knots, f, y, w, P);

    case 'logit'
      Bias = 0;
      for k = 1:50
        [~, dL] = linearLoss (f + Bias, y, P.LossFunction, P.Epsilon);
        g = sum (w .* dL);
        s = 1 ./ (1 + exp (-(y .* (f + Bias))));
        h = sum (w .* s .* (1 - s));
        if (h <= 0)
          break;
        endif
        step = g / h;
        Bias -= step;
        if (abs (step) <= 1e-12)
          break;
        endif
      endfor

  endswitch

endfunction

function b = bestOf (knots, f, y, w, P)

  best = Inf;
  b = 0;
  for k = 1:numel (knots)
    L = linearLoss (f + knots(k), y, P.LossFunction, P.Epsilon);
    v = sum (w .* L);
    if (v < best)
      best = v;
      b = knots(k);
    endif
  endfor

endfunction

## Map what the engine reports onto MATLAB's TerminationCode.  The engine
## hands out a stable token rather than a sentence, because the wording is
## the caller's: fitcnet says 'Relative gradient tolerance reached.' where
## fitclinear says 'Tolerance on gradient satisfied.' for the same event.  A
## build that offers no token has nothing to map, and reports the iteration
## limit.
##
## Measured on R2024a fitclinear: the coefficient tolerance is code 1, the
## gradient tolerance 2, the iteration limit 0, and a line search that cannot
## improve the objective -11, worded 'Unable to find a step decreasing the
## objective.'  'loss' and 'step' are unreachable from here: LossTolerance is
## -Inf and StepTolerance 0, and a failed line search returns before the step
## test is ever reached, so a step of exactly zero cannot arise.
function [code, status] = terminationOf (info)

  token = '';
  if (isfield (info, 'Criterion'))
    token = info.Criterion;
  endif

  switch (token)
    case 'beta'
      code = 1;
    case {'gradient', 'step'}
      code = 2;
    case 'linesearch'
      code = -11;
    otherwise
      code = 0;
  endswitch
  [code, status] = terminationOfCode (code);

endfunction

function [code, status] = terminationOfCode (code)

  switch (code)
    case 1
      status = {'Tolerance on coefficients satisfied.'};
    case 2
      status = {'Tolerance on gradient satisfied.'};
    case 4
      status = {'Tolerance on the complementarity gap satisfied.'};
    case -11
      status = {'Unable to find a step decreasing the objective.'};
    otherwise
      code = 0;
      status = {'Iteration limit exceeded.'};
  endswitch

endfunction
