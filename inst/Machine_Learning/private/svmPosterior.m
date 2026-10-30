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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{st}, @var{ST}, @var{separated}] =} svmPosterior (@var{s}, @var{k}, @var{prior}, @var{classnames})
##
## Fit the score-to-posterior transformation of a two-class SVM.
##
## @var{s} holds the scores of the second class, @code{@var{classnames}(2)},
## and @var{k} the class index, 1 or 2, of each observation; rows whose score
## is @qcode{NaN} or whose index is neither are left out.  @var{prior} is the
## model's prior.
##
## A zero in @var{prior} gives the constant transformation, posterior 1 for
## the other class, and @var{s} and @var{k} are not read.  Otherwise, scores
## of the two classes that do not overlap give the step function, 0 below the
## largest score of the first class, 1 above the smallest score of the second
## and @code{@var{prior}(2)} between, and @var{separated} is true.  Any other
## data give the sigmoid @code{1 / (1 + exp (A * s + B))}, fitted by
## maximum likelihood with @code{glmfit} to the class of each observation,
## unweighted.
##
## @var{st} is the transformation as the @qcode{ScoreTransform} property
## takes it, each parameter written in the fewest significant digits, 15 to
## 17, that read back to the same double, and @var{ST} the structure
## @code{fitPosterior} returns.
## @end deftypefn

function [st, ST, separated] = svmPosterior (s, k, prior, classnames)

  separated = false;

  ## One class cannot occur: every observation is the other
  if (any (prior == 0))
    c = find (prior > 0, 1);
    st = sprintf ('@(S)constant(S,%d)', c);
    ST = struct ('Type', 'constant', 'PredictedClass', {classnames(c)});
    return;
  endif

  keep = ! isnan (s(:)) & (k(:) == 1 | k(:) == 2);
  s = s(keep);
  y = (k(keep) == 2);

  ## Classes the scores separate: a step function
  if (any (y) && any (! y) && max (s(! y)) < min (s(y)))
    separated = true;
    lb = max (s(! y));
    ub = min (s(y));
    pp = prior(2);
    st = sprintf ('@(S)step(S,%s,%s,%s)', exact (lb), exact (ub), exact (pp));
    ST = struct ('Type', 'step', 'LowerBound', lb, 'UpperBound', ub, ...
                 'PositiveClassProbability', pp);
    return;
  endif

  ## Overlapping classes: logistic regression of the class on the score
  b = glmfit (s, double (y), 'binomial');
  A = -b(2);
  B = -b(1);
  st = sprintf ('@(S)sigmoid(S,%s,%s)', exact (A), exact (B));
  ST = struct ('Type', 'sigmoid', 'Slope', A, 'Intercept', B);

endfunction

## X in the fewest significant digits, from 15 to 17, that read back to X
## exactly, so the text is the parameter itself and stays short where it can.
function t = exact (x)
  for d = 15:16
    t = sprintf ('%.*g', d, x);
    if (str2double (t) == x)
      return;
    endif
  endfor
  t = sprintf ('%.17g', x);
endfunction
