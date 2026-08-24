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
##

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{L} =} nbLogLik (@var{X}, @var{D}, @var{DP}, @var{k})
##
## The class conditional log-likelihood of every observation under a fitted
## naive Bayes model.
##
## @var{L} is an observation-by-class matrix.  The predictors are conditionally
## independent given the class, so a row is the sum over the predictors of each
## one's log density, which is where the model's cost in the number of
## predictors stays linear.
##
## The sum is taken in logs rather than as a product of densities: with even a
## moderate number of predictors the product underflows to zero for every
## class at once, and the posterior that follows would be 0/0 rather than a
## well determined ratio.
##
## A row can come back as @math{-Inf} for every class, when a categorical
## predictor takes a level the model never saw.  That is a real answer, not a
## failure, and the caller reads it as such: the observation carries no
## information about the class and its posterior is the prior.
##
## A @qcode{NaN} predictor is skipped rather than propagated.  The predictors
## are conditionally independent given the class, so an observation missing
## one of them is still described by the others, and only its own term drops
## out.  An observation missing every predictor contributes nothing to any
## class and falls back to the prior, by the same arithmetic.  Measured
## against R2024a.
##
## @seealso{ClassificationNaiveBayes, fitcnb}
## @end deftypefn

function L = nbLogLik (X, D, DP, k, CL)

  n = rows (X);
  p = columns (X);
  L = zeros (n, k);

  ## The multinomial reads a row as token counts and scores it by how likely
  ## the class was to spend its tokens that way.  The multinomial coefficient
  ## depends on the row alone, not on the class, so it cancels in the
  ## posterior and is not computed.
  if (ischar (D))
    for i = 1:k
      for j = 1:p
        t = X(:,j) * log (DP{i,j});
        t(isnan (X(:,j))) = 0;
        L(:,i) += t;
      endfor
    endfor
    return;
  endif

  for i = 1:k
    for j = 1:p

      switch (D{j})

        case 'mvmn'
          ## An observation is scored by the probability of the level it
          ## takes.  A level the model never saw carries none, and makes the
          ## observation impossible under every class alike.
          lev = CL{j};
          pr = DP{i,j};
          lp = -Inf (n, 1);
          for l = 1:numel (lev)
            lp(X(:,j) == lev(l)) = log (pr(l));
          endfor
          lp(isnan (X(:,j))) = 0;
          L(:,i) += lp;

        case 'normal'
          mu = DP{i,j}(1);
          sigma = DP{i,j}(2);
          if (sigma > 0)
            z = (X(:,j) - mu) / sigma;
            t = -0.5 * z .^ 2 - log (sigma) - 0.5 * log (2 * pi);
            t(isnan (X(:,j))) = 0;
            L(:,i) += t;
          else
            ## A class holding one distinct value has no spread: the density
            ## is a spike, so an observation matching it is certain and every
            ## other one impossible.
            L(X(:,j) != mu, i) = -Inf;
          endif

        case 'kernel'
          t = log (pdf (DP{i,j}, X(:,j)));
          t(isnan (X(:,j))) = 0;
          L(:,i) += t;

      endswitch

    endfor
  endfor

endfunction
