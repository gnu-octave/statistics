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
## @deftypefn {Private Function} {[@var{DP}, @var{K}, @var{S}, @var{W}] =} nbFit (@var{X}, @var{gY}, @var{k}, @var{D}, @var{w}, @var{Kernel}, @var{Support}, @var{Width}, @var{classname})
##
## Fit one univariate density per class and per predictor for a naive Bayes
## classifier.
##
## @var{X} holds the retained observations, @var{gY} indexes the class of each
## and @var{k} counts them, @var{D} names the distribution of each predictor
## and @var{w} weights each observation.  @var{CN} and @var{PN} name the
## classes and the predictors, and are used only to say which combination a
## refusal is about.  @var{Kernel}, @var{Support} and
## @var{Width} are the kernel arguments as the caller was given them, each
## possibly empty.
##
## @var{DP} is the class-by-predictor cell of fitted parameters: a column
## vector holding the mean and the standard deviation for a @qcode{'normal'}
## predictor, and a @code{prob.KernelDistribution} for a @qcode{'kernel'} one.
## @var{K}, @var{S} and @var{W} report the kernel, the support and the
## bandwidth actually used, empty where a predictor is not a kernel one.
##
## The standard deviation is the weighted, bias-corrected one, which reduces
## to the sample standard deviation when the weights within a class are equal,
## as they are whenever the prior is spread evenly over a class's members.
##
## @seealso{ClassificationNaiveBayes, fitcnb}
## @end deftypefn

function [DP, K, S, W, CL] = nbFit (X, gY, k, D, w, Kernel, Support, ...
                                    Width, classname, CN, PN)

  p = columns (X);
  DP = cell (k, p);
  K = cell (1, p);
  S = cell (1, p);
  CL = cell (1, p);
  W = [];

  ## The multinomial is one distribution over the whole predictor vector: a
  ## row is a vector of token counts, and a class is described by how it
  ## spends its tokens across the predictors.  Each class's counts are pooled
  ## and smoothed by one, so a token a class never spent is improbable rather
  ## than impossible, which a product over predictors would otherwise make it.
  if (ischar (D))
    for i = 1:k
      xk = X(gY == i, :);
      cnt = sum (xk, 1);
      tot = sum (cnt) + p;
      for j = 1:p
        DP{i,j} = (cnt(j) + 1) / tot;
      endfor
    endfor
    return;
  endif

  anykernel = iscell (D) && any (strcmp (D, 'kernel'));
  if (anykernel)
    [Kname, Sname, Wgiven] = nbKernelArgs (Kernel, Support, Width, k, p, ...
                                           classname);
    W = zeros (k, p);
  elseif (any (strcmp (D, 'mvmn')))
    ## A categorical predictor has no bandwidth, and MATLAB reports the
    ## absence as a matrix of NaN rather than as an empty one.
    W = NaN (k, p);
  endif

  for j = 1:p
    for i = 1:k

      idx = gY == i;
      xij = X(idx, j);
      wij = w(idx);

      if (isempty (xij))
        error (strcat (classname, ": class", sprintf (" %d", i), ...
                       " holds no observation."));
      endif

      switch (D{j})

        case 'mvmn'
          ## One categorical distribution per class and predictor, over the
          ## levels the predictor takes across the whole training set.  Every
          ## level's count is raised by one before it is normalized, so a
          ## level a class never took keeps a small probability instead of
          ## ruling that class out on this predictor alone.
          if (isempty (CL{j}))
            CL{j} = nbLevels (X(:,j));
          endif
          nlev = numel (CL{j});
          cnt = zeros (nlev, 1);
          for l = 1:nlev
            cnt(l) = sum (xij == CL{j}(l));
          endfor
          DP{i,j} = (cnt + 1) / (numel (xij) + nlev);

        case 'normal'
          sw = sum (wij);
          mu = sum (wij .* xij) / sw;
          ## The bias correction of a reliability weight: with equal weights
          ## the denominator is the count less one.
          den = sw - sum (wij .^ 2) / sw;
          sigma = 0;
          if (den > 0)
            sigma = sqrt (sum (wij .* (xij - mu) .^ 2) / den);
          endif
          ## A predictor that does not vary within a class has no normal
          ## density to fit: every observation sits on the mean, and the
          ## limit is a spike carrying no scale.  Refuse rather than answer,
          ## as MATLAB does, naming the combination it could not fit.  A
          ## kernel or a multivariate multinomial has no such difficulty and
          ## is not refused, which is a route left open to the caller.
          ## Test the data rather than the computed spread: the weighted
          ## mean of identical values need not land exactly on them, and a
          ## variance of 1e-15 where the true one is zero would slip a
          ## degenerate fit past a test of sigma alone.  The second half
          ## catches values that differ by so little that the squared
          ## deviation underflows.
          if (all (xij == xij(1)) || ! (sigma > 0))
            error (strcat (classname, ": a normal distribution cannot be", ...
                           " fit for the combination of class", ...
                           sprintf (" %s and predictor %s.", ...
                                    nbLabelName (CN, i), PN{j}), ...
                           " The data has zero variance."));
          endif
          DP{i,j} = [mu; sigma];

        case 'kernel'
          ## The density fits itself, so the default bandwidth rule and the
          ## transform a bounded support needs are ksdensity's own, not a
          ## second copy of them here.
          if (isempty (Wgiven))
            pd = prob.KernelDistribution.fit (xij, Kname{j}, Sname{j});
          else
            pd = prob.KernelDistribution.fit (xij, Kname{j}, Sname{j}, ...
                                              Wgiven(i,j));
          endif
          W(i,j) = pd.Bandwidth;
          DP{i,j} = pd;
          K{j} = Kname{j};
          S{j} = Sname{j};

      endswitch

    endfor
  endfor

endfunction
