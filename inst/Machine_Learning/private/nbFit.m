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
##

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{DP}, @var{K}, @var{S}, @var{W}, @var{CL}] =} nbFit (@var{X}, @var{gY}, @var{k}, @var{D}, @var{w}, @var{Kernel}, @var{Support}, @var{Width}, @var{classname}, @var{CN}, @var{PN})
##
## Fit one univariate density per class and per predictor for a naive Bayes
## classifier.
##
## @var{X} holds the retained observations, @var{gY} indexes the class of each
## and @var{k} counts them, @var{D} names the distribution of each predictor
## and @var{w} weights each observation.  @var{CN} and @var{PN} name the
## classes and the predictors, and are used only to say which combination a
## refusal is about.  @var{classname} prefixes any error raised.  @var{Kernel},
## @var{Support} and @var{Width} are the kernel arguments as the caller was
## given them, each possibly empty.
##
## @var{DP} is the class-by-predictor cell of fitted parameters: a column
## vector holding the mean and the standard deviation for a @qcode{'normal'}
## predictor, and a @code{prob.KernelDistribution} for a @qcode{'kernel'} one.
## @var{K}, @var{S} and @var{W} report the kernel, the support and the
## bandwidth actually used, empty where a predictor is not a kernel one.
## @var{CL} is a cell with the levels of each @qcode{'mvmn'} predictor in
## ascending order, empty for every other predictor.
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
      if (isempty (xk))
        continue;
      endif
      ## The multinomial is fitted on the rows whose every count is known, as
      ## R2026b fits it; a row missing a count still counts in the prior
      wk = w(gY == i);
      full = ! any (isnan (xk), 2);
      if (! any (full))
        error ("%s: class %s has no observation with every count known.", ...
               classname, nbLabelName (CN, i));
      endif
      xk = xk(full,:);
      wk = wk(full);
      ## Weighted counts, the weights scaled to the class's size so that the
      ## smoothing weighs as it does without them, as R2024a counts
      wk = wk(:) / sum (wk) * rows (xk);
      cnt = sum (wk .* xk, 1);
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
  else
    ## A predictor that is not kernel smoothed has no bandwidth, and MATLAB
    ## reports the absence as a matrix of NaN rather than as an empty one.
    ## It declares the matrix at its full size whatever the distributions.
    W = NaN (k, p);
  endif

  for j = 1:p
    for i = 1:k

      idx = gY == i;
      xij = X(idx, j);
      wij = w(idx);

      ## Only a class removed for its prior of 0 holds no observation, and
      ## its parameters stay empty, as MATLAB leaves them
      if (isempty (xij))
        continue;
      endif

      ## Each predictor is fitted on the observations it is known for, so a
      ## row missing one predictor still informs the others.  A class in
      ## which a predictor is never known has nothing to fit it on.
      obs = ! isnan (xij);
      if (! any (obs))
        error ("%s: predictor %s has no known value in class %s.", ...
               classname, PN{j}, nbLabelName (CN, i));
      endif
      xij = xij(obs);
      wij = wij(obs);

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
          ## Weighted counts, the weights scaled to the class's size so that
          ## the smoothing weighs as it does without them, as R2024a counts
          wk = wij(:) / sum (wij) * numel (xij);
          cnt = zeros (nlev, 1);
          for l = 1:nlev
            cnt(l) = sum (wk(xij == CL{j}(l)));
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
          ## Weights shape the density but not its bandwidth, which R2024a
          ## chooses from the observations alone.
          if (isempty (Wgiven))
            bw = prob.KernelDistribution.fit (xij, Kname{j}, ...
                                              Sname{j}).Bandwidth;
          else
            bw = Wgiven(i,j);
          endif
          freq = [];
          if (any (wij != wij(1)))
            freq = wij(:);
          endif
          pd = prob.KernelDistribution.fit (xij, Kname{j}, Sname{j}, bw, ...
                                            freq);
          W(i,j) = pd.Bandwidth;
          DP{i,j} = pd;
          K{j} = Kname{j};
          S{j} = Sname{j};

      endswitch

    endfor
  endfor

endfunction

## Resolve the kernel arguments of a naive Bayes classifier into one kernel and
## one support per predictor, and a bandwidth per class and predictor.
##
## Each of Kernel and Support may be given once for every predictor or once per
## predictor.  W is empty when no bandwidth was given, which leaves each
## density to choose its own.
function [K, S, W] = nbKernelArgs (Kernel, Support, Width, k, p, classname)

  known = {'box', 'epanechnikov', 'normal', 'triangle'};

  ## Kernel
  if (isempty (Kernel))
    K = repmat ({'normal'}, 1, p);
  elseif (ischar (Kernel) && isrow (Kernel))
    kname = lower (Kernel);
    K = repmat ({kname}, 1, p);
  elseif (iscellstr (Kernel) && numel (Kernel) == p)
    K = cellfun (@lower, Kernel(:)', 'UniformOutput', false);
  else
    error (strcat (classname, ": 'Kernel' must be a character vector or a", ...
                   " cell array of character vectors, one per predictor."));
  endif
  bad = ! ismember (K, known);
  if (any (bad))
    error (strcat (classname, ": unsupported kernel", ...
                   sprintf (" '%s'.", K{find (bad, 1)})));
  endif

  ## Support
  if (isempty (Support))
    S = repmat ({'unbounded'}, 1, p);
  elseif (ischar (Support) && isrow (Support))
    sname = lower (Support);
    S = repmat ({sname}, 1, p);
  elseif (isnumeric (Support) && numel (Support) == 2)
    sbnds = sort (Support(:)');
    S = repmat ({sbnds}, 1, p);
  elseif (iscell (Support) && numel (Support) == p)
    S = Support(:)';
  else
    error (strcat (classname, ": 'Support' must be 'unbounded',", ...
                   " 'positive', a two element numeric vector, or a cell", ...
                   " array holding one of those per predictor."));
  endif
  for j = 1:p
    if (ischar (S{j}) && ! any (strcmpi (S{j}, {'positive', 'unbounded'})))
      error (strcat (classname, ": a character vector 'Support' must be", ...
                     " 'unbounded' or 'positive'."));
    endif
  endfor

  ## Width
  if (isempty (Width))
    W = [];
    return;
  endif
  if (! (isnumeric (Width) && isreal (Width) && all (Width(:) > 0)))
    error (strcat (classname, ": 'Width' must be positive and real."));
  endif
  if (isscalar (Width))
    W = repmat (Width, k, p);
  elseif (isequal (size (Width), [k, p]))
    W = Width;
  elseif (isrow (Width) && numel (Width) == p)
    W = repmat (Width(:)', k, 1);
  elseif (iscolumn (Width) && numel (Width) == k)
    W = repmat (Width(:), 1, p);
  else
    error (strcat (classname, ": 'Width' must be a scalar, one value per", ...
                   " predictor, one per class, or a matrix of one per", ...
                   " class and predictor."));
  endif

endfunction

## One class's name, as a message should print it.
##
## C is the model's ClassNames in any of the forms a response may take, and i
## indexes it.  A textual class prints as its name and every other kind as its
## value, which is what MATLAB's own messages do: a numeric response names
## class 2 where a cellstr one names class beta.
function s = nbLabelName (C, i)

  if (iscellstr (C))
    s = C{i};
  elseif (ischar (C))
    s = strtrim (C(i,:));
  elseif (islogical (C))
    s = mat2str (C(i));
  else
    s = num2str (C(i));
  endif

endfunction

## The distinct levels of one categorical predictor, in ascending order.  The
## order is the one every per-level quantity is reported in, so it is fixed
## here once rather than recovered wherever a level is looked up.
function L = nbLevels (x)

  L = unique (x(:));
  L = L(! isnan (L));

endfunction
