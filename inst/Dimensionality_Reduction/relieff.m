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
## @deftypefn  {statistics} {@var{idx} =} relieff (@var{X}, @var{Y}, @var{k})
## @deftypefnx {statistics} {@var{idx} =} relieff (@var{X}, @var{Y}, @var{k}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{idx}, @var{weights}] =} relieff (@dots{})
##
## Rank the importance of predictors with the ReliefF or RReliefF algorithm.
##
## @code{@var{idx} = relieff (@var{X}, @var{Y}, @var{k})} ranks the
## predictors, the columns of the matrix @var{X}, by how well they tell apart
## observations with different responses @var{Y}, judged by the @var{k}
## nearest neighbors of each observation.  @var{idx} is a row vector of
## column indices, from the most important predictor to the least.
##
## @var{X} is a matrix of type double or single, with one row per observation.
## @var{Y} holds one response per row of @var{X}: a numeric vector, a logical
## vector, a categorical vector, a cell array of character vectors, a string
## array, or a character matrix with one row per observation.  A numeric
## @var{Y} is a regression problem and any other type a classification
## problem, unless @qcode{'method'} says otherwise.  Rows with a missing value
## (@qcode{NaN}) in @var{X} or a missing response (@qcode{NaN}, an empty or
## blank name, or an undefined category) are left out.  @var{k} is a positive
## integer.
##
## @code{[@var{idx}, @var{weights}] = relieff (@dots{})} also returns the
## weight of each predictor, a row vector in the order of the columns of
## @var{X}, of type double.  A larger weight marks a more important predictor,
## and weights range from -1 to 1.  @var{idx} sorts @var{weights} in
## descending order.
##
## Predictors are compared on their range: for continuous predictors, the
## difference between two observations in predictor @math{j} is
## @math{|x_rj - x_qj|} divided by the range of predictor @math{j}, and the
## distance between them is the sum of these differences over all predictors
## (city block distance).  A predictor holding a single value is left out of
## the distance and has weight @qcode{NaN}, ranked last.  For categorical
## predictors the difference is 0 when two observations share a value and 1
## otherwise.  The neighbors of an observation are the @var{k} observations
## nearest to it, ties between equal distances going to the observation that
## comes first in @var{X}.  Each neighbor counts with a weight that falls with
## its rank @math{i} among them as @math{exp (-(i / sigma)^2)}, normalized to
## sum to 1 over the neighbors.
##
## For classification (ReliefF), each observation @math{r} has @var{k}
## nearest neighbors in its own class (hits) and @var{k} in each of the other
## classes (misses), or the whole class when it holds fewer.  The weight of
## predictor @math{j} falls by the weighted difference to the hits and rises
## by the weighted difference to the misses of class @math{c}, scaled by
## @math{p_c / (1 - p_r)}, @math{p} being the prior probabilities of the
## classes; each observation adds @math{1/m} of its own share, @math{m} being
## the number of observations used.  For regression (RReliefF), with
## @math{dy} the difference of two responses divided by the range of @var{Y}
## and @math{dj} the difference in predictor @math{j}, the weighted sums of
## @math{dy}, @math{dj} and @math{dy dj} over the @var{k} nearest neighbors
## of each observation are accumulated as @math{W_dy}, @math{W_dj} and
## @math{W_dydj}, and the weight of predictor @math{j} is
## @math{W_dydj / W_dy - (W_dj - W_dydj) / (m - W_dy)}.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.2 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'method'} @tab @qcode{'classification'} or
## @qcode{'regression'}.  The default is @qcode{'regression'} for a numeric
## @var{Y} and @qcode{'classification'} otherwise.  A numeric @var{Y} may be
## classified, its distinct values being the classes, but only a numeric
## @var{Y} can be regressed.
##
## @item @qcode{'prior'} @tab The prior probabilities of the classes, for
## classification only: @qcode{'empirical'}, the class frequencies in @var{Y},
## which is the default; @qcode{'uniform'}, equal for every class; a numeric
## vector with one value per class; or a structure with the fields
## @qcode{group}, the class names, and @qcode{prob}, their probabilities.  A
## numeric vector follows the order of the classes in @var{Y}: sorted values
## for a numeric, logical or categorical @var{Y}, and order of first
## appearance for text.  A structure must name every class of @var{Y} and may
## name others, which are ignored.  The probabilities are normalized to sum
## to 1, and the observations of a class with prior probability 0 are left
## out.
##
## @item @qcode{'updates'} @tab The number of observations, drawn at random
## without replacement, for which neighbors are sought: a positive integer, or
## @qcode{'all'}, the default, for every observation.  A number at least as
## large as the number of observations uses them all.
##
## @item @qcode{'categoricalx'} @tab @qcode{'on'} to treat every predictor as
## categorical, or @qcode{'off'}, the default, to treat every predictor as
## continuous.
##
## @item @qcode{'sigma'} @tab The distance scaling factor of the rank weights
## of the neighbors, a positive scalar.  The default is @code{Inf} for
## classification, weighing all neighbors alike, and 50 for regression.
## @end multitable
##
## When @var{X} is single, the weights are computed in single precision and
## returned as double.  The option names are not case sensitive but must be
## given in full, where MATLAB accepts any shortened name.  Where MATLAB
## accepts a non-integer @var{k} or @qcode{'updates'}, answers for complex or
## infinite values in @var{X}, gives weights for a single class of responses,
## and returns @qcode{NaN} for a constant response, @code{relieff} refuses
## them.  It also accepts numeric class names in a @qcode{'prior'} structure,
## which MATLAB refuses.
##
## References:
## @enumerate
## @item
## I. Kononenko, E. @v{S}imec and M. Robnik-@v{S}ikonja.  Overcoming the
## myopia of inductive learning algorithms with RELIEFF.
## @emph{Applied Intelligence}, 7(1):39--55, 1997.
## @item
## M. Robnik-@v{S}ikonja and I. Kononenko.  Theoretical and empirical
## analysis of ReliefF and RReliefF.  @emph{Machine Learning},
## 53(1-2):23--69, 2003.
## @end enumerate
##
## @seealso{pca, fitcknn, knnsearch}
## @end deftypefn

function [idx, weights] = relieff (X, Y, k, varargin)

  if (nargin < 3)
    print_usage ();
  endif
  if (! (isfloat (X) && isreal (X) && ismatrix (X) && ! isempty (X)))
    error ("relieff: X must be a real matrix of type double or single.");
  endif
  if (ischar (Y))
    nY = rows (Y);
  elseif ((isnumeric (Y) || islogical (Y) || iscategorical (Y)
           || iscellstr (Y) || isstring (Y)) && isvector (Y))
    nY = numel (Y);
  else
    error ("relieff: Y must be a vector of responses or class names.");
  endif
  if (rows (X) != nY)
    error ("relieff: X must have as many rows as Y has observations.");
  endif
  if (! (isnumeric (k) && isreal (k) && isscalar (k) && isfinite (k)
         && k > 0 && k == fix (k)))
    error ("relieff: K must be a positive integer scalar.");
  endif

  ## Parse optional Name-Value paired arguments
  optNames = {'method', 'prior', 'updates', 'categoricalx', 'sigma'};
  dfValues = {[], [], 'all', 'off', []};
  [method, prior, updates, catX, sigma, args] = ...
                    parsePairedArguments (optNames, dfValues, varargin(:));
  if (! isempty (args))
    if (mod (numel (args), 2) == 0 && ischar (args{1}))
      error ("relieff: unknown parameter name '%s'.", args{1});
    endif
    error ("relieff: optional arguments must be Name-Value pairs.");
  endif
  if (isempty (method))
    if (isnumeric (Y))
      method = 'regression';
    else
      method = 'classification';
    endif
  elseif (! ((ischar (method) && isrow (method))
             || (isstring (method) && isscalar (method))))
    error ("relieff: 'method' must be a character vector.");
  endif
  method = lower (char (method));
  if (! any (strcmp (method, {'classification', 'regression'})))
    error ("relieff: unknown method '%s'.", method);
  endif
  isreg = strcmp (method, 'regression');
  if (isreg && ! isnumeric (Y))
    error ("relieff: regression needs a numeric Y.");
  endif
  if (isreg && ! isempty (prior))
    error ("relieff: 'prior' applies to classification only.");
  endif
  if (isempty (prior))
    prior = 'empirical';
  endif
  if (isstring (updates) && isscalar (updates))
    updates = char (updates);
  endif
  if (ischar (updates))
    if (! strcmpi (updates, 'all'))
      error ("relieff: 'updates' must be 'all' or a positive integer.");
    endif
    updates = Inf;
  elseif (! (isnumeric (updates) && isreal (updates) && isscalar (updates)
             && isfinite (updates) && updates > 0
             && updates == fix (updates)))
    error ("relieff: 'updates' must be 'all' or a positive integer.");
  endif
  if (isstring (catX) && isscalar (catX))
    catX = char (catX);
  endif
  if (! (ischar (catX) && any (strcmpi (catX, {'on', 'off'}))))
    error ("relieff: 'categoricalx' must be 'on' or 'off'.");
  endif
  catX = strcmpi (catX, 'on');
  if (isempty (sigma))
    if (isreg)
      sigma = 50;
    else
      sigma = Inf;
    endif
  elseif (! (isnumeric (sigma) && isreal (sigma) && isscalar (sigma)
             && sigma > 0))
    error ("relieff: 'sigma' must be a positive scalar.");
  endif

  ## Responses as values or class indices, with the missing ones as NaN
  if (isreg)
    y = double (Y(:));
    if (any (isinf (y)))
      error ("relieff: Y must not contain infinite values.");
    endif
  else
    [y, gn, gl] = grp2idx (Y);
  endif
  if (any (isinf (X(:))))
    error ("relieff: X must not contain infinite values.");
  endif
  keep = ! (any (isnan (X), 2) | isnan (y));

  if (! isreg)
    ## Classes present in the data, in the order of grp2idx
    held = unique (y(keep));
    [~, y] = ismember (y, held);
    gn = gn(held);
    gl = gl(held);
    nc = numel (held);
    if (ischar (prior) || (isstring (prior) && isscalar (prior)))
      prior = lower (char (prior));
      if (strcmp (prior, 'empirical'))
        prior = accumarray (y(keep), 1, [nc, 1])';
      elseif (strcmp (prior, 'uniform'))
        prior = ones (1, nc);
      else
        error ("relieff: unknown 'prior' value '%s'.", prior);
      endif
    elseif (isstruct (prior))
      if (! (isscalar (prior) && all (isfield (prior, {'group', 'prob'}))))
        error (strcat ("relieff: a 'prior' structure must have the", ...
                       " fields 'group' and 'prob'."));
      endif
      group = prior.group;
      prob = prior.prob;
      if (! (isnumeric (prob) && isreal (prob) && isvector (prob)))
        error ("relieff: 'prob' in 'prior' must be a numeric vector.");
      endif
      if ((isnumeric (group) || islogical (group))
          && (isnumeric (Y) || islogical (Y)))
        group = double (group(:));
        gl = double (gl);
      elseif (ischar (group))
        group = cellstr (group);
        gl = gn;
      elseif (iscellstr (group) || isstring (group) || iscategorical (group))
        group = cellstr (group(:));
        gl = gn;
      else
        error (strcat ("relieff: 'group' in 'prior' must name the", ...
                       " classes of Y."));
      endif
      if (numel (group) != numel (prob))
        error (strcat ("relieff: 'group' and 'prob' in 'prior' must", ...
                       " have the same length."));
      endif
      [found, at] = ismember (gl, group);
      if (! all (found))
        error ("relieff: 'prior' has no probability for class '%s'.", ...
               gn{find (! found, 1)});
      endif
      prior = prob(at);
    elseif (! (isnumeric (prior) && isreal (prior) && isvector (prior)
               && numel (prior) == nc))
      error ("relieff: 'prior' must be a vector of %d probabilities.", nc);
    endif
    prior = double (prior(:)');
    if (! (all (isfinite (prior)) && all (prior >= 0) && any (prior > 0)))
      error (strcat ("relieff: prior probabilities must be finite,", ...
                     " non-negative and not all zero."));
    endif
    prior = prior / sum (prior);
    ## The observations of a class with prior probability 0 are left out
    keep &= ismember (y, find (prior > 0));
    if (sum (prior > 0) < 2)
      error (strcat ("relieff: Y must hold at least two classes with", ...
                     " positive prior probability."));
    endif
  endif

  X = X(keep,:);
  y = y(keep);
  [n, p] = size (X);
  if (n < 2)
    error ("relieff: X and Y must hold at least two complete observations.");
  endif

  ## Differences between observations, scaled by the range of each predictor
  if (catX)
    diffs = @(r) cast (X != X(r,:), class (X));
    used = true (1, p);
  else
    w = 1 ./ (max (X, [], 1) - min (X, [], 1));
    diffs = @(r) abs (X - X(r,:)) .* w;
    used = isfinite (w);
  endif

  ## Observations whose neighbors are sought
  if (updates >= n)
    upd = 1:n;
  else
    upd = randperm (n, updates);
  endif
  m = numel (upd);

  if (isreg)
    yrange = max (y) - min (y);
    if (yrange == 0)
      error ("relieff: Y must not be constant.");
    endif
    Wdy = 0;
    Wdj = zeros (1, p, class (X));
    Wdydj = zeros (1, p, class (X));
    for r = upd
      D = diffs (r);
      cand = [1:r-1, r+1:n]';
      [~, o] = sort (sum (D(cand, used), 2));
      nb = cand(o(1:min (k, n - 1)));
      ## Normalized rank weights of the neighbors
      dw = exp (- ((1:numel (nb))' / sigma) .^ 2);
      dw /= sum (dw);
      dy = abs (y(nb) - y(r)) / yrange;
      Wdy += sum (dy .* dw);
      Wdj += sum (D(nb,:) .* dw, 1);
      Wdydj += sum (D(nb,:) .* (dy .* dw), 1);
    endfor
    weights = Wdydj / Wdy - (Wdj - Wdydj) / (m - Wdy);
  else
    weights = zeros (1, p, class (X));
    members = arrayfun (@(c) find (y == c), 1:nc, 'UniformOutput', false);
    for r = upd
      D = diffs (r);
      for c = 1:nc
        cand = members{c};
        cand(cand == r) = [];
        if (isempty (cand))
          continue;
        endif
        [~, o] = sort (sum (D(cand, used), 2));
        nb = cand(o(1:min (k, numel (cand))));
        dw = exp (- ((1:numel (nb))' / sigma) .^ 2);
        share = sum (D(nb,:) .* (dw / sum (dw)), 1) / m;
        if (c == y(r))
          weights -= share;
        else
          weights += prior(c) / (1 - prior(y(r))) * share;
        endif
      endfor
    endfor
  endif

  weights = double (weights);
  ## Most important first, a NaN weight last
  [~, idx] = sort (- weights);

endfunction

%!demo
%! ## Rank the four measurements of the iris flowers by how well they tell
%! ## the three species apart
%! load fisheriris
%! [idx, weights] = relieff (meas, species, 10);
%! bar (weights(idx));
%! xticks (1:4);
%! xticklabels (cellstr (num2str (idx')));
%! xlabel ("Predictor (column of meas), most important first");
%! ylabel ("Predictor importance weight");

## Test output
## Values checked against MATLAB R2024a
%!test
%! load fisheriris
%! [idx, w] = relieff (meas, species, 10);
%! assert_equal (idx, [4, 3, 1, 2]);
%! assert_equal (w, [0.13990740740740726, 0.12258333333333336, ...
%!                   0.35898870056497195, 0.37541666666666723], 1e-14);
%!test
%! load fisheriris
%! [~, w] = relieff (meas, species, 60);
%! assert_equal (w, [0.15579758125472401, 0.067501360544217723, ...
%!                   0.38647167070217914, 0.40577256235827669], 1e-14);
%!test
%! load fisheriris
%! y = strcmp (species(51:150), 'virginica');
%! [~, w] = relieff (meas(51:150,:), y, 10);
%! assert_equal (w, [0.045166666666666674, 0.017999999999999992, ...
%!                   0.14182051282051272, 0.2323333333333332], 1e-14);
%!test
%! load fisheriris
%! [~, w] = relieff (meas, char (species), 10);
%! assert_equal (w, [0.13990740740740731, 0.12258333333333328, ...
%!                   0.35898870056497167, 0.37541666666666712], 1e-14);
%!test
%! load fisheriris
%! [~, w] = relieff (meas, string (species), 10);
%! assert_equal (w, [0.13990740740740737, 0.12258333333333331, ...
%!                   0.35898870056497195, 0.37541666666666701], 1e-14);
%!test
%! load fisheriris
%! s = [1:53, 101:150];
%! P = struct ('group', {{'virginica', 'setosa', 'versicolor'}}, ...
%!             'prob', [0.5, 0.2, 0.3]);
%! [~, w] = relieff (meas(s,:), species(s), 5, 'prior', P);
%! assert_equal (w, [0.23459855139466784, 0.12308504602193938, ...
%!                   0.41972930722395935, 0.45912390198797975], 1e-14);
%!test
%! load fisheriris
%! [~, w] = relieff (meas(1:20,2:4), meas(1:20,1), 3);
%! assert_equal (w, [0.039264168421865836, 0.009725964602631948, ...
%!                   -0.0081765673957266277], 1e-14);
%!test
%! load carbig
%! X = [grp2idx(Mfg), grp2idx(Model), grp2idx(Origin)];
%! [idx, w] = relieff (X, MPG, 10, 'categoricalx', 'on');
%! assert_equal (idx, [2, 3, 1]);
%! assert_equal (w, [-0.0018716957152787317, 0.050096007760873751, ...
%!                   0.011404762142927625], 1e-14);
%!test
%! load carbig
%! X = [grp2idx(Mfg), grp2idx(Model), grp2idx(Origin)];
%! [~, w] = relieff (X(1:60,:), Origin(1:60,:), 5, 'categoricalx', 'on');
%! assert_equal (w, [0.73055555555555629, 0.077777777777777751, ...
%!                   1.0000000000000013], 1e-13);
%!shared X, y
%! X = [0, 0; 0.5, 0.5; 0.9, 0; 1, 1; 0.2, 0.9; 0.8, 0.3; 0.1, 0.6; 0.6, 0.9];
%! y = [1; 1; 2; 2; 1; 2; 1; 2];
%!test
%! [~, w] = relieff (X, y, 1, 'method', 'classification');
%! assert_equal (w, [0.34999999999999998, -0.16249999999999998], 1e-15);
%!test
%! [~, w] = relieff (X, y, 2, 'method', 'classification');
%! assert_equal (w, [0.29375000000000001, -0.17499999999999999], 1e-15);
%!test
%! [~, w] = relieff (X, y, 3, 'method', 'classification');
%! assert_equal (w, [0.34583333333333333, -0.21666666666666667], 1e-15);
%!test
%! [~, w] = relieff (X, y, 2, 'method', 'classification', 'sigma', 1);
%! assert_equal (w, [0.34466458926752375, -0.16368564682943915], 1e-15);
%!test
%! X = [1, 2; 1, 2; 2, 3; 3, 3; 1, 1; 2, 2; 3, 1; 1, 3];
%! [idx, w] = relieff (X, [1; 1; 1; 1; 2; 2; 2; 2] > 1, 2);
%! assert_equal (idx, [2, 1]);
%! assert_equal (w, [-0.25, 0.0625], 1e-15);
%!shared XA, yA, gA, XB, yB, gB
%! XA = [0.22665944695472717, 0.048830259591341019, 0.25411459803581238; ...
%!       0.2114872932434082, 0.77445060014724731, 0.53388506174087524; ...
%!       0.33685088157653809, 0.34994497895240784, 0.47933411598205566; ...
%!       0.056142527610063553, 0.34475591778755188, 0.27004486322402954; ...
%!       0.6719900369644165, 0.21916352212429047, 0.965262770652771; ...
%!       0.43299207091331482, 0.28173792362213135, 0.96468651294708252; ...
%!       0.59951025247573853, 0.25370854139328003, 0.90369051694869995; ...
%!       0.61985927820205688, 0.29071399569511414, 0.25213038921356201; ...
%!       0.008284449577331543, 0.14870639145374298, 0.8450893759727478; ...
%!       0.55917549133300781, 0.58620524406433105, 0.29522070288658142; ...
%!       0.57303112745285034, 0.71820592880249023, 0.43747693300247192; ...
%!       0.38736355304718018, 0.98570030927658081, 0.9433753490447998; ...
%!       0.72225898504257202, 0.32017061114311218, 0.2949681282043457; ...
%!       0.23183286190032959, 0.81152999401092529, 0.79241269826889038; ...
%!       0.085462123155593872, 0.64490348100662231, 0.016459200531244278; ...
%!       0.12298344820737839, 0.52148056030273438, 0.55653321743011475; ...
%!       0.55263209342956543, 0.87742847204208374, 0.13704623281955719; ...
%!       0.40711727738380432, 0.19621142745018005, 0.90734457969665527; ...
%!       0.61790573596954346, 0.8348698616027832, 0.016985762864351273; ...
%!       0.76110851764678955, 0.46426743268966675, 0.20533750951290131; ...
%!       0.90961301326751709, 0.54733520746231079, 0.81192338466644287; ...
%!       0.31863677501678467, 0.44641050696372986, 0.61500537395477295; ...
%!       0.23273782432079315, 0.73465442657470703, 0.74490165710449219; ...
%!       0.18065030872821808, 0.39856430888175964, 0.81240725517272949; ...
%!       0.1189243420958519, 0.20111806690692902, 0.59000658988952637; ...
%!       0.32951122522354126, 0.56738686561584473, 0.39614564180374146; ...
%!       0.91682279109954834, 0.52208811044692993, 0.78768277168273926; ...
%!       0.1428496241569519, 0.74294054508209229, 0.99268949031829834; ...
%!       0.97257781028747559, 0.54679542779922485, 0.62209063768386841; ...
%!       0.48548838496208191, 0.87875497341156006, 0.41486015915870667; ...
%!       0.39223554730415344, 0.62517315149307251, 0.44906029105186462; ...
%!       0.29200407862663269, 0.85335928201675415, 0.31995701789855957; ...
%!       0.77067995071411133, 0.64253586530685425, 0.18106210231781006; ...
%!       0.54787242412567139, 0.67721915245056152, 0.72045046091079712; ...
%!       0.12052302062511444, 0.69361585378646851, 0.053895667195320129; ...
%!       0.66175878047943115, 0.90023434162139893, 0.12591096758842468; ...
%!       0.32526770234107971, 0.26533603668212891, 0.059136319905519485; ...
%!       0.92174309492111206, 0.10128919035196304, 0.5817381739616394; ...
%!       0.96077883243560791, 0.85762202739715576, 0.90069860219955444; ...
%!       0.26504355669021606, 0.51390761137008667, 0.94491356611251831];
%! yA = [4.9158596992492676; 5.7109177112579346; 0.63560187816619873; ...
%!       0.31373497098684311; 4.1956847906112671; 1.6540294885635376; ...
%!       9.9510228633880615; 7.8833210468292236; 7.1206676959991455; ...
%!       6.3018262386322021; 7.5598502159118652; 8.6826205253601074; ...
%!       9.9508476257324219; 9.4992786645889282; 2.8919506072998047; ...
%!       3.5654160380363464; 1.446286141872406; 4.3914371728897095; ...
%!       4.0492373704910278; 1.2478844821453094; 4.1050773859024048; ...
%!       4.5266264677047729; 1.7514462769031525; 3.5887551307678223; ...
%!       3.4448394179344177; 1.9565299153327942; 1.4055664837360382; ...
%!       2.7674135565757751; 0.063248542137444019; 2.8539523482322693; ...
%!       5.6255972385406494; 9.9377048015594482; 6.2569814920425415; ...
%!       3.4007266163825989; 4.4449138641357422; 8.6494183540344238; ...
%!       9.2109048366546631; 2.5252458453178406; 0.74270054697990417; ...
%!       4.1986265778541565];
%! gA = mod ((1:40)', 3) + 1;
%! XB = [0, 8; 8, 0; 1, 2; 0, 7; 7, 1; 5, 3; 4, 4; 4, 8; 5, 7; 3, 6; ...
%!       1, 2; 3, 4; 3, 3; 4, 8; 8, 1; 2, 8; 7, 1; 7, 1; 0, 2; 7, 7; ...
%!       7, 0; 6, 4; 7, 6; 1, 3];
%! yB = [0; 16; 3; 12; 6; 0; 4; 11; 11; 5; 14; 6; 4; 5; 15; 8; 13; 3; ...
%!       13; 15; 2; 3; 4; 3];
%! gB = mod ((1:24)', 2) + 1;
%!test
%! [~, w] = relieff (XA, yA, 5);
%! assert_equal (w, [0.01388223198910471, -0.0092574193195629761, ...
%!                   -0.012105763609689901], 1e-15);
%!test
%! [~, w] = relieff (XA, yA, 5, 'sigma', 2);
%! assert_equal (w, [-0.0041728199442636477, -0.014809858819979496, ...
%!                   -0.010613408476370462], 1e-15);
%!test
%! [~, w] = relieff (XA, yA, 2);
%! assert_equal (w, [-0.0042147304959920528, -0.01845014821865984, ...
%!                   -0.0021988211394621537], 1e-15);
%!test
%! [~, w] = relieff (XA(1:6,:), yA(1:6), 10);
%! assert_equal (w, [-0.053827585946063905, 0.079324770706013403, ...
%!                   -0.13727902710184209], 1e-15);
%!test
%! [~, w] = relieff (XA, yA, 5, 'updates', 100);
%! assert_equal (w, [0.01388223198910471, -0.0092574193195630039, ...
%!                   -0.012105763609689846], 1e-15);
%!test
%! [~, w] = relieff (XA, gA, 5, 'method', 'classification');
%! assert_equal (w, [0.0021386292049166023, 0.026522833631013005, ...
%!                   -0.033770330323774092], 1e-15);
%!test
%! [~, w] = relieff (XA, gA, 5, 'method', 'classification', 'sigma', 2);
%! assert_equal (w, [0.0024473915815581577, 0.0047736034538774725, ...
%!                   -0.028887639163964078], 1e-15);
%!test
%! g = gA;
%! g(5) = 4;
%! [~, w] = relieff (XA, g, 5, 'method', 'classification');
%! assert_equal (w, [0.0085623364208488022, 0.039399610632522876, ...
%!                   -0.0025562241901765361], 1e-15);
%!test
%! [~, w] = relieff (XA, gA, 5, 'method', 'classification', ...
%!                   'prior', [0, 0.5, 0.5]);
%! assert_equal (w, [0.019765058234956511, 0.035556998164938798, ...
%!                   -0.064207251594005171], 1e-15);
%!test
%! [~, w] = relieff (XA, gA, 5, 'method', 'classification', ...
%!                   'prior', [2, 3, 5]);
%! assert_equal (w, [-0.00017399953766534648, 0.028017907734860202, ...
%!                   -0.037361123235406815], 1e-15);
%!test
%! P = struct ('group', [3, 1, 2], 'prob', [0.5, 0.2, 0.3]);
%! [~, w] = relieff (XA, gA, 5, 'method', 'classification', 'prior', P);
%! assert_equal (w, [-0.00017399953766534387, 0.028017907734860195, ...
%!                   -0.037361123235406808], 1e-15);
%!test
%! P = struct ('group', {{'3', '1', '2', '7'}}, 'prob', [0.5, 0.2, 0.3, 0.4]);
%! [~, w] = relieff (XA, cellstr (num2str (gA)), 5, 'prior', P);
%! assert_equal (w, [-0.00017399953766534387, 0.028017907734860188, ...
%!                   -0.037361123235406787], 1e-15);
%!test
%! P = struct ('group', {{'3', '1', '2'}}, 'prob', [5, 2, 3]);
%! [~, w] = relieff (XA, cellstr (num2str (gA)), 5, 'prior', P);
%! assert_equal (w, [-0.00017399953766534561, 0.028017907734860195, ...
%!                   -0.037361123235406801], 1e-15);
%!test
%! [~, w] = relieff (XA, cellstr (num2str (gA)), 5, 'prior', [0.2; 0.3; 0.5]);
%! assert_equal (w, [0.002374388536664899, 0.029954823922265453, ...
%!                   -0.028297242346431235], 1e-15);
%!test
%! [~, w] = relieff (XA, categorical (gA, 1:4), 5);
%! assert_equal (w, [0.002138629204916604, 0.026522833631012988, ...
%!                   -0.033770330323774085], 1e-15);
%!test
%! [~, w] = relieff (XA, categorical (gA, 1:4), 5, 'prior', 'uniform');
%! assert_equal (w, [0.0020399923619240142, 0.02693029952279568, ...
%!                   -0.033439605784320946], 1e-15);
%!test
%! [~, w] = relieff (XA, categorical (gA, 1:4), 5, 'prior', [0.2, 0.3, 0.5]);
%! assert_equal (w, [-0.00017399953766534474, 0.028017907734860192, ...
%!                   -0.037361123235406801], 1e-15);
%!test
%! X = XA;
%! X([3, 7], 2) = NaN;
%! y = yA;
%! y(10) = NaN;
%! [~, w] = relieff (X, y, 5);
%! assert_equal (w, [0.0096723878723309897, -0.0050429741167224706, ...
%!                   -0.014019835775775619], 1e-15);
%!test
%! X = XA;
%! X([3, 7], 2) = NaN;
%! y = cellstr (num2str (gA));
%! y{12} = '';
%! [~, w] = relieff (X, y, 5);
%! assert_equal (w, [-0.0048661229809348686, 0.045375869235758988, ...
%!                   -0.019229982110778903], 1e-15);
%!test
%! [idx, w] = relieff ([XA(:,1:2), 3 * ones(40, 1)], yA, 5);
%! assert_equal (idx, [1, 2, 3]);
%! assert_equal (w, [0.0060575278856132564, -0.011272620704472172, NaN], ...
%!               1e-15);
%!test
%! [idx, w] = relieff ([XA(:,1:2), 3 * ones(40, 1)], gA, 5, ...
%!                     'method', 'classification');
%! assert_equal (idx, [1, 2, 3]);
%! assert_equal (w, [0.015455906238505155, -0.0058458986874681575, NaN], ...
%!               1e-15);
%!test
%! [idx, w] = relieff (single (XA), yA, 5);
%! assert_equal (class (idx), 'double');
%! assert_equal (class (w), 'double');
%! assert_equal (w, [0.013882221976705744, -0.0092574315350100744, ...
%!                   -0.012105737249683085], 1e-6);
%!test
%! [~, w] = relieff (single (XA), gA, 5, 'method', 'classification');
%! assert_equal (w, [0.0021386349108070135, 0.026522833853960037, ...
%!                   -0.033770337700843811], 1e-6);
%!test
%! [~, w] = relieff (XB, yB, 4);
%! assert_equal (w, [-0.015833802328636792, -0.017033417536778783], 1e-15);
%!test
%! [~, w] = relieff (XB, yB, 4, 'sigma', 1);
%! assert_equal (w, [-0.0062979535245773766, -0.026549633818311213], 1e-15);
%!test
%! [~, w] = relieff (XB, yB, 4, 'sigma', Inf);
%! assert_equal (w, [-0.015822377582084102, -0.017017026629918081], 1e-15);
%!test
%! [~, w] = relieff (XB, gB, 3, 'method', 'classification');
%! assert_equal (w, [-0.020833333333333318, -0.012152777777777783], 1e-15);
%!test
%! [~, w] = relieff (XB, gB, 3, 'method', 'classification', 'sigma', 1);
%! assert_equal (w, [-0.064720206045157824, -0.0019737829821077271], 1e-15);
%!test
%! [~, w] = relieff (XB, yB, 4, 'categoricalx', 'on');
%! assert_equal (w, [-0.0021677478940700379, 0.036294175269704521], 1e-15);
%!test
%! [~, w] = relieff (XB, yB, 4, 'categoricalx', 'on', 'sigma', 1);
%! assert_equal (w, [-0.095389338861245093, -0.003446685860238996], 1e-15);
%!test
%! [~, w] = relieff (XB, gB, 3, 'method', 'classification', ...
%!                   'categoricalx', 'on');
%! assert_equal (w, [-0.1666666666666666, -0.013888888888888895], 1e-15);
%!test
%! rand ('seed', 1);
%! [idx, w] = relieff (XA, yA, 5, 'updates', 10);
%! assert_equal (sort (idx), [1, 2, 3]);
%! assert_equal (all (isfinite (w)), true);

## Test input validation
%!error <Invalid call to relieff> relieff (XA, yA)
%!error <relieff: X must be a real matrix of type double or single.> ...
%! relieff (int8 (XB), yB, 4)
%!error <relieff: X must be a real matrix of type double or single.> ...
%! relieff (XA + 1i, yA, 5)
%!error <relieff: Y must be a vector of responses or class names.> ...
%! relieff (XA, [yA, yA], 5)
%!error <relieff: X must have as many rows as Y has observations.> ...
%! relieff (XA, yA(1:30), 5)
%!error <relieff: K must be a positive integer scalar.> relieff (XA, yA, 0)
%!error <relieff: K must be a positive integer scalar.> relieff (XA, yA, 1.5)
%!error <relieff: K must be a positive integer scalar.> relieff (XA, yA, [5, 6])
%!error <relieff: unknown parameter name 'foo'.> relieff (XA, yA, 5, 'foo', 1)
%!error <relieff: optional arguments must be Name-Value pairs.> ...
%! relieff (XA, yA, 5, 'method')
%!error <relieff: 'method' must be a character vector.> ...
%! relieff (XA, yA, 5, 'method', 1)
%!error <relieff: unknown method 'foo'.> relieff (XA, yA, 5, 'method', 'foo')
%!error <relieff: regression needs a numeric Y.> ...
%! relieff (XA, cellstr (num2str (gA)), 5, 'method', 'regression')
%!error <relieff: 'prior' applies to classification only.> ...
%! relieff (XA, yA, 5, 'prior', 'uniform')
%!error <relieff: 'updates' must be 'all' or a positive integer.> ...
%! relieff (XA, yA, 5, 'updates', 0)
%!error <relieff: 'updates' must be 'all' or a positive integer.> ...
%! relieff (XA, yA, 5, 'updates', 1.5)
%!error <relieff: 'updates' must be 'all' or a positive integer.> ...
%! relieff (XA, yA, 5, 'updates', 'some')
%!error <relieff: 'categoricalx' must be 'on' or 'off'.> ...
%! relieff (XA, yA, 5, 'categoricalx', 'yes')
%!error <relieff: 'categoricalx' must be 'on' or 'off'.> ...
%! relieff (XA, yA, 5, 'categoricalx', true)
%!error <relieff: 'sigma' must be a positive scalar.> ...
%! relieff (XA, yA, 5, 'sigma', 0)
%!error <relieff: 'sigma' must be a positive scalar.> ...
%! relieff (XA, yA, 5, 'sigma', [1, 2])
%!error <relieff: Y must not contain infinite values.> ...
%! relieff (XA, [Inf; yA(2:end)], 5)
%!error <relieff: X must not contain infinite values.> ...
%! relieff ([Inf, 0, 0; XA(2:end,:)], yA, 5)
%!error <relieff: unknown 'prior' value 'foo'.> ...
%! relieff (XA, gA, 5, 'method', 'classification', 'prior', 'foo')
%!error <relieff: a 'prior' structure must have the fields 'group' and 'prob'.> ...
%! relieff (XA, gA, 5, 'method', 'classification', 'prior', struct ('a', 1))
%!error <relieff: 'prob' in 'prior' must be a numeric vector.> ...
%! relieff (XA, gA, 5, 'method', 'classification', 'prior', ...
%!          struct ('group', [1, 2, 3], 'prob', 'abc'))
%!error <relieff: 'group' in 'prior' must name the classes of Y.> ...
%! relieff (XA, cellstr (num2str (gA)), 5, 'prior', ...
%!          struct ('group', [1, 2, 3], 'prob', [1, 1, 1]))
%!error <relieff: 'group' and 'prob' in 'prior' must have the same length.> ...
%! relieff (XA, gA, 5, 'method', 'classification', 'prior', ...
%!          struct ('group', [1, 2, 3], 'prob', [1, 1]))
%!error <relieff: 'prior' has no probability for class '2'.> ...
%! relieff (XA, cellstr (num2str (gA)), 5, 'prior', ...
%!          struct ('group', {{'3', '1'}}, 'prob', [0.5, 0.2]))
%!error <relieff: 'prior' must be a vector of 3 probabilities.> ...
%! relieff (XA, gA, 5, 'method', 'classification', 'prior', [0.5, 0.5])
%!error <relieff: prior probabilities must be finite, non-negative and not all zero.> ...
%! relieff (XA, gA, 5, 'method', 'classification', 'prior', [0.2, 0.3, NaN])
%!error <relieff: prior probabilities must be finite, non-negative and not all zero.> ...
%! relieff (XA, gA, 5, 'method', 'classification', 'prior', [-1, 1, 1])
%!error <relieff: prior probabilities must be finite, non-negative and not all zero.> ...
%! relieff (XA, gA, 5, 'method', 'classification', 'prior', [0, 0, 0])
%!error <relieff: Y must hold at least two classes with positive prior probability.> ...
%! relieff (XA, gA, 5, 'method', 'classification', 'prior', [1, 0, 0])
%!error <relieff: Y must hold at least two classes with positive prior probability.> ...
%! relieff (XA, ones (40, 1), 5, 'method', 'classification')
%!error <relieff: X and Y must hold at least two complete observations.> ...
%! relieff (XA(1,:), yA(1), 5)
%!error <relieff: Y must not be constant.> relieff (XA, ones (40, 1), 5)
