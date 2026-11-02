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
## @deftypefn  {statistics} {@var{idx} =} fsrftest (@var{Tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{idx} =} fsrftest (@var{Tbl}, @var{formula})
## @deftypefnx {statistics} {@var{idx} =} fsrftest (@var{Tbl}, @var{Y})
## @deftypefnx {statistics} {@var{idx} =} fsrftest (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{idx} =} fsrftest (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{idx}, @var{scores}] =} fsrftest (@dots{})
##
## Rank the predictors for regression with F-tests.
##
## @code{@var{idx} = fsrftest (@var{X}, @var{Y})} tests each predictor, a
## column of the matrix @var{X}, with a one-way analysis of variance: whether
## the response @var{Y} has the same mean in every group of observations the
## predictor forms.  It ranks the predictors by the result.  @var{idx} is a row
## vector of column indices, from the most important predictor to the least;
## predictors of equal score keep their order.
##
## @var{X} is a real matrix of type double or single with one row per
## observation, and @var{Y} a real vector of type double or single with one
## response per row of @var{X}.
##
## @code{@var{idx} = fsrftest (@var{Tbl}, @var{ResponseVarName})} takes the
## predictors and the response from the table @var{Tbl}: @var{ResponseVarName}
## names the variable holding the response and every other variable is a
## predictor.  @code{@var{idx} = fsrftest (@var{Tbl}, @var{formula})} names
## them with a model formula instead, @qcode{'Y ~ x1 + x2'}, which holds main
## effects only, and ranks only the predictors it names, @var{idx} indexing
## them in the order the table holds them.  @code{@var{idx} = fsrftest
## (@var{Tbl}, @var{Y})} takes the predictors from the table and the response
## from @var{Y}.  A logical, categorical, text or string variable of the table
## is a categorical predictor; an ordinal categorical one is continuous, its
## categories in their order.
##
## @code{[@var{idx}, @var{scores}] = fsrftest (@dots{})} also returns the
## score of each predictor, a row vector of type double in the order of the
## predictors, so that @code{@var{scores}(@var{idx})} is in decreasing order.
## The score is @math{-log (p)}, @math{p} being the p-value of the predictor's
## F-test, so a larger score marks a more important predictor; it is
## @code{Inf} where @math{p} is 0.  A predictor that forms a single group, or
## as many groups as observations, scores 0, and so does every predictor when
## the response is constant.
##
## A continuous predictor is first binned.  Its @math{n} known values are
## sorted, the upper edge of bin @math{k} is the value at position
## @code{ceil (@var{k} * @var{n} / @var{NumBins})}, and a value goes to the
## first bin whose edge it does not exceed.  The bins so hold equal numbers of
## observations as nearly as the ties allow, the first ones taking any
## remainder, and equal values always share a bin, which may leave fewer bins
## than asked for.  A categorical predictor takes its categories as they
## are.  The F-test weighs each observation by its weight, scaling the
## weights not mattering.
##
## Each predictor is tested on the observations it is known for.  An
## observation whose response is missing is left out of every test.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.3 0.7
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'CategoricalPredictors'} @tab The categorical predictors: a
## vector of their indices, a logical vector with one value per predictor, or
## @qcode{'all'}.  With a table they are added to those the table makes
## categorical, and they may also be named.  The default is none for a matrix.
##
## @item @qcode{'NumBins'} @tab The number of bins of a continuous predictor, a
## positive integer.  The default is 10.
##
## @item @qcode{'UseMissing'} @tab @qcode{false}, the default, to test each
## predictor on the observations it is known for, or @qcode{true} to test it
## on every observation, the missing values of a categorical predictor forming
## one more category and those of a continuous predictor one more bin.  MATLAB
## R2024a and R2026b document this option but ignore it in @code{fsrftest},
## scoring every predictor alike either way; Octave does what they document.
##
## @item @qcode{'Weights'} @tab The observation weights, a nonnegative finite
## vector of type double or single with one value per row, or with a table the
## name of one of its variables, which is then no predictor.  The default is
## equal weights.
## @end multitable
##
## The option names are not case sensitive but must be given in full.  With a
## matrix, the categorical predictors cannot be named, the predictors having no
## names, as in MATLAB.  The test is computed in double precision whatever the
## type of the data; MATLAB R2024a gives a different score for a response of
## type single than for the same values as double, which Octave does not.
##
## @seealso{fscchi2, relieff, sequentialfs, anova1}
## @end deftypefn

function [idx, scores] = fsrftest (X, Y, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  ## Every option is named; a table brings its own predictor names, which the
  ## call may not set
  optNames = {'CategoricalPredictors', 'NumBins', 'UseMissing', 'Weights'};
  if (mod (numel (varargin), 2) != 0)
    error ("fsrftest: optional arguments must be Name-Value pairs.");
  endif
  for k = 1:2:numel (varargin)
    if (! (ischar (varargin{k}) && any (strcmpi (varargin{k}, optNames))))
      if (ischar (varargin{k}))
        error ("fsrftest: unknown parameter name '%s'.", varargin{k});
      endif
      error ("fsrftest: optional arguments must be Name-Value pairs.");
    endif
  endfor

  ## A table is resolved into the predictor matrix, its categorical columns
  ## and the response; an ordinal column is continuous, as in MATLAB
  catGiven = [];
  if (istable (X))
    tnames = X.Properties.VariableNames;
    [X, Y, args, ~, errmsg] = tableFrame (X, Y, varargin, true);
    if (! isempty (errmsg))
      error ("fsrftest: %s", errmsg);
    endif
    [pnames, ~, catGiven, args] = parsePairedArguments ( ...
                   {'PredictorNames', 'ResponseName', ...
                    'CategoricalPredictors'}, {[], [], []}, args);
    ## The predictors are ranked in the order the table holds them, a
    ## formula naming only which, as in MATLAB
    [~, pos] = ismember (pnames, tnames);
    [~, order] = sort (pos);
    X = X(:,order);
    iscatT = false (1, numel (pnames));
    iscatT(catGiven) = true;
    catGiven = find (iscatT(order));
  else
    args = varargin;
  endif
  dfValues = {[], 10, false, []};
  [CatPreds, NumBins, UseMissing, Weights] = ...
                   parsePairedArguments (optNames, dfValues, args);

  ## Predictors and response
  if (! (isfloat (X) && isreal (X) && ismatrix (X) && ! isempty (X)))
    error (strcat ("fsrftest: X must be a nonempty real matrix of type", ...
                   " double or single."));
  endif
  X = double (X);
  [n, p] = size (X);
  if (! (isfloat (Y) && isreal (Y) && isvector (Y)))
    error ("fsrftest: Y must be a real vector of type double or single.");
  endif
  if (numel (Y) != n)
    error ("fsrftest: X and Y must have the same number of observations.");
  endif
  Y = double (Y(:));

  ## Options
  if (! (isnumeric (NumBins) && isreal (NumBins) && isscalar (NumBins)
         && NumBins >= 1 && NumBins == fix (NumBins)))
    error ("fsrftest: 'NumBins' must be a positive integer.");
  endif
  if (! ((islogical (UseMissing) || isnumeric (UseMissing))
         && isscalar (UseMissing) && any (UseMissing == [0, 1])))
    error ("fsrftest: 'UseMissing' must be a logical scalar.");
  endif
  UseMissing = logical (UseMissing);
  if (isempty (Weights))
    Weights = ones (n, 1);
  else
    errmsg = weightsClass (Weights);
    if (! isempty (errmsg))
      error ("fsrftest: %s", errmsg);
    endif
    if (! (isvector (Weights) && numel (Weights) == n))
      error ("fsrftest: 'Weights' must have one value per observation.");
    endif
    Weights = double (Weights(:));
    if (! (all (isfinite (Weights)) && all (Weights >= 0)
           && sum (Weights) > 0))
      error (strcat ("fsrftest: 'Weights' must be nonnegative, finite and", ...
                     " not all zero."));
    endif
  endif

  ## The categorical predictors: those a table gave and those the call named,
  ## the call naming them by index, by a logical vector or as 'all'
  iscat = false (1, p);
  iscat(catGiven) = true;
  if (! isempty (CatPreds))
    if (ischar (CatPreds) && isrow (CatPreds) && strcmpi (CatPreds, 'all'))
      iscat(:) = true;
    elseif (islogical (CatPreds) && isvector (CatPreds))
      if (numel (CatPreds) != p)
        error (strcat ("fsrftest: 'CategoricalPredictors' as a logical", ...
                       " vector must have one value per predictor."));
      endif
      iscat |= CatPreds(:)';
    elseif (isnumeric (CatPreds) && isvector (CatPreds))
      if (! (isreal (CatPreds) && all (CatPreds == fix (CatPreds))
             && all (CatPreds >= 1) && all (CatPreds <= p)))
        error (strcat ("fsrftest: 'CategoricalPredictors' must index the", ...
                       " predictors."));
      endif
      iscat(CatPreds) = true;
    else
      error (strcat ("fsrftest: 'CategoricalPredictors' must be indices, a", ...
                     " logical vector or 'all'; naming the predictors", ...
                     " needs a table."));
    endif
  endif

  ## A missing response leaves its observation out
  used = ! isnan (Y);
  X = X(used,:);
  Y = Y(used);
  w = Weights(used);

  ## One test per predictor
  scores = zeros (1, p);
  for j = 1:p
    x = X(:,j);
    known = ! isnan (x);
    if (UseMissing)
      inTest = true (rows (x), 1);
    else
      inTest = known;
    endif
    if (! any (known & inTest))
      continue;
    endif

    ## Bins of a continuous predictor, bounded by the values at equal steps
    ## through the sorted values, a value going to the first bin whose upper
    ## edge it does not exceed; categories of a categorical one.  The missing
    ## values, when they are used, are one more group.
    xk = x(known);
    if (iscat(j))
      [~, ~, g] = unique (xk);
    else
      sx = sort (xk);
      edges = sx(ceil ((1:NumBins-1) * numel (xk) / NumBins));
      g = 1 + sum (xk > edges(:)', 2);
    endif
    group = zeros (rows (x), 1);
    group(known) = g;
    group(! known) = max (g) + 1;

    ## The weighted one-way analysis of variance of the response across the
    ## groups, the empty ones left out
    gt = group(inTest);
    yt = Y(inTest);
    wt = w(inTest);
    keep = wt > 0;
    gt = gt(keep);
    yt = yt(keep);
    wt = wt(keep);
    [~, ~, gt] = unique (gt);
    K = max ([0; gt]);
    nt = numel (yt);
    if (K < 2 || nt - K < 1)
      continue;
    endif
    m = sum (wt .* yt) / sum (wt);
    gw = accumarray (gt, wt);
    gm = accumarray (gt, wt .* yt) ./ gw;
    ssb = sum (gw .* (gm - m) .^ 2);
    ssw = sum (wt .* (yt - gm(gt)) .^ 2);
    if (! (ssb > 0))
      continue;
    endif
    F = (ssb / (K - 1)) / (ssw / (nt - K));
    scores(j) = -log (fcdf (F, K - 1, nt - K, 'upper'));
  endfor

  [~, idx] = sort (scores, 'descend');

endfunction

## Test output: the binning, measured on R2024a
%!shared xb, yb, w
%! xb = ((1:40) .^ 2)';
%! yb = mod ((1:40)', 7) + (1:40)' / 10;
%! w = 1 + mod ((1:40)', 3);
%!test  # ten bins, the default
%! [idx, s] = fsrftest (xb, yb);
%! assert_equal (idx, 1);
%! assert_equal (s, 3.28972104684277, 1e-13);
%!test  # the remainder goes to the first bins
%! [~, s] = fsrftest (xb, yb, 'NumBins', 3);
%! assert_equal (s, 3.96652214599653, 1e-13);
%!test
%! [~, s] = fsrftest (xb, yb, 'NumBins', 7);
%! assert_equal (s, 2.73641811567188, 1e-13);
%!test  # ties share a bin
%! xt = [ones(15, 1); 2 * ones(10, 1); (3:17)'];
%! [~, s] = fsrftest (xt, yb);
%! assert_equal (s, 3.63139675868502, 1e-13);
%!test  # equal values share a bin, here the two smallest
%! [~, s] = fsrftest (mod (xb, 7), yb, 'NumBins', 4);
%! assert_equal (s, 1.65333228183337, 1e-13);
%!test  # a categorical predictor takes its values as they are
%! [~, s] = fsrftest ([xb, mod((1:40)', 4)], yb, 'CategoricalPredictors', 2);
%! assert_equal (s, [3.28972104684277, 0.0395057912929183], 1e-13);
%! [~, s] = fsrftest ([xb, mod((1:40)', 4)], yb, ...
%!                    'CategoricalPredictors', [false, true]);
%! assert_equal (s, [3.28972104684277, 0.0395057912929183], 1e-13);

## Weights and missing values, measured on R2024a
%!test  # weights, whose scale does not matter
%! [~, s] = fsrftest (xb, yb, 'Weights', w);
%! assert_equal (s, 3.40381130018073, 1e-13);
%! assert_equal (nthargout (2, @fsrftest, xb, yb, 'Weights', 3 * w), s, 1e-14);
%!test  # a predictor is tested on the observations it is known for
%! xn = xb;
%! xn([3, 17, 29]) = NaN;
%! [~, s] = fsrftest (xn, yb);
%! assert_equal (s, 3.39160408073151, 1e-13);
%! [~, s] = fsrftest (xn, yb, 'Weights', w);
%! assert_equal (s, 4.6530885484303, 1e-13);
%!test  # an observation of missing response is left out
%! yn = yb;
%! yn([5, 11]) = NaN;
%! [~, s] = fsrftest (xb, yn);
%! assert_equal (s, 2.46318862209335, 1e-13);
%! [~, s] = fsrftest (xb, yn, 'Weights', w);
%! assert_equal (s, 3.68126942373109, 1e-13);

## Octave tests a predictor's missing values as one more group when asked,
## as MATLAB documents; R2024a and R2026b ignore the option
%!test  # three groups: F = 64 on 2 and 3 degrees of freedom
%! x = [1; 1; 2; 2; NaN; NaN];
%! y = [1; 2; 5; 6; 9; 10];
%! [~, s] = fsrftest (x, y, 'CategoricalPredictors', 1, 'UseMissing', true);
%! assert_equal (s, -log (fcdf (64, 2, 3, 'upper')), 1e-14);
%!test  # without it, two: F = 32 on 1 and 2
%! x = [1; 1; 2; 2; NaN; NaN];
%! y = [1; 2; 5; 6; 9; 10];
%! [~, s] = fsrftest (x, y, 'CategoricalPredictors', 1);
%! assert_equal (s, -log (fcdf (32, 1, 2, 'upper')), 1e-14);
%!test  # a continuous predictor's missing values form one more bin
%! x = [1; 2; 3; 4; NaN; NaN];
%! y = [1; 2; 5; 6; 9; 10];
%! [~, s] = fsrftest (x, y, 'NumBins', 2, 'UseMissing', true);
%! assert_equal (s, -log (fcdf (64, 2, 3, 'upper')), 1e-14);

## Scores of 0 and Inf, measured on R2024a
%!test  # ranking keeps the order of equal scores
%! [idx, s] = fsrftest ([xb, ones(40, 1), NaN(40, 1), xb], yb);
%! assert_equal (idx, [1, 4, 2, 3]);
%! assert_equal (s, [3.28972104684277, 0, 0, 3.28972104684277], 1e-13);
%!test  # a single bin
%! [~, s] = fsrftest (xb, yb, 'NumBins', 1);
%! assert_equal (s, 0);
%!test  # as many groups as observations
%! [~, s] = fsrftest (xb, yb, 'NumBins', 100);
%! assert_equal (s, 0);
%!test  # a constant response
%! [~, s] = fsrftest (xb, ones (40, 1));
%! assert_equal (s, 0);
%!test  # groups that separate the response wholly
%! x = (1:400)';
%! [~, s] = fsrftest (x, [zeros(200, 1); ones(200, 1)] + mod (x, 3) / 100);
%! assert_equal (s, Inf);
%!test  # single predictors are tested in double
%! [~, s] = fsrftest (single (xb), yb);
%! assert_equal (class (s), 'double');
%! assert_equal (s, 3.28972104684277, 1e-13);

## A table, measured on R2024a
%!test  # a categorical variable with an undefined value
%! k = categorical (repmat ({'u'; 'v'; 'w'; 'u'}, 10, 1));
%! k([2, 9, 30]) = missing;
%! [~, s] = fsrftest (table (k, xb, yb), 'yb');
%! assert_equal (s, [0.189033357981702, 3.28972104684277], 1e-13);
%!test  # the response given beside the table
%! k = categorical (repmat ({'u'; 'v'; 'w'; 'u'}, 10, 1));
%! k([2, 9, 30]) = missing;
%! [~, s] = fsrftest (table (k, xb), yb);
%! assert_equal (s, [0.189033357981702, 3.28972104684277], 1e-13);
%!test  # a formula ranks the predictors it names, in the table's order
%! k = categorical (repmat ({'u'; 'v'; 'w'; 'u'}, 10, 1));
%! k([2, 9, 30]) = missing;
%! [idx, s] = fsrftest (table (xb, k, yb), 'yb ~ k + xb');
%! assert_equal (idx, [1, 2]);
%! assert_equal (s, [3.28972104684277, 0.189033357981702], 1e-13);
%!test  # weights named in the table are no predictor
%! x1 = (1:12)';
%! wt = [1; 1; 2; 2; 1; 1; 3; 3; 1; 1; 2; 2];
%! y12 = mod ((1:12)', 5) + (1:12)' / 7;
%! [idx, s] = fsrftest (table (x1, wt, y12), 'y12', 'Weights', 'wt');
%! assert_equal (idx, 1);
%! assert_equal (s, 1.0443158836381, 1e-13);

## Test input validation
%!error<Invalid call to fsrftest> fsrftest (1)
%!error<fsrftest: optional arguments must be Name-Value pairs.> ...
%! fsrftest ([1; 2], [1; 2], 'NumBins')
%!error<fsrftest: unknown parameter name 'Prior'.> ...
%! fsrftest ([1; 2], [1; 2], 'Prior', 'uniform')
%!error<fsrftest: X must be a nonempty real matrix of type double or single.> ...
%! fsrftest (int8 ([1; 2]), [1; 2])
%!error<fsrftest: Y must be a real vector of type double or single.> ...
%! fsrftest ([1; 2], [true; false])
%!error<fsrftest: Y must be a real vector of type double or single.> ...
%! fsrftest ([1; 2], {'a'; 'b'})
%!error<fsrftest: X and Y must have the same number of observations.> ...
%! fsrftest ([1; 2; 3], [1; 2])
%!error<fsrftest: 'NumBins' must be a positive integer.> ...
%! fsrftest ([1; 2], [1; 2], 'NumBins', 0)
%!error<fsrftest: 'UseMissing' must be a logical scalar.> ...
%! fsrftest ([1; 2], [1; 2], 'UseMissing', 2)
%!error<fsrftest: 'Weights' must have one value per observation.> ...
%! fsrftest ([1; 2], [1; 2], 'Weights', [1; 1; 1])
%!error<fsrftest: 'Weights' must be nonnegative, finite and not all zero.> ...
%! fsrftest ([1; 2], [1; 2], 'Weights', [1; -1])
%!error<fsrftest: 'CategoricalPredictors' as a logical vector must have one value per predictor.> ...
%! fsrftest ([1; 2], [1; 2], 'CategoricalPredictors', [true, false])
%!error<fsrftest: 'CategoricalPredictors' must index the predictors.> ...
%! fsrftest ([1; 2], [1; 2], 'CategoricalPredictors', 2)
%!error<fsrftest: 'CategoricalPredictors' must be indices, a logical vector or 'all'; naming the predictors needs a table.> ...
%! fsrftest ([1; 2], [1; 2], 'CategoricalPredictors', 'x1')

