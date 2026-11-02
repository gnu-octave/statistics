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
## @deftypefn  {statistics} {@var{idx} =} fscchi2 (@var{Tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{idx} =} fscchi2 (@var{Tbl}, @var{formula})
## @deftypefnx {statistics} {@var{idx} =} fscchi2 (@var{Tbl}, @var{Y})
## @deftypefnx {statistics} {@var{idx} =} fscchi2 (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{idx} =} fscchi2 (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{idx}, @var{scores}] =} fscchi2 (@dots{})
##
## Rank the predictors for classification with chi-square tests.
##
## @code{@var{idx} = fscchi2 (@var{X}, @var{Y})} tests each predictor, a
## column of the matrix @var{X}, for independence from the class labels
## @var{Y} with a chi-square test, and ranks the predictors by the result.
## @var{idx} is a row vector of column indices, from the most important
## predictor to the least; predictors of equal score keep their order.
##
## @var{X} is a real matrix of type double or single with one row per
## observation.  @var{Y} holds one class label per row of @var{X}: a numeric
## or logical vector, a categorical vector, a cell array of character vectors,
## a string array, or a character matrix with one row per observation.
##
## @code{@var{idx} = fscchi2 (@var{Tbl}, @var{ResponseVarName})} takes the
## predictors and the class labels from the table @var{Tbl}:
## @var{ResponseVarName} names the variable holding the labels and every other
## variable is a predictor.  @code{@var{idx} = fscchi2 (@var{Tbl},
## @var{formula})} names them with a model formula instead,
## @qcode{'Y ~ x1 + x2'}, which holds main effects only, and ranks only the
## predictors it names, @var{idx} indexing them in the order the table holds
## them.  @code{@var{idx} = fscchi2 (@var{Tbl}, @var{Y})} takes the predictors
## from the table and the labels from @var{Y}.  A logical, categorical, text
## or string variable of the table is a categorical predictor; an ordinal
## categorical one is continuous, its categories in their order.
##
## @code{[@var{idx}, @var{scores}] = fscchi2 (@dots{})} also returns the score
## of each predictor, a row vector in the order of the predictors, so that
## @code{@var{scores}(@var{idx})} is in decreasing order.  The score is
## @math{-log (p)}, @math{p} being the p-value of the predictor's test, so a
## larger score marks a more important predictor; it is @code{Inf} where
## @math{p} is 0.  A predictor that cannot separate the classes, being
## constant or wholly missing, scores 0, and so does every predictor when the
## labels hold a single class.
##
## A continuous predictor is first binned.  Its @math{n} known values are
## sorted, the upper edge of bin @math{k} is the value at position
## @code{ceil (@var{k} * @var{n} / @var{NumBins})}, and a value goes to the
## first bin whose edge it does not exceed.  The bins so hold equal numbers of
## observations as nearly as the ties allow, the first ones taking any
## remainder, and equal values always share a bin, which may leave fewer bins
## than asked for.  A categorical predictor takes its categories as they
## are.  The test is Pearson's chi-square test of
## independence on the table of weighted counts of the bins or categories
## against the classes, the empty ones left out.
##
## Each predictor is tested on the observations it is known for.  Their
## weights are scaled within each class so that the class carries its prior
## probability, and over all so that they sum to the number of those
## observations.  The empirical prior, the default, is the weighted share of
## each class among those observations, which leaves the weights as given.  An
## observation whose class label is missing is left out of every test.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.25 0.75
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'CategoricalPredictors'} @tab The categorical predictors: a
## vector of their indices, a logical vector with one value per predictor, or
## @qcode{'all'}.  With a table they are added to those the table makes
## categorical, and they may also be named.  The default is none for a matrix.
##
## @item @qcode{'ClassNames'} @tab The classes to use, in the type of @var{Y}.
## Observations of other classes are left out, and the order of the names is
## the order of a @qcode{'Prior'} vector.  The default is every class of
## @var{Y}, sorted.
##
## @item @qcode{'NumBins'} @tab The number of bins of a continuous predictor, a
## positive integer.  The default is 10.
##
## @item @qcode{'Prior'} @tab The prior probabilities of the classes:
## @qcode{'empirical'}, the default; @qcode{'uniform'}; a vector with one
## nonnegative value per class, in the order of the classes; or a structure
## with the fields @qcode{ClassNames} and @qcode{ClassProbs}.  The values are
## normalized to sum to 1.
##
## @item @qcode{'UseMissing'} @tab @qcode{false}, the default, to test each
## predictor on the observations it is known for, or @qcode{true} to test it
## on every observation, the missing values of a categorical predictor forming
## one more category and those of a continuous predictor one more bin.  MATLAB
## R2024a and R2026b document this option but ignore it in @code{fscchi2},
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
## names, as in MATLAB.
##
## @seealso{relieff, sequentialfs, crosstab}
## @end deftypefn

function [idx, scores] = fscchi2 (X, Y, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  ## Every option is named; a table brings its own predictor names, which the
  ## call may not set
  optNames = {'CategoricalPredictors', 'ClassNames', 'NumBins', 'Prior', ...
              'UseMissing', 'Weights'};
  if (mod (numel (varargin), 2) != 0)
    error ("fscchi2: optional arguments must be Name-Value pairs.");
  endif
  for k = 1:2:numel (varargin)
    if (! (ischar (varargin{k}) && any (strcmpi (varargin{k}, optNames))))
      if (ischar (varargin{k}))
        error ("fscchi2: unknown parameter name '%s'.", varargin{k});
      endif
      error ("fscchi2: optional arguments must be Name-Value pairs.");
    endif
  endfor

  ## A table is resolved into the predictor matrix, its categorical columns
  ## and the response; an ordinal column is continuous, as in MATLAB
  catGiven = [];
  if (istable (X))
    tnames = X.Properties.VariableNames;
    [X, Y, args, ~, errmsg] = tableFrame (X, Y, varargin, true);
    if (! isempty (errmsg))
      error ("fscchi2: %s", errmsg);
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
  dfValues = {[], [], 10, 'empirical', false, []};
  [CatPreds, ClassNames, NumBins, Prior, UseMissing, Weights] = ...
                   parsePairedArguments (optNames, dfValues, args);

  ## Predictors and response
  if (! (isfloat (X) && isreal (X) && ismatrix (X) && ! isempty (X)))
    error (strcat ("fscchi2: X must be a nonempty real matrix of type", ...
                   " double or single."));
  endif
  X = double (X);
  [n, p] = size (X);
  if (isempty (Y) || rows (Y) != n || (! ischar (Y) && ! isvector (Y)))
    error ("fscchi2: X and Y must have the same number of observations.");
  endif

  ## Options
  if (! (isnumeric (NumBins) && isreal (NumBins) && isscalar (NumBins)
         && NumBins >= 1 && NumBins == fix (NumBins)))
    error ("fscchi2: 'NumBins' must be a positive integer.");
  endif
  if (! ((islogical (UseMissing) || isnumeric (UseMissing))
         && isscalar (UseMissing) && any (UseMissing == [0, 1])))
    error ("fscchi2: 'UseMissing' must be a logical scalar.");
  endif
  UseMissing = logical (UseMissing);
  if (isempty (Weights))
    Weights = ones (n, 1);
  else
    errmsg = weightsClass (Weights);
    if (! isempty (errmsg))
      error ("fscchi2: %s", errmsg);
    endif
    if (! (isvector (Weights) && numel (Weights) == n))
      error ("fscchi2: 'Weights' must have one value per observation.");
    endif
    Weights = double (Weights(:));
    if (! (all (isfinite (Weights)) && all (Weights >= 0)
           && sum (Weights) > 0))
      error (strcat ("fscchi2: 'Weights' must be nonnegative, finite and", ...
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
        error (strcat ("fscchi2: 'CategoricalPredictors' as a logical", ...
                       " vector must have one value per predictor."));
      endif
      iscat |= CatPreds(:)';
    elseif (isnumeric (CatPreds) && isvector (CatPreds))
      if (! (isreal (CatPreds) && all (CatPreds == fix (CatPreds))
             && all (CatPreds >= 1) && all (CatPreds <= p)))
        error (strcat ("fscchi2: 'CategoricalPredictors' must index the", ...
                       " predictors."));
      endif
      iscat(CatPreds) = true;
    else
      error (strcat ("fscchi2: 'CategoricalPredictors' must be indices, a", ...
                     " logical vector or 'all'; naming the predictors", ...
                     " needs a table."));
    endif
  endif

  ## The classes, sorted or in the order ClassNames gives them.  A missing
  ## label leaves its observation out, and so does a class ClassNames omits.
  rowsKept = ! isnan (grp2idx (Y));
  Yk = Y(rowsKept,:);
  if (! isempty (ClassNames))
    [~, errmsg] = namedClasses (uniqueLabels (Yk), ClassNames);
    if (! isempty (errmsg))
      error ("fscchi2: %s", errmsg);
    endif
  endif
  [C, gY] = classOrder (Yk, ClassNames);
  K = classCount (C);
  inC = gY > 0;
  used = find (rowsKept);
  used = used(inC);
  gY = gY(inC);
  X = X(used,:);
  w = Weights(used);

  ## The prior: empirical is taken per predictor, from the observations its
  ## test uses; any other is the same for every predictor
  if (isstruct (Prior))
    Prior = priorFromStruct (Prior, C, 'fscchi2');
  elseif (ischar (Prior) && isrow (Prior) && strcmpi (Prior, 'uniform'))
    Prior = ones (1, K) / K;
  elseif (ischar (Prior) && isrow (Prior) && strcmpi (Prior, 'empirical'))
    Prior = [];
  elseif (isnumeric (Prior) && isvector (Prior))
    if (numel (Prior) != K)
      error ("fscchi2: 'Prior' must have one value per class.");
    endif
    if (! (isreal (Prior) && all (isfinite (Prior)) && all (Prior >= 0)
           && sum (Prior) > 0))
      error (strcat ("fscchi2: 'Prior' must be nonnegative, finite and not", ...
                     " all zero."));
    endif
    Prior = double (Prior(:)') / sum (Prior);
  else
    error (strcat ("fscchi2: 'Prior' must be 'empirical', 'uniform', a", ...
                   " numeric vector or a structure."));
  endif

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

    ## The weights of the observations tested, each class carrying its prior
    ## and the whole summing to their number
    gt = gY(inTest);
    wt = w(inTest);
    if (! isempty (Prior))
      for k = 1:K
        m = gt == k;
        if (any (m) && sum (wt(m)) > 0)
          wt(m) = Prior(k) * wt(m) / sum (wt(m));
        endif
      endfor
    endif
    if (! (sum (wt) > 0))
      continue;
    endif
    wt = wt / sum (wt) * numel (wt);

    ## Pearson's chi-square test on the table of weighted counts, the empty
    ## groups and classes left out
    O = accumarray ([group(inTest), gt], wt);
    O = O(any (O, 2), any (O, 1));
    df = (rows (O) - 1) * (columns (O) - 1);
    if (df < 1)
      continue;
    endif
    E = sum (O, 2) * sum (O, 1) / sum (O(:));
    stat = sum ((O(:) - E(:)) .^ 2 ./ E(:));
    scores(j) = -log (chi2cdf (stat, df, 'upper'));
  endfor

  [~, idx] = sort (scores, 'descend');

endfunction

## Test output: the binning, measured on R2024a
%!shared xb, c2, c3, w
%! xb = ((1:40) .^ 2)';
%! c2 = [repmat({'p'}, 12, 1); repmat({'q'}, 10, 1); ...
%!       repmat({'p'}, 8, 1); repmat({'q'}, 10, 1)];
%! c3 = repmat ({'a'; 'b'; 'c'; 'a'; 'b'}, 8, 1);
%! w = 1 + mod ((1:40)', 3);
%!test  # ten bins, the default
%! [idx, s] = fscchi2 (xb, c2);
%! assert_equal (idx, 1);
%! assert_equal (s, 8.52155983968549, 1e-13);
%!test  # the remainder goes to the first bins
%! [~, s] = fscchi2 (xb, c2, 'NumBins', 3);
%! assert_equal (s, 5.8021978021978, 1e-13);
%!test
%! [~, s] = fscchi2 (xb, c2, 'NumBins', 7);
%! assert_equal (s, 11.6724134015208, 1e-13);
%!test  # ties share a bin
%! xt = [ones(15, 1); 2 * ones(10, 1); (3:17)'];
%! [~, s] = fscchi2 (xt, c2, 'NumBins', 3);
%! assert_equal (s, 4.75128205128205, 1e-13);
%!test  # few distinct values, kept apart by these edges
%! xu = [repmat([1; 2; 3], 13, 1); 1];
%! [~, s] = fscchi2 (xu, c2);
%! assert_equal (s, 0.219780219780219, 1e-13);
%!test  # equal values share a bin, here the two smallest
%! [~, s] = fscchi2 (mod (xb, 7), c3, 'NumBins', 4);
%! assert_equal (s, 0.0761861991737764, 1e-14);
%!test  # and so does naming the predictor categorical
%! xu = [repmat([1; 2; 3], 13, 1); 1];
%! [~, s] = fscchi2 (xu, c2, 'CategoricalPredictors', 1);
%! assert_equal (s, 0.219780219780219, 1e-13);

## Weights and the prior, measured on R2024a
%!test  # weights
%! [~, s] = fscchi2 (xb, c2, 'Weights', w);
%! assert_equal (s, 8.57212519846196, 1e-13);
%!test  # the uniform prior
%! [~, s] = fscchi2 (xb, c3, 'Prior', 'uniform');
%! assert_equal (s, 0.00381024182234503, 1e-15);
%!test  # a given prior
%! [~, s] = fscchi2 (xb, c3, 'Prior', [0.6, 0.3, 0.1]);
%! assert_equal (s, 0.00021661173625569, 1e-15);
%!test  # a prior as a structure
%! P = struct ('ClassNames', {{'c', 'a', 'b'}}, 'ClassProbs', [0.1, 0.6, 0.3]);
%! [~, s] = fscchi2 (xb, c3, 'Prior', P);
%! assert_equal (s, 0.00021661173625569, 1e-15);
%!test  # weights and a prior
%! [~, s] = fscchi2 (xb, c3, 'Weights', w, 'Prior', 'uniform');
%! assert_equal (s, 0.100179060091692, 1e-13);

## Missing values, measured on R2024a
%!test  # a predictor is tested on the observations it is known for
%! xn = xb;
%! xn([3, 17, 29]) = NaN;
%! [~, s] = fscchi2 (xn, c2);
%! assert_equal (s, 6.95245500311045, 1e-13);
%!test  # with weights, the empirical prior taken over those observations
%! xn = xb;
%! xn([3, 17, 29]) = NaN;
%! [~, s] = fscchi2 (xn, c2, 'Weights', w);
%! assert_equal (s, 7.41096325113626, 1e-13);
%!test  # a given prior is carried by those observations
%! xn = xb;
%! xn([3, 17, 29]) = NaN;
%! [~, s] = fscchi2 (xn, c3, 'Prior', [0.6, 0.3, 0.1]);
%! assert_equal (s, 0.000358852938452796, 1e-15);
%!test  # an observation of missing class is left out
%! cm = c3;
%! cm([4, 9]) = {''};
%! [~, s] = fscchi2 (xb, cm);
%! assert_equal (s, 0.00110436213552064, 1e-15);
%!test  # and its weight with it
%! cm = c3;
%! cm([4, 9]) = {''};
%! [~, s] = fscchi2 (xb, cm, 'Weights', w);
%! assert_equal (s, 0.162559576661729, 1e-13);

## Octave tests a predictor's missing values as one more group when asked,
## as MATLAB documents; R2024a and R2026b ignore the option
%!test  # a categorical predictor: chi-square 6 on 2 degrees of freedom
%! x = [1; 1; 2; 2; NaN; NaN];
%! y = {'a'; 'a'; 'b'; 'b'; 'b'; 'b'};
%! [~, s] = fscchi2 (x, y, 'CategoricalPredictors', 1, 'UseMissing', true);
%! assert_equal (s, 3, 1e-14);
%!test  # without it, chi-square 4 on 1 degree of freedom
%! x = [1; 1; 2; 2; NaN; NaN];
%! y = {'a'; 'a'; 'b'; 'b'; 'b'; 'b'};
%! [~, s] = fscchi2 (x, y, 'CategoricalPredictors', 1);
%! assert_equal (s, -log (chi2cdf (4, 1, 'upper')), 1e-14);
%!test  # a continuous predictor's missing values form one more bin
%! x = [1; 2; 3; 4; NaN; NaN];
%! y = {'a'; 'a'; 'b'; 'b'; 'b'; 'b'};
%! [~, s] = fscchi2 (x, y, 'NumBins', 2, 'UseMissing', true);
%! assert_equal (s, 3, 1e-14);

## What cannot separate the classes scores 0, measured on R2024a
%!test  # ranking keeps the order of equal scores
%! [idx, s] = fscchi2 ([xb, ones(40, 1), NaN(40, 1), xb], c2);
%! assert_equal (idx, [1, 4, 2, 3]);
%! assert_equal (s, [8.52155983968549, 0, 0, 8.52155983968549], 1e-13);
%!test
%! [idx, s] = fscchi2 ([xb, xb, mod(xb, 7)], c3, 'NumBins', 4);
%! assert_equal (idx, [3, 1, 2]);
%! assert_equal (s, [0, 0, 0.0761861991737764], 1e-14);
%!test  # a single bin
%! [~, s] = fscchi2 (xb, c3, 'NumBins', 1);
%! assert_equal (s, 0);
%!test  # a single class
%! [~, s] = fscchi2 (xb, repmat ({'a'}, 40, 1));
%! assert_equal (s, 0);
%!test  # more bins than distinct ranks
%! [~, s] = fscchi2 (xb, c3, 'NumBins', 100);
%! assert_equal (s, 0.877012088871555, 1e-13);

## The response, measured on R2024a
%!test  # the type of the labels does not matter
%! [~, s] = fscchi2 (xb, c3);
%! assert_equal (s, 0.001140903415465, 1e-15);
%! assert_equal (nthargout (2, @fscchi2, xb, categorical (c3)), s);
%! assert_equal (nthargout (2, @fscchi2, xb, char (c3)), s);
%! assert_equal (nthargout (2, @fscchi2, xb, double (categorical (c3))), s);
%!test  # logical labels
%! [~, s] = fscchi2 (xb, mod ((1:40)', 2) == 0);
%! assert_equal (s, 0);
%!test  # a subset of the classes
%! [~, s] = fscchi2 (xb, c3, 'ClassNames', {'a', 'b'});
%! assert_equal (s, 0.0242309863346623, 1e-15);
%!test  # categorical predictors given by index and as a logical vector
%! X = [xb, mod((1:40)', 4)];
%! [~, s] = fscchi2 (X, c3, 'CategoricalPredictors', 2);
%! assert_equal (s, [0.001140903415465, 0], 1e-15);
%! [~, s] = fscchi2 (X, c3, 'CategoricalPredictors', [false, true]);
%! assert_equal (s, [0.001140903415465, 0], 1e-15);

## A table, measured on R2024a
%!test  # a categorical variable with an undefined value
%! k = categorical (repmat ({'u'; 'v'; 'w'; 'u'}, 10, 1));
%! k([2, 9, 30]) = missing;
%! [~, s] = fscchi2 (table (k, xb, c3), 'c3');
%! assert_equal (s, [0.0566685023986232, 0.001140903415465], 1e-15);
%!test  # a formula ranks the predictors it names, in the table's order
%! k = categorical (repmat ({'u'; 'v'; 'w'; 'u'}, 10, 1));
%! k([2, 9, 30]) = missing;
%! [idx, s] = fscchi2 (table (xb, k, c3), 'c3 ~ k + xb');
%! assert_equal (idx, [2, 1]);
%! assert_equal (s, [0.001140903415465, 0.0566685023986232], 1e-15);
%!test  # weights named in the table are no predictor
%! x1 = (1:12)';
%! w = [1; 1; 2; 2; 1; 1; 3; 3; 1; 1; 2; 2];
%! lab = repmat ({'p'; 'q'}, 6, 1);
%! [idx, s] = fscchi2 (table (x1, w, lab), 'lab', 'Weights', 'w');
%! assert_equal (idx, 1);
%! assert_equal (s, 0.484012943062914, 1e-13);
%!test  # an ordinal variable is continuous, unless named categorical
%! o12 = categorical ((1:12)', 1:12, 'Ordinal', true);
%! n12 = categorical ((1:12)');
%! v = (1:12)';
%! lab = repmat ({'p'; 'q'; 'q'; 'p'}, 3, 1);
%! T = table (o12, n12, v, lab);
%! [~, s] = fscchi2 (T, 'lab');
%! assert_equal (s, [0.627085663064649, 1.01158205543499, ...
%!                   0.627085663064649], 1e-13);
%! [~, s] = fscchi2 (T, 'lab', 'CategoricalPredictors', 3);
%! assert_equal (s, [0.627085663064649, 1.01158205543499, ...
%!                   1.01158205543499], 1e-13);
%! [~, s] = fscchi2 (T, 'lab', 'CategoricalPredictors', 'all');
%! assert_equal (s, [1.01158205543499, 1.01158205543499, ...
%!                   1.01158205543499], 1e-13);

## Test input validation
%!error<Invalid call to fscchi2> fscchi2 (1)
%!error<fscchi2: optional arguments must be Name-Value pairs.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'NumBins')
%!error<fscchi2: unknown parameter name 'PredictorNames'.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'PredictorNames', {'x'})
%!error<fscchi2: X must be a nonempty real matrix of type double or single.> ...
%! fscchi2 (int8 ([1; 2]), {'a'; 'b'})
%!error<fscchi2: X and Y must have the same number of observations.> ...
%! fscchi2 ([1; 2; 3], {'a'; 'b'})
%!error<fscchi2: 'NumBins' must be a positive integer.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'NumBins', 0)
%!error<fscchi2: 'NumBins' must be a positive integer.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'NumBins', 2.5)
%!error<fscchi2: 'UseMissing' must be a logical scalar.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'UseMissing', 2)
%!error<fscchi2: 'Weights' must have one value per observation.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'Weights', [1; 1; 1])
%!error<fscchi2: 'Weights' must be nonnegative, finite and not all zero.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'Weights', [1; -1])
%!error<fscchi2: 'Weights' must be nonnegative, finite and not all zero.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'Weights', [1; Inf])
%!error<fscchi2: 'CategoricalPredictors' as a logical vector must have one value per predictor.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'CategoricalPredictors', [true, false])
%!error<fscchi2: 'CategoricalPredictors' must index the predictors.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'CategoricalPredictors', 2)
%!error<fscchi2: 'CategoricalPredictors' must be indices, a logical vector or 'all'; naming the predictors needs a table.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'CategoricalPredictors', 'x1')
%!error<fscchi2: 'Prior' must have one value per class.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'Prior', [1, 2, 3])
%!error<fscchi2: 'Prior' must be nonnegative, finite and not all zero.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'Prior', [1, -1])
%!error<fscchi2: 'Prior' must be 'empirical', 'uniform', a numeric vector or a structure.> ...
%! fscchi2 ([1; 2], {'a'; 'b'}, 'Prior', 'other')
