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
## @deftypefn  {statistics} {@var{kappa} =} cohenKappa (@var{x}, @var{y})
## @deftypefnx {statistics} {@var{kappa} =} cohenKappa (@var{T})
## @deftypefnx {statistics} {@var{kappa} =} cohenKappa (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{kappa}, @var{ci}] =} cohenKappa (@dots{})
## @deftypefnx {statistics} {[@var{kappa}, @var{ci}, @var{stats}] =} cohenKappa (@dots{})
##
## Cohen's kappa, the agreement between two raters beyond chance.
##
## @code{@var{kappa} = cohenKappa (@var{x}, @var{y})} returns Cohen's kappa
## for two raters who each assigned the same subjects to categories, the
## ratings of the first in the vector @var{x} and of the second in @var{y}.
## The ratings may be numeric or logical, categorical, or text: a character
## matrix with one rating per row, a cell array of character vectors, or a
## string array, the two raters' of the same kind.  A subject whose rating by
## either rater is missing, @qcode{NaN}, an undefined category, an empty text
## or a missing string, is left out.
##
## @code{@var{kappa} = cohenKappa (@var{T})} takes the agreement table
## instead, a square matrix of counts whose element @code{@var{T}(i,j)} is
## the number of subjects the first rater assigned to category @math{i} and
## the second to category @math{j}.
##
## With @code{po} the proportion of subjects on whom the raters agree and
## @code{pe} the proportion expected by chance from each rater's own
## frequencies of the categories,
## @tex
## $$ \kappa = \frac{p_o - p_e}{1 - p_e}, $$
##
## @end tex
## @ifnottex
## @code{kappa = (po - pe) / (1 - pe)},
## @end ifnottex
## 1 for perfect agreement, 0 for agreement no better than chance, and
## negative for less.  Weighted kappa counts partial agreement between
## ordered categories, weighting each pair of categories @math{i} and
## @math{j} by @code{w(i,j)}, 1 when they are the same, in both @code{po}
## and @code{pe}.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.2 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Weights'} @tab @qcode{'none'} (default), unweighted kappa;
## @qcode{'linear'}, @code{w(i,j) = 1 - abs (i - j) / (c - 1)};
## @qcode{'quadratic'}, @code{w(i,j) = 1 - (i - j)^2 / (c - 1)^2}, for
## @math{c} categories in their order; or a @math{c}-by-@math{c} matrix of
## agreement weights between 0 and 1, with ones on its diagonal.
##
## @item @qcode{'Categories'} @tab The categories, in their order: numeric
## for numeric ratings, text for text or categorical ratings.  The default is
## the sorted values of numeric ratings, the categories of categorical
## ratings in their own order, unused ones included, and the sorted texts of
## text ratings.  A category neither rater used counts towards @math{c}
## and spaces the others: placed between used categories it changes linear
## and quadratic kappa, placed after them it does not.  Not accepted with
## @var{T}, whose rows and columns are the categories.
##
## @item @qcode{'Alpha'} @tab The significance level of the confidence
## interval, a scalar in the range @math{(0,1)}.  The default is
## @math{0.05}.
## @end multitable
##
## @code{[@var{kappa}, @var{ci}] = cohenKappa (@dots{})} also returns the
## @math{100 (1 - alpha)} percent confidence interval of kappa, from its
## large-sample standard error (Fleiss, Cohen and Everitt, 1969), clipped to
## @math{[-1, 1]}.
##
## @code{[@var{kappa}, @var{ci}, @var{stats}] = cohenKappa (@dots{})} also
## returns a structure with the fields
##
## @multitable @columnfractions 0.2 0.8
## @item @qcode{po} @tab the observed agreement, weighted if weights are
## given
## @item @qcode{pe} @tab the agreement expected by chance
## @item @qcode{se} @tab the standard error of kappa that @var{ci} is based
## on
## @item @qcode{zval} @tab the statistic of the test that kappa is 0, kappa
## over its standard error when there is no agreement beyond chance
## @item @qcode{p} @tab the two-sided p-value of that test
## @item @qcode{n} @tab the number of subjects
## @item @qcode{table} @tab the agreement table
## @item @qcode{categories} @tab the categories, a numeric column or a cell
## array of character vectors
## @item @qcode{weights} @tab the matrix of weights used
## @end multitable
##
## Kappa is undefined when chance alone predicts complete agreement, as when
## both raters use a single category, and such ratings are refused.
## Agreement among more than two raters is measured by @code{fleissKappa}.
## MATLAB has no Cohen's kappa, so @code{cohenKappa} is specific to Octave.
##
## References:
## @enumerate
## @item
## J. Cohen.  A coefficient of agreement for nominal scales.
## @emph{Educational and Psychological Measurement}, 20(1):37--46, 1960.
## @item
## J. Cohen.  Weighted kappa: nominal scale agreement with provision for
## scaled disagreement or partial credit.  @emph{Psychological Bulletin},
## 70(4):213--220, 1968.
## @item
## J. L. Fleiss, J. Cohen and B. S. Everitt.  Large sample standard errors of
## kappa and weighted kappa.  @emph{Psychological Bulletin}, 72(5):323--327,
## 1969.
## @end enumerate
##
## @seealso{fleissKappa, crosstab, limitsOfAgreement}
## @end deftypefn

function [kappa, ci, stats] = cohenKappa (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## An agreement table, or the two raters' ratings
  T = varargin{1};
  istab = (isnumeric (T) && ! isvector (T)
           && (nargin == 1 || ischar (varargin{2}) || isstring (varargin{2})));
  if (istab)
    args = varargin(2:end);
  else
    if (nargin < 2)
      print_usage ();
    endif
    args = varargin(3:end);
  endif

  ## Parse optional Name-Value paired arguments
  optNames = {'Weights', 'Categories', 'Alpha'};
  dfValues = {'none', [], 0.05};
  [wt, cats, alpha, rest] = parsePairedArguments (optNames, dfValues, args(:));
  if (! isempty (rest))
    if (mod (numel (rest), 2) == 0 && ischar (rest{1}))
      error ("cohenKappa: unknown parameter name '%s'.", rest{1});
    endif
    error ("cohenKappa: optional arguments must be Name-Value pairs.");
  endif
  if (! (isfloat (alpha) && isreal (alpha) && isscalar (alpha)
         && alpha > 0 && alpha < 1))
    error ("cohenKappa: 'Alpha' must be a scalar between 0 and 1.");
  endif

  if (istab)
    if (! isempty (cats))
      error ("cohenKappa: 'Categories' cannot be given with a table.");
    endif
    if (! (isreal (T) && rows (T) == columns (T) && ndims (T) == 2
           && all (isfinite (T(:))) && all (T(:) >= 0)))
      error (strcat ("cohenKappa: T must be a square matrix of", ...
                     " non-negative counts."));
    endif
    T = double (T);
    c = rows (T);
    cats = (1:c)';
  else
    [codes, cats, errmsg] = ratingcodes (varargin(1:2), cats);
    if (! isempty (errmsg))
      error ("cohenKappa: %s", errmsg);
    endif
    codes = codes(! any (isnan (codes), 2),:);
    c = numel (cats);
    T = zeros (c, c);
    if (! isempty (codes))
      T = accumarray (codes, 1, [c, c]);
    endif
  endif
  n = sum (T(:));
  if (n == 0)
    error ("cohenKappa: no subject was rated by both raters.");
  endif

  ## The agreement weights
  if (ischar (wt) || (isstring (wt) && isscalar (wt)))
    switch (lower (char (wt)))
      case 'none'
        W = eye (c);
      case 'linear'
        W = 1 - abs ((1:c)' - (1:c)) / max (c - 1, 1);
      case 'quadratic'
        W = 1 - ((1:c)' - (1:c)) .^ 2 / max (c - 1, 1) ^ 2;
      otherwise
        error (strcat ("cohenKappa: 'Weights' must be 'none', 'linear',", ...
                       " 'quadratic' or a matrix of weights."));
    endswitch
  elseif (isnumeric (wt) && isreal (wt) && isequal (size (wt), [c, c])
          && all (wt(:) >= 0 & wt(:) <= 1) && all (diag (wt) == 1))
    W = double (wt);
  else
    error (strcat ("cohenKappa: a matrix of 'Weights' must be square, one", ...
                   " row per category, between 0 and 1, with ones on its", ...
                   " diagonal."));
  endif

  ## Agreement observed and expected by chance
  P = T / n;
  pr = sum (P, 2);
  pc = sum (P, 1);
  po = sum (W(:) .* P(:));
  E = pr * pc;
  pe = sum (W(:) .* E(:));
  if (pe >= 1)
    error (strcat ("cohenKappa: kappa is undefined when chance predicts", ...
                   " complete agreement, as when both raters use a single", ...
                   " category."));
  endif
  kappa = (po - pe) / (1 - pe);

  ## Large-sample variances of Fleiss, Cohen and Everitt, with the weights'
  ## means over each row and column taken against the other rater's margin
  wr = W * pc';
  wc = pr' * W;
  D = wr + wc;
  v = (sum (sum (P .* (W - D * (1 - kappa)) .^ 2)) ...
       - (kappa - pe * (1 - kappa)) ^ 2) / (n * (1 - pe) ^ 2);
  v0 = (sum (sum (E .* (W - D) .^ 2)) - pe ^ 2) / (n * (1 - pe) ^ 2);
  se = sqrt (max (v, 0));
  se0 = sqrt (max (v0, 0));

  ci = min (max (kappa + [-1, 1] * norminv (1 - alpha / 2) * se, -1), 1);
  zval = kappa / se0;
  p = 2 * normcdf (-abs (zval));

  stats = struct ('po', po, 'pe', pe, 'se', se, 'zval', zval, 'p', p, ...
                  'n', n, 'table', T, 'categories', {cats}, 'weights', W);

endfunction

%!demo
%! ## Two raters who sorted 20 subjects into three ordered grades
%! x = [1, 2, 3, 1, 2, 3, 2, 2, 1, 3, 3, 1, 2, 1, 3, 2, 2, 3, 1, 1];
%! y = [1, 2, 2, 1, 3, 3, 2, 1, 1, 3, 2, 1, 2, 2, 3, 2, 3, 3, 1, 2];
%! [kappa, ci] = cohenKappa (x, y)
%! ## Partial agreement between neighbouring grades counted
%! [kappa, ci] = cohenKappa (x, y, 'Weights', 'quadratic')

## Test output
## Values checked against R 4.5.0, irr 0.85 (kappa2) and psych (cohen.kappa)
%!shared x, y
%! x = [1, 2, 3, 1, 2, 3, 2, 2, 1, 3, 3, 1, 2, 1, 3, 2, 2, 3, 1, 1];
%! y = [1, 2, 2, 1, 3, 3, 2, 1, 1, 3, 2, 1, 2, 2, 3, 2, 3, 3, 1, 2];
%!test
%! [k, ci, s] = cohenKappa (x, y);
%! assert_equal (k, 0.473684210526316, 1e-14);
%! assert_equal (ci, [0.156627045921977, 0.790741375130655], 1e-14);
%! assert_equal (s.zval, 3.00357781880274, 1e-13);
%! assert_equal (s.p, 0.00266825301537432, 1e-14);
%! assert_equal (s.n, 20);
%! assert_equal (s.categories, [1; 2; 3]);
%! assert_equal (s.table, [5, 2, 0; 1, 4, 2; 0, 2, 4]);
%!test
%! [k, ci, s] = cohenKappa (x, y, 'Weights', 'linear');
%! assert_equal (k, 0.593023255813953, 1e-14);
%! assert_equal (ci, [0.330633074635, 0.855413436992907], 1e-14);
%! assert_equal (s.zval, 3.42637031579504, 1e-13);
%! assert_equal (s.p, 0.000611705786193184, 1e-14);
%!test
%! [k, ci, s] = cohenKappa (x, y, 'Weights', 'quadratic');
%! assert_equal (k, 0.72, 1e-14);
%! assert_equal (ci, [0.521837223191117, 0.918162776808883], 1e-14);
%! assert_equal (s.zval, 3.22873201500379, 1e-13);
%! assert_equal (s.p, 0.00124340351789742, 1e-14);
%!test
%! W = 1 - abs ((1:3)' - (1:3)) / 2;
%! assert_equal (cohenKappa (x, y, 'Weights', W), ...
%!               cohenKappa (x, y, 'Weights', 'linear'), 1e-15);
%!test
%! ## A category after the used ones leaves weighted kappa unchanged
%! assert_equal (cohenKappa (x, y, 'Weights', 'linear', 'Categories', 1:4), ...
%!               0.593023255813954, 1e-14);
%!test
%! ## One between them does not
%! x4 = x;
%! x4(x == 3) = 4;
%! y4 = y;
%! y4(y == 3) = 4;
%! [k, ci] = cohenKappa (x4, y4, 'Weights', 'linear', 'Categories', 1:4);
%! assert_equal (k, 0.5703125, 1e-14);
%! assert_equal (ci, [0.28133843832902, 0.85928656167098], 1e-14);
%!test
%! [~, ~, s] = cohenKappa (x, y);
%! [k, ci, s2] = cohenKappa (s.table);
%! assert_equal (k, 0.473684210526316, 1e-14);
%! assert_equal (ci, [0.156627045921977, 0.790741375130655], 1e-14);
%! assert_equal (s2.zval, s.zval, 1e-14);
%!test
%! a = {'cat', 'dog', 'dog', 'bird', '', 'cat', 'bird', 'dog', 'cat', 'dog'};
%! b = {'cat', 'dog', 'cat', 'bird', 'dog', 'cat', 'dog', 'dog', 'bird', ''};
%! [k, ci, s] = cohenKappa (a, b);
%! assert_equal (k, 0.428571428571429, 1e-14);
%! assert_equal (ci, [-0.0756238505834626, 0.93276670772632], 1e-14);
%! assert_equal (s.zval, 1.69705627484771, 1e-13);
%! assert_equal (s.p, 0.0896860217703646, 1e-14);
%! assert_equal (s.n, 8);
%! assert_equal (s.categories, {'bird'; 'cat'; 'dog'});
%!test
%! a = {'cat', 'dog', 'dog', 'bird', '', 'cat', 'bird', 'dog', 'cat', 'dog'};
%! b = {'cat', 'dog', 'cat', 'bird', 'dog', 'cat', 'dog', 'dog', 'bird', ''};
%! k = cohenKappa (char (a), char (b));
%! assert_equal (k, 0.428571428571429, 1e-14);
%!test
%! a = {'cat', 'dog', 'dog', 'bird', '', 'cat', 'bird', 'dog', 'cat', 'dog'};
%! b = {'cat', 'dog', 'cat', 'bird', 'dog', 'cat', 'dog', 'dog', 'bird', ''};
%! k = cohenKappa (categorical (a), categorical (b));
%! assert_equal (k, 0.428571428571429, 1e-14);
%!test
%! k = cohenKappa (logical (x == 1), logical (y == 1));
%! assert_equal (k, cohenKappa (double (x == 1), double (y == 1)));
%!test
%! [k, ci, s] = cohenKappa ([1, 2, 1, 2], [1, 2, 1, 2]);
%! assert_equal (k, 1);
%! assert_equal (ci, [1, 1]);
%! assert_equal (s.zval, 2);
%!test
%! [~, ci] = cohenKappa (x, y, 'Alpha', 0.1);
%! assert_equal (ci, 0.473684210526316 + [-1, 1] * norminv (0.95) ...
%!                   * 0.161766831995509, 1e-12);

## Test input validation
%!error <Invalid call to cohenKappa> cohenKappa ()
%!error <Invalid call to cohenKappa> cohenKappa ([1, 2, 3])
%!error <cohenKappa: unknown parameter name 'Foo'.> cohenKappa (x, y, 'Foo', 1)
%!error <cohenKappa: optional arguments must be Name-Value pairs.> ...
%! cohenKappa (x, y, 'Alpha')
%!error <cohenKappa: 'Alpha' must be a scalar between 0 and 1.> ...
%! cohenKappa (x, y, 'Alpha', 1)
%!error <cohenKappa: 'Weights' must be 'none', 'linear', 'quadratic' or a matrix of weights.> ...
%! cohenKappa (x, y, 'Weights', 'cubic')
%!error <cohenKappa: a matrix of 'Weights' must be square, one row per category, between 0 and 1, with ones on its diagonal.> ...
%! cohenKappa (x, y, 'Weights', eye (2))
%!error <cohenKappa: a matrix of 'Weights' must be square, one row per category, between 0 and 1, with ones on its diagonal.> ...
%! cohenKappa (x, y, 'Weights', zeros (3))
%!error <cohenKappa: 'Categories' cannot be given with a table.> ...
%! cohenKappa (eye (2), 'Categories', [1, 2])
%!error <cohenKappa: T must be a square matrix of non-negative counts.> ...
%! cohenKappa ([1, 2, 3; 4, 5, 6])
%!error <cohenKappa: T must be a square matrix of non-negative counts.> ...
%! cohenKappa ([1, -2; 3, 4])
%!error <cohenKappa: ratings of every rater must be of the same kind.> ...
%! cohenKappa ([1, 2], {'a', 'b'})
%!error <cohenKappa: every rater must rate the same number of subjects.> ...
%! cohenKappa ([1, 2, 3], [1, 2])
%!error <cohenKappa: ratings must be numeric, logical, categorical or text vectors.> ...
%! cohenKappa ([1, 2] + 1i, [1, 2])
%!error <cohenKappa: a rating is not one of 'Categories'.> ...
%! cohenKappa ([1, 2, 3], [1, 2, 2], 'Categories', [1, 2])
%!error <cohenKappa: 'Categories' must not repeat a category.> ...
%! cohenKappa ([1, 2, 3], [1, 2, 2], 'Categories', [1, 2, 2, 3])
%!error <cohenKappa: 'Categories' must be a numeric vector for numeric ratings.> ...
%! cohenKappa ([1, 2, 3], [1, 2, 2], 'Categories', {'a'})
%!error <cohenKappa: 'Categories' must be text for text or categorical ratings.> ...
%! cohenKappa ({'a', 'b'}, {'a', 'a'}, 'Categories', [1, 2])
%!error <cohenKappa: kappa is undefined when chance predicts complete agreement, as when both raters use a single category.> ...
%! cohenKappa ([1, 1, 1], [1, 1, 1])
%!error <cohenKappa: no subject was rated by both raters.> ...
%! cohenKappa ([1, NaN], [NaN, 2])
