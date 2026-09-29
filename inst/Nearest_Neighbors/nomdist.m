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
## @deftypefn  {statistics} {@var{D} =} nomdist (@var{X})
## @deftypefnx {statistics} {@var{D} =} nomdist (@var{X}, @var{measure})
## @deftypefnx {statistics} {@var{D} =} nomdist (@dots{}, @var{Name}, @var{Value})
##
## Dissimilarity between every pair of rows of nominal data.
##
## @code{@var{D} = nomdist (@var{X})} returns the Goodall 3 dissimilarity
## between each pair of rows of @var{X}, whose columns are nominal variables
## and whose values are their levels.  @var{D} is a row vector in the order
## @code{pdist} returns, pairs @math{(1,2)}, @math{(1,3)}, @dots{},
## @math{(2,3)}, @dots{}, so @code{squareform} turns it into a matrix and
## @code{linkage} takes it as it stands.  An @var{X} of one row gives a
## @math{1*0} @var{D}.
##
## @var{X} is a numeric, logical, categorical, string or cellstr matrix, each
## distinct value being a level, or a table whose every variable is one column
## of these.  A missing value, @code{NaN}, an empty character vector, a
## missing string or an undefined category, is refused.
##
## The measures weigh a match or a mismatch by how often the levels involved
## occur, counted over the rows of @var{X} unless @qcode{'Reference'} names
## another sample.
##
## @code{@var{D} = nomdist (@var{X}, @var{measure})} uses the measure named,
## in any case, from those of Boriah, Chandola and Kumar (2008) and of Sulc and
## Rezankova (2019), as the R package @code{nomclust} 2.8.1 implements them.
## In the table, @math{N} is the number of rows the frequencies are counted
## over, @math{f} the count of a level, @math{p = f/N} its frequency,
## @math{p2 = f(f-1)/(N(N-1))} the chance of drawing it twice, and @math{n}
## the number of levels the variable holds.  A subscript names the level of
## either row, and a variable scores 0 wherever its entry names no score.
##
## @multitable @columnfractions 0.16 0.02 0.82
## @headitem @var{measure} @tab @tab Score of one variable
## @item @qcode{'anderberg'} @tab @tab
## The similarity is
## @tex
## \(\frac{\sum_{x = y} c/p_x^2}
##       {\sum_{x = y} c/p_x^2 + \sum_{x \ne y} c/(2 p_x p_y)},\)
## @end tex
## @ifnottex
## the sum over the matches of @math{c/p^2}, divided by itself plus the sum over
## the mismatches of @math{c/(2 p_x p_y)},
## @end ifnottex
## with
## @tex
## \(c = 2/(n(n+1)).\)
## @end tex
## @ifnottex
## @math{c = 2/(n(n+1))}.
## @end ifnottex
## @item @qcode{'burnaby'} @tab @tab
## 1 on a match; on a mismatch
## @tex
## \(\frac{B}{\log \frac{p_x p_y}{(1-p_x)(1-p_y)} + B},\)
## @end tex
## @ifnottex
## @math{B/(log (p_x p_y/((1-p_x)(1-p_y))) + B)},
## @end ifnottex
## with
## @tex
## \(B = \sum_q 2 \log (1 - p_q).\)
## @end tex
## @ifnottex
## @math{B} the sum of @math{2 log (1-p)} over the levels.
## @end ifnottex
## @item @qcode{'eskin'} @tab @tab
## 1 on a match,
## @tex
## \(\frac{n^2}{n^2 + 2}\)
## @end tex
## @ifnottex
## @math{n^2/(n^2+2)}
## @end ifnottex
## on a mismatch.
## @item @qcode{'gambaryan'} @tab @tab
## @tex
## \(-\left(p_x \log_2 p_x + (1-p_x) \log_2 (1-p_x)\right)\)
## @end tex
## @ifnottex
## @math{-(p log2 (p) + (1-p) log2 (1-p))}
## @end ifnottex
## on a match, which peaks where @math{p} is 0.5.
## @item @qcode{'goodall1'} @tab @tab
## On a match
## @tex
## \(1 - \sum_{q:\, p_q \le p_x} \mathit{p2}_q.\)
## @end tex
## @ifnottex
## 1 less the sum of @math{p2} over the levels at most as frequent as the
## one matched.
## @end ifnottex
## @item @qcode{'goodall2'} @tab @tab
## On a match
## @tex
## \(1 - \sum_{q:\, p_q \ge p_x} \mathit{p2}_q.\)
## @end tex
## @ifnottex
## 1 less the sum of @math{p2} over the levels at least as frequent as the
## one matched.
## @end ifnottex
## @item @qcode{'goodall3'} @tab @tab
## @tex
## \(1 - \mathit{p2}_x\)
## @end tex
## @ifnottex
## @math{1 - p2}
## @end ifnottex
## on a match.
## @item @qcode{'goodall4'} @tab @tab
## @tex
## \(\mathit{p2}_x\)
## @end tex
## @ifnottex
## @math{p2}
## @end ifnottex
## on a match.
## @item @qcode{'iof'} @tab @tab
## 1 on a match,
## @tex
## \(\frac{1}{1 + \log f_x \log f_y}\)
## @end tex
## @ifnottex
## @math{1/(1 + log (f_x) log (f_y))}
## @end ifnottex
## on a mismatch.
## @item @qcode{'lin'} @tab @tab
## @tex
## \(2 \log p_x\)
## @end tex
## @ifnottex
## @math{2 log (p)}
## @end ifnottex
## on a match,
## @tex
## \(2 \log (p_x + p_y)\)
## @end tex
## @ifnottex
## @math{2 log (p_x + p_y)}
## @end ifnottex
## on a mismatch, over the sum of
## @tex
## \(\log p_x + \log p_y.\)
## @end tex
## @ifnottex
## @math{log (p_x) + log (p_y)}.
## @end ifnottex
## @item @qcode{'lin1'} @tab @tab
## On a match
## @tex
## \(\log \sum_{q:\, p_q = p_x} p_q;\)
## @end tex
## @ifnottex
## the log of the summed @math{p} of the levels as frequent as the one matched;
## @end ifnottex
## on a mismatch
## @tex
## \(2 \log \sum_{q:\, \min (p_x, p_y) \le p_q \le \max (p_x, p_y)} p_q;\)
## @end tex
## @ifnottex
## twice the log of the summed @math{p} of the levels whose frequency lies from
## @math{p_x} to @math{p_y};
## @end ifnottex
## over the same sum as @qcode{'lin'}.
## @item @qcode{'of'} @tab @tab
## 1 on a match,
## @tex
## \(\frac{1}{1 + \log (N/f_x) \log (N/f_y)}\)
## @end tex
## @ifnottex
## @math{1/(1 + log (N/f_x) log (N/f_y))}
## @end ifnottex
## on a mismatch.
## @item @qcode{'sm'} @tab @tab
## 1 on a match, the simple matching coefficient.
## @item @qcode{'smirnov'} @tab @tab
## On a match
## @tex
## \(2 + \frac{N - f_x}{f_x} + \sum_{q \ne x} \frac{f_q}{N - f_q};\)
## @end tex
## @ifnottex
## @math{2 + (N-f)/f} plus the sum of @math{f/(N-f)} over the other levels;
## @end ifnottex
## on a mismatch
## @tex
## \(\sum_{q \ne x, y} \frac{f_q}{N - f_q}.\)
## @end tex
## @ifnottex
## that sum over the levels other than both.
## @end ifnottex
## @item @qcode{'ve'} @tab @tab
## On a match the entropy of the variable,
## @tex
## \(-\frac{1}{\log n} \sum_q p_q \log p_q.\)
## @end tex
## @ifnottex
## @math{-sum (p log (p))/log (n)}.
## @end ifnottex
## @item @qcode{'vm'} @tab @tab
## On a match the Gini index of the variable,
## @tex
## \(\frac{n}{n-1} \left(1 - \sum_q p_q^2\right).\)
## @end tex
## @ifnottex
## @math{(1 - sum (p^2)) n/(n-1)}.
## @end ifnottex
## @end multitable
##
## The scores of a pair are summed over the variables, each times its weight,
## and divided by the sum of the weights, except where the table names another
## divisor, and @qcode{'gambaryan'} and @qcode{'smirnov'} divide by the number
## of levels over every variable.  With @math{S} that similarity,
## @qcode{'eskin'}, @qcode{'iof'}, @qcode{'lin'}, @qcode{'lin1'} and
## @qcode{'of'} return @math{1/S - 1}, @qcode{'smirnov'} returns
## @math{1/(1+S)}, and every other measure @math{1 - S}@.  Several measures,
## the four Goodall ones among them, place two identical rows at a nonzero
## dissimilarity, a match on a common level saying less than one on a rare
## level.  @qcode{'lin'} and @qcode{'lin1'} return @code{Inf} for a pair of
## similarity 0, such as two rows differing in every variable where each
## variable holds two levels.
##
## @code{@var{D} = nomdist (@dots{}, @var{Name}, @var{Value})} takes the
## following options.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
## @item @qcode{'Weights'} @tab @tab One weight per variable, each from 0 to
## 1 and at least one positive; a variable of weight 0 is left out.  Not taken
## by @qcode{'anderberg'}, @qcode{'gambaryan'} and @qcode{'smirnov'}.
## @item @qcode{'Reference'} @tab @tab A sample to count level frequencies
## over in place of @var{X}, of the same kind: a matrix with as many columns,
## the same type in each, or a table holding every variable of @var{X}, read
## by name.  It must hold every level @var{X} holds.  Given the whole data,
## it measures a subsample by how rare its levels are in the whole.
## @end multitable
##
## Where @code{nomclust} answers otherwise, the answer here is the limit of
## the measure's formula.  @qcode{'lin'} and @qcode{'lin1'} sum the counts
## of levels, so a sum over every level is exactly 1 and a similarity of 0 is
## exactly 0, where @code{nomclust} rounds it to about @math{1e16} or
## replaces it with one more than the largest finite dissimilarity.
## @qcode{'gambaryan'} scores a level every row holds as 0, where
## @code{nomclust} returns @code{NaN}.  Where every variable holds one level,
## @qcode{'lin'} returns 0 and @qcode{'lin1'} 1, where @code{nomclust}
## returns @code{NaN}.
##
## References:
## @itemize
## @item Boriah, S., Chandola, V. and Kumar, V. (2008).  Similarity measures
## for categorical data: a comparative evaluation.  Proceedings of the 8th
## SIAM International Conference on Data Mining, 243-254.
## @item Sulc, Z. and Rezankova, H. (2019).  Comparison of similarity
## measures for categorical data in hierarchical clustering.  Journal of
## Classification, 36(1), 58-72.
## @end itemize
##
## @seealso{nomdist2, pdist, squareform, linkage}
## @end deftypefn

function D = nomdist (X, varargin)

  ## Input validation
  if (nargin < 1)
    error ("nomdist: too few input arguments.");
  endif
  [C, measure, w, errmsg] = __nomprep__ (X, [], false, varargin);
  if (! isempty (errmsg))
    error ("nomdist: %s", errmsg);
  endif

  D = __nomdist__ (C{1}, [], C{3}, measure, w);

endfunction

%!shared X
%! X = [1, 1, 1; 1, 2, 1; 1, 1, 2; 2, 2, 2; 2, 1, 1; 3, 2, 2; 4, 1, 3];
%!test
%! ## The default is Goodall 3; values from nomclust 2.8.1
%! D = [0.4285714285714286, 0.4761904761904762, 1, 0.4761904761904762, 1, ...
%!      0.7619047619047619, 0.7142857142857143, 0.7142857142857143, ...
%!      0.7142857142857143, 0.7142857142857143, 1, 0.7142857142857143, ...
%!      0.7619047619047619, 0.7142857142857143, 0.7619047619047619, ...
%!      0.6825396825396826, 0.4285714285714286, 1, 1, 0.7619047619047619, ...
%!      1];
%! assert_equal (nomdist (X), D, -1e-14);
%!test
%! D = [0.3333333333333334, 0.3333333333333334, 1, 0.3333333333333334, 1, ...
%!      0.6666666666666667, 0.6666666666666667, 0.6666666666666667, ...
%!      0.6666666666666667, 0.6666666666666667, 1, 0.6666666666666667, ...
%!      0.6666666666666667, 0.6666666666666667, 0.6666666666666667, ...
%!      0.6666666666666667, 0.3333333333333334, 1, 1, 0.6666666666666667, ...
%!      1];
%! assert_equal (nomdist (X, 'SM'), D, -1e-14);
%!assert_equal (nomdist (X, "lin"), nomdist (X, 'lin'))
%!test
%! ## Weighted Lin; values from nomclust 2.8.1
%! D = nomdist (X, 'lin', 'Weights', [0.7, 1, 0.4]);
%! assert_equal (D(1:3), [0.7547596113257813, 0.2283122503928992, ...
%!                        4.98065967640019], -1e-14);
%!assert_equal (nomdist (X, 'Weights', [1, 1, 1]), nomdist (X))
%!assert_equal (nomdist (X, 'Weights', []), nomdist (X))
%!assert_equal (nomdist (X, 'lin', 'Weights', [1, 0, 1]), ...
%!              nomdist (X(:,[1, 3]), 'lin'))
%!assert_equal (nomdist (10 * X - 25), nomdist (X))
%!assert_equal (nomdist (X == 1), nomdist (double (X == 1)))
%!assert_equal (nomdist (int8 (X)), nomdist (X))
%!test
%! C = {'a', 'b', 'c', 'd'};
%! assert_equal (nomdist (C(X)), nomdist (X));
%!test
%! C = {'a', 'b', 'c', 'd'};
%! assert_equal (nomdist (string (C(X))), nomdist (X));
%!test
%! C = {'a', 'b', 'c', 'd'};
%! assert_equal (nomdist (categorical (C(X))), nomdist (X));
%!test
%! C = {'a', 'b', 'c', 'd'};
%! T = table (categorical (C(X(:,1))'), X(:,2) == 1, C(X(:,3))', ...
%!            'VariableNames', {'A', 'B', 'C'});
%! assert_equal (nomdist (T, 'of'), nomdist (X, 'of'));
%!assert_equal (nomdist ([1, 2, 3]), zeros (1, 0))
%!assert_equal (nomdist ([1; 1; 2], 'goodall3'), [1/3, 1, 1], -1e-14)
%!assert_equal (nomdist ([1; 1; 2], 'goodall3', 'Reference', ...
%!                      [1; 1; 1; 1; 2]), [0.6, 1, 1], -1e-14)
%!assert_equal (nomdist ([1; 1; 2], 'Reference', []), nomdist ([1; 1; 2]))
%!test
%! T = table ([1; 1; 2], {'a'; 'b'; 'a'}, 'VariableNames', {'A', 'B'});
%! R = table ({'a'; 'b'; 'a'; 'a'}, [1; 1; 2; 1], 'VariableNames', {'B', 'A'});
%! assert_equal (nomdist (T, 'Reference', R), ...
%!               nomdist ([1, 1; 1, 2; 2, 1], 'Reference', ...
%!                        [1, 1; 1, 2; 2, 1; 1, 1]));

%!error<nomdist: too few input arguments.> nomdist ()
%!error<nomdist: MEASURE must be one of 'anderberg', 'burnaby', 'eskin', 'gambaryan', 'goodall1', 'goodall2', 'goodall3', 'goodall4', 'iof', 'lin', 'lin1', 'of', 'sm', 'smirnov', 've', 'vm'.> ...
%! nomdist ([1; 2], 'hamming')
%!error<nomdist: MEASURE must be one of 'anderberg', 'burnaby', 'eskin', 'gambaryan', 'goodall1', 'goodall2', 'goodall3', 'goodall4', 'iof', 'lin', 'lin1', 'of', 'sm', 'smirnov', 've', 'vm'.> ...
%! nomdist ([1; 2], 1)
%!error<nomdist: invalid optional paired argument.> ...
%! nomdist ([1; 2], 'Scale', 1)
%!error<nomdist: X must be a numeric, logical, categorical, string or cellstr matrix, or a table.> ...
%! nomdist (['ab'; 'cd'])
%!error<nomdist: X must be a numeric, logical, categorical, string or cellstr matrix, or a table.> ...
%! nomdist (ones (2, 2, 2))
%!error<nomdist: X must not be empty.> nomdist (zeros (0, 2))
%!error<nomdist: every variable of X must be one column of levels.> ...
%! nomdist (table ([1, 2; 3, 4]))
%!error<nomdist: X must not hold missing values.> nomdist ([1; NaN])
%!error<nomdist: X must not hold missing values.> nomdist ({'a'; ''})
%!error<nomdist: 'Reference' must be a table, as X is.> ...
%! nomdist (table ([1; 2]), 'Reference', [1; 2])
%!error<nomdist: 'Reference' must hold every variable of X.> ...
%! nomdist (table ([1; 2]), 'Reference', table ([1; 2], 'VariableNames', {'Z'}))
%!error<nomdist: 'Reference' must be a matrix, as X is.> ...
%! nomdist ([1; 2], 'Reference', table ([1; 2]))
%!error<nomdist: 'Reference' must be a numeric, logical, categorical, string or cellstr matrix, or a table.> ...
%! nomdist ([1; 2], 'Reference', ['a'; 'b'])
%!error<nomdist: 'Reference' must have as many columns as X.> ...
%! nomdist ([1; 2], 'Reference', [1, 1; 2, 2])
%!error<nomdist: 'Reference' must hold the same type as X in every variable.> ...
%! nomdist ([1; 2], 'Reference', {'a'; 'b'})
%!error<nomdist: 'Reference' must not hold missing values.> ...
%! nomdist ([1; 2], 'Reference', [1; 2; NaN])
%!error<nomdist: 'Reference' must hold every level that X holds.> ...
%! nomdist ([1; 2], 'Reference', [1; 1])
%!error<nomdist: 'Weights' does not apply to the 'smirnov' measure.> ...
%! nomdist ([1; 2], 'smirnov', 'Weights', 1)
%!error<nomdist: 'Weights' must be a vector of 2 values between 0 and 1.> ...
%! nomdist ([1, 1; 2, 1], 'Weights', 1)
%!error<nomdist: 'Weights' must be a vector of 2 values between 0 and 1.> ...
%! nomdist ([1, 1; 2, 1], 'Weights', [1, 2])
%!error<nomdist: 'Weights' must hold at least one positive value.> ...
%! nomdist ([1, 1; 2, 1], 'Weights', [0, 0])
