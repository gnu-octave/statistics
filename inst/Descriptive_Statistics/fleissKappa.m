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
## @deftypefn  {statistics} {@var{kappa} =} fleissKappa (@var{X})
## @deftypefnx {statistics} {@var{kappa} =} fleissKappa (@var{X}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{kappa}, @var{stats}] =} fleissKappa (@dots{})
##
## Fleiss' kappa, the agreement among several raters beyond chance.
##
## @code{@var{kappa} = fleissKappa (@var{X})} returns Fleiss' kappa for
## @math{m} raters who each assigned the same @math{n} subjects to
## categories, @var{X} holding one row per subject and one column per rater.
## The ratings may be numeric or logical, categorical, or text in a cell array
## of character vectors or a string array.  A subject with a missing rating,
## @qcode{NaN}, an undefined category, an empty text or a missing string, is
## left out, every subject having to be rated by all @math{m} raters.
##
## With @code{N(i,j)} the number of raters who put subject @math{i} in
## category @math{j}, @code{p(j)} the proportion of all ratings in category
## @math{j}, @code{Pbar} the mean over the subjects of
## @code{(sum (N(i,:) .^ 2) - m) / (m * (m - 1))}, the proportion of pairs
## of raters who agree, and @code{Pe = sum (p .^ 2)},
## @tex
## $$ \kappa = \frac{\bar P - P_e}{1 - P_e}, $$
##
## @end tex
## @ifnottex
## @code{kappa = (Pbar - Pe) / (1 - Pe)},
## @end ifnottex
## 1 for perfect agreement and 0 for agreement no better than chance.  With
## two raters it is Scott's pi, not Cohen's kappa, chance being estimated
## from the ratings of all raters together; @code{cohenKappa} measures the
## agreement of two raters from their own frequencies.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.2 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Categories'} @tab The categories: numeric for numeric
## ratings, text for text or categorical ratings.  The default is the sorted
## values of numeric ratings, the categories of categorical ratings, unused
## ones included, and the sorted texts of text ratings.  A category no rater
## used has a kappa of @qcode{NaN} in @var{stats}.
##
## @item @qcode{'Exact'} @tab @qcode{false} (default) for Fleiss' kappa, or
## @qcode{true} for the variant of Conger (1980), which takes each rater's own
## frequencies into account and equals Cohen's kappa for two raters.  It has
## no test, and @var{stats} then holds @qcode{NaN} for its statistic and
## p-value.
## @end multitable
##
## @code{[@var{kappa}, @var{stats}] = fleissKappa (@dots{})} also returns a
## structure with the fields
##
## @multitable @columnfractions 0.2 0.8
## @item @qcode{po} @tab the observed agreement @code{Pbar}
## @item @qcode{pe} @tab the agreement expected by chance
## @item @qcode{se} @tab the standard error of kappa when there is no
## agreement beyond chance (Fleiss, Nee and Landis, 1979)
## @item @qcode{zval} @tab the statistic of the test that kappa is 0,
## @code{@var{kappa} / se}
## @item @qcode{p} @tab the two-sided p-value of that test
## @item @qcode{n} @tab the number of subjects
## @item @qcode{raters} @tab the number of raters @math{m}
## @item @qcode{counts} @tab the @math{n}-by-@math{c} matrix @code{N}
## @item @qcode{categories} @tab the categories, a numeric column or a cell
## array of character vectors
## @item @qcode{categorykappa} @tab the kappa of each category, the agreement
## on whether a subject belongs to it
## @item @qcode{categoryzval} @tab the statistic of the test that each
## category's kappa is 0
## @item @qcode{categoryp} @tab the two-sided p-value of that test
## @end multitable
##
## No confidence interval is given, there being no standard large-sample
## standard error for kappa away from 0.  Kappa is undefined when chance
## alone predicts complete agreement, as when every rating falls in one
## category, and such ratings are refused.  MATLAB has no Fleiss' kappa, so
## @code{fleissKappa} is specific to Octave.
##
## References:
## @enumerate
## @item
## J. L. Fleiss.  Measuring nominal scale agreement among many raters.
## @emph{Psychological Bulletin}, 76(5):378--382, 1971.
## @item
## J. L. Fleiss, J. C. M. Nee and J. R. Landis.  Large sample variance of
## kappa in the case of different sets of raters.  @emph{Psychological
## Bulletin}, 86(5):974--977, 1979.
## @item
## A. J. Conger.  Integration and generalization of kappas for multiple
## raters.  @emph{Psychological Bulletin}, 88(2):322--328, 1980.
## @end enumerate
##
## @seealso{cohenKappa, friedman}
## @end deftypefn

function [kappa, stats] = fleissKappa (X, varargin)

  if (nargin < 1)
    print_usage ();
  endif
  if (! (ismatrix (X) && ndims (X) == 2 && ! ischar (X)))
    error ("fleissKappa: X must be a matrix with one column per rater.");
  endif

  ## Parse optional Name-Value paired arguments
  optNames = {'Categories', 'Exact'};
  dfValues = {[], false};
  [cats, exact, rest] = parsePairedArguments (optNames, dfValues, varargin(:));
  if (! isempty (rest))
    if (mod (numel (rest), 2) == 0 && ischar (rest{1}))
      error ("fleissKappa: unknown parameter name '%s'.", rest{1});
    endif
    error ("fleissKappa: optional arguments must be Name-Value pairs.");
  endif
  if (! (isscalar (exact) && (islogical (exact) || isnumeric (exact))
         && any (exact == [0, 1])))
    error ("fleissKappa: 'Exact' must be a logical scalar.");
  endif

  m = columns (X);
  if (m < 2)
    error ("fleissKappa: X must hold the ratings of at least two raters.");
  endif
  R = cell (1, m);
  for j = 1:m
    R{j} = X(:,j);
  endfor
  [codes, cats, errmsg] = ratingcodes (R, cats);
  if (! isempty (errmsg))
    error ("fleissKappa: %s", errmsg);
  endif
  codes = codes(! any (isnan (codes), 2),:);
  n = rows (codes);
  if (n == 0)
    error ("fleissKappa: no subject was rated by every rater.");
  endif

  ## Count of raters per subject and category, and agreement
  c = numel (cats);
  N = zeros (n, c);
  for j = 1:m
    N += (codes(:,j) == (1:c));
  endfor
  pj = sum (N, 1) / (n * m);
  po = mean ((sum (N .^ 2, 2) - m) / (m * (m - 1)));
  pe = sum (pj .^ 2);
  if (exact)
    ## Conger: less the variance over raters of each rater's proportions
    Q = zeros (m, c);
    for j = 1:m
      Q(j,:) = sum (codes(:,j) == (1:c), 1) / n;
    endfor
    pe -= sum (var (Q, 1, 1)) / (m - 1);
  endif
  if (pe >= 1)
    error (strcat ("fleissKappa: kappa is undefined when chance predicts", ...
                   " complete agreement, as when every rating falls in one", ...
                   " category."));
  endif
  kappa = (po - pe) / (1 - pe);

  ## Tests under no agreement beyond chance (Fleiss, Nee and Landis, 1979)
  qj = 1 - pj;
  if (exact)
    se = NaN;
  else
    spq = sum (pj .* qj);
    se = sqrt (2 * (spq ^ 2 - sum (pj .* qj .* (qj - pj))) ...
               / (spq ^ 2 * n * m * (m - 1)));
  endif
  zval = kappa / se;
  p = 2 * normcdf (-abs (zval));

  ## Each category's kappa, against all the others
  pjk = (sum (N .^ 2, 1) - n * m * pj) ./ (n * m * (m - 1) * pj);
  ck = (pjk - pj) ./ qj;
  cz = ck / sqrt (2 / (n * m * (m - 1)));
  cp = 2 * normcdf (-abs (cz));

  stats = struct ('po', po, 'pe', pe, 'se', se, 'zval', zval, 'p', p, ...
                  'n', n, 'raters', m, 'counts', N, ...
                  'categories', {cats}, 'categorykappa', ck(:), ...
                  'categoryzval', cz(:), 'categoryp', cp(:));

endfunction

%!demo
%! ## Six psychiatrists diagnosed 30 patients (R's irr::diagnoses)
%! D = [4 4 4 4 4 4; 2 2 2 5 5 5; 2 3 3 3 3 5; 5 5 5 5 5 5; 2 2 2 4 4 4;
%!      1 1 3 3 3 3; 3 3 3 3 5 5; 1 1 3 3 3 4; 1 1 4 4 4 4; 5 5 5 5 5 5;
%!      1 4 4 4 4 4; 1 2 4 4 4 4; 2 2 2 3 3 3; 1 4 4 4 4 4; 2 2 4 4 4 5;
%!      3 3 3 3 3 5; 1 1 1 4 5 5; 1 1 1 1 1 2; 2 2 4 4 4 4; 1 3 3 5 5 5;
%!      5 5 5 5 5 5; 2 4 4 4 4 4; 2 2 4 5 5 5; 1 1 4 4 4 4; 1 4 4 4 4 5;
%!      2 2 2 2 2 4; 1 1 1 1 5 5; 2 2 4 4 4 4; 1 3 3 3 3 3; 5 5 5 5 5 5];
%! [kappa, stats] = fleissKappa (D)

## Test output
## Values checked against R 4.5.0, irr 0.85 (kappam.fleiss)
%!shared D
%! D = [4 4 4 4 4 4; 2 2 2 5 5 5; 2 3 3 3 3 5; 5 5 5 5 5 5; 2 2 2 4 4 4;
%!      1 1 3 3 3 3; 3 3 3 3 5 5; 1 1 3 3 3 4; 1 1 4 4 4 4; 5 5 5 5 5 5;
%!      1 4 4 4 4 4; 1 2 4 4 4 4; 2 2 2 3 3 3; 1 4 4 4 4 4; 2 2 4 4 4 5;
%!      3 3 3 3 3 5; 1 1 1 4 5 5; 1 1 1 1 1 2; 2 2 4 4 4 4; 1 3 3 5 5 5;
%!      5 5 5 5 5 5; 2 4 4 4 4 4; 2 2 4 5 5 5; 1 1 4 4 4 4; 1 4 4 4 4 5;
%!      2 2 2 2 2 4; 1 1 1 1 5 5; 2 2 4 4 4 4; 1 3 3 3 3 3; 5 5 5 5 5 5];
%!test
%! [k, s] = fleissKappa (D);
%! assert_equal (k, 0.430244520060141, 1e-14);
%! assert_equal (s.zval, 17.6518305829914, 1e-12);
%! assert_equal (s.p, 2 * normcdf (-17.6518305829914), 1e-80);
%! assert_equal (s.n, 30);
%! assert_equal (s.raters, 6);
%! assert_equal (s.categories, (1:5)');
%!test
%! [~, s] = fleissKappa (D);
%! assert_equal (s.categorykappa, [0.244755244755245; 0.244755244755245; ...
%!                                 0.52; 0.471127272727273; ...
%!                                 0.566117806823969], 1e-14);
%! assert_equal (s.categoryzval, [5.1920427989222; 5.1920427989222; ...
%!                                11.0308657865101; 9.99411868042136; ...
%!                                12.0091722046705], 1e-12);
%!test
%! [k, s] = fleissKappa (D, 'Exact', true);
%! assert_equal (k, 0.441808540329333, 1e-14);
%! assert_equal (s.zval, NaN);
%! assert_equal (s.p, NaN);
%!test
%! L = {'Depression', 'Personality', 'Schizophrenia', 'Neurosis', 'Other'};
%! assert_equal (fleissKappa (L(D)), 0.430244520060141, 1e-14);
%!test
%! L = {'Depression', 'Personality', 'Schizophrenia', 'Neurosis', 'Other'};
%! assert_equal (fleissKappa (categorical (L(D))), 0.430244520060141, 1e-14);
%!test
%! [k, s] = fleissKappa ([D; NaN, 1, 1, 1, 1, 1]);
%! assert_equal (k, 0.430244520060141, 1e-14);
%! assert_equal (s.n, 30);
%!test
%! [~, s] = fleissKappa (D, 'Categories', 1:6);
%! assert_equal (s.categorykappa(6), NaN);
%! assert_equal (size (s.counts), [30, 6]);
%!test
%! x = [1, 2, 3, 1, 2, 3, 2, 2, 1, 3, 3, 1, 2, 1, 3, 2, 2, 3, 1, 1]';
%! y = [1, 2, 2, 1, 3, 3, 2, 1, 1, 3, 2, 1, 2, 2, 3, 2, 3, 3, 1, 2]';
%! assert_equal (fleissKappa ([x, y], 'Exact', true), cohenKappa (x, y), 1e-14);

## Test input validation
%!error <Invalid call to fleissKappa> fleissKappa ()
%!error <fleissKappa: X must be a matrix with one column per rater.> ...
%! fleissKappa (['ab'; 'cd'])
%!error <fleissKappa: X must hold the ratings of at least two raters.> ...
%! fleissKappa ([1; 2; 3])
%!error <fleissKappa: unknown parameter name 'Foo'.> fleissKappa (D, 'Foo', 1)
%!error <fleissKappa: optional arguments must be Name-Value pairs.> ...
%! fleissKappa (D, 'Exact')
%!error <fleissKappa: 'Exact' must be a logical scalar.> ...
%! fleissKappa (D, 'Exact', 2)
%!error <fleissKappa: a rating is not one of 'Categories'.> ...
%! fleissKappa (D, 'Categories', 1:4)
%!error <fleissKappa: kappa is undefined when chance predicts complete agreement, as when every rating falls in one category.> ...
%! fleissKappa (ones (4, 3))
%!error <fleissKappa: no subject was rated by every rater.> ...
%! fleissKappa ([1, NaN; NaN, 2])
