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
## @deftypefn  {statistics} {@var{mmdval} =} mmdtest (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{mmdval} =} mmdtest (@var{X}, @var{Y}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{mmdval}, @var{p}] =} mmdtest (@dots{})
## @deftypefnx {statistics} {[@var{mmdval}, @var{p}, @var{h}] =} mmdtest (@dots{})
##
## Two-sample multivariate test on the maximum mean discrepancy.
##
## @code{@var{mmdval} = mmdtest (@var{X}, @var{Y})} returns the squared
## maximum mean discrepancy between the samples @var{X} and @var{Y}, whose rows
## are observations and whose columns are variables.  It measures how far
## apart the two distributions lie, 0 where the two samples are alike.  With
## @math{m} rows in @var{X}, @math{n} in @var{Y} and @math{k} the kernel, it is
## the biased estimate
## @tex
## $$\mathrm{MMD}^2 = {1 \over m^2} \sum_{i,j=1}^{m} k(x_i, x_j)
## - {2 \over m n} \sum_{i=1}^{m} \sum_{j=1}^{n} k(x_i, y_j)
## + {1 \over n^2} \sum_{i,j=1}^{n} k(y_i, y_j).$$
## @end tex
## @ifnottex
## @math{MMD^2 = sum (k(x_i, x_j)) / m^2 - 2 sum (k(x_i, y_j)) / (m n) +
## sum (k(y_i, y_j)) / n^2}.
## @end ifnottex
##
## @var{X} and @var{Y} are both numeric matrices with the same number of
## columns, or both tables, read over the variables they share, which must be
## all the variables of one of them.  An observation holding a missing value in
## a variable used is left out.
##
## The kernel is Gaussian, @math{k(u, v) = exp (-d^2 / s)}.  Each continuous
## variable is standardised by the mean and standard deviation of @var{X}, a
## standard deviation of 0 being taken as 1, and contributes its squared
## difference to @math{d^2}; a variable holding levels contributes 2 where the
## two observations differ, as one column per level would.  The scale @math{s}
## is the median of @math{d^2} over the pairs of observations in @var{X}, taken
## as 1 where that median is 0.
##
## @code{[@var{mmdval}, @var{p}, @var{h}] = mmdtest (@dots{})} also returns the
## p-value of the permutation test of the null hypothesis that @var{X} and
## @var{Y} come from the same distribution, and @var{h}, which is 1 where the
## null hypothesis is rejected at the significance level @qcode{'Alpha'} and 0
## otherwise.  Each permutation deals the pooled observations out again as
## @var{X} and @var{Y}, standardises and scales them afresh, and computes the
## statistic; @var{p} is the share of the permutations whose statistic is at
## least the observed one.  The permutations are drawn with @code{randperm},
## so the state of @code{rand} decides them.  They are drawn only where
## @var{p} or @var{h} is asked for.
##
## @code{@var{mmdval} = mmdtest (@dots{}, @var{Name}, @var{Value})} takes the
## following options.
##
## @multitable @columnfractions 0.26 0.02 0.72
## @headitem @var{Name} @tab @tab @var{Value}
## @item @qcode{'Alpha'} @tab @tab The significance level, a scalar between 0
## and 1, 0.05 by default.
## @item @qcode{'NumPermutations'} @tab @tab The number of permutations, a
## positive integer, 1000 by default.
## @item @qcode{'VariableNames'} @tab @tab The variables to use, among those
## the two tables share, as a character vector, a string array or a cell array
## of character vectors.  All the shared variables by default.
## @item @qcode{'CategoricalVariables'} @tab @tab The variables holding
## levels: @qcode{'all'}, their indices, a logical vector over the variables,
## or, for tables, their names.  A table variable holding logical values, an
## unordered @code{categorical} array, a @code{string} array or a cell array of
## character vectors holds levels whether named here or not.  A matrix holds
## none unless named.
## @item @qcode{'Options'} @tab @tab A structure as @code{statset} returns.
## Parallel computing and random streams are not implemented, and are refused
## where asked for.
## @end multitable
##
## Reference: A. Gretton, K. M. Borgwardt, M. J. Rasch, B. Schoelkopf and
## A. Smola (2012).  A kernel two-sample test.  Journal of Machine Learning
## Research, 13, 723-773.
##
## @seealso{knntest, kstest2}
## @end deftypefn

function [mmdval, p, h] = mmdtest (X, Y, varargin)

  ## Input validation
  if (nargin < 2)
    error ("mmdtest: too few input arguments.");
  endif
  optNames = {'Alpha', 'NumPermutations', 'VariableNames', ...
              'CategoricalVariables', 'Options'};
  dfValues = {0.05, 1000, [], [], []};
  [alpha, nperm, vnames, catvars, opts, args] = ...
                        parsePairedArguments (optNames, dfValues, varargin(:));
  if (! isempty (args))
    error ("mmdtest: invalid optional paired argument.");
  endif
  if (! (isnumeric (alpha) && isreal (alpha) && isscalar (alpha)
         && alpha > 0 && alpha < 1))
    error ("mmdtest: 'Alpha' must be a scalar between 0 and 1.");
  endif
  if (! (isnumeric (nperm) && isreal (nperm) && isscalar (nperm)
         && isfinite (nperm) && nperm >= 1 && nperm == fix (nperm)))
    error ("mmdtest: 'NumPermutations' must be a positive integer.");
  endif
  if (! isempty (opts))
    if (! isstruct (opts))
      error ("mmdtest: 'Options' must be a structure.");
    endif
    if (isfield (opts, 'UseParallel') && ! isempty (opts.UseParallel)
        && ! strcmpi (opts.UseParallel, 'never')
        && ! isequal (opts.UseParallel, false))
      error ("mmdtest: parallel computing is not implemented.");
    endif
    if (isfield (opts, 'Streams') && ! isempty (opts.Streams))
      error ("mmdtest: random streams are not implemented.");
    endif
  endif

  ## The pooled sample coded one column per variable, rows holding a missing
  ## value left out, and which variables hold levels
  [Z, iscat, mx, my, errmsg] = __twosample__ (X, Y, vnames, catvars);
  if (! isempty (errmsg))
    error ("mmdtest: %s", errmsg);
  endif
  m = mx + my;

  ## A difference in levels does not depend on how the sample is dealt out;
  ## a continuous difference is rescaled by each dealing's own X
  Dlev = zeros (m);
  for j = find (iscat)
    Dlev += 2 * (Z(:,j) != Z(:,j)');
  endfor
  cont = find (! iscat);
  G = cell (1, numel (cont));
  for j = 1:numel (cont)
    G{j} = (Z(:,cont(j)) - Z(:,cont(j))') .^ 2;
  endfor

  mmdval = mmdStat (Z(:,cont), G, Dlev, 1:mx, mx+1:m);
  if (nargout < 2)
    return;
  endif

  ## Permutations as extreme as the observed statistic, to within rounding
  count = 0;
  tol = 1e-12 * abs (mmdval);
  for i = 1:nperm
    r = randperm (m);
    count += mmdStat (Z(:,cont), G, Dlev, r(1:mx), r(mx+1:end)) >= mmdval - tol;
  endfor
  p = count / nperm;
  h = double (p <= alpha);

endfunction

## The squared maximum mean discrepancy with the rows A as X and B as Y.
function v = mmdStat (C, G, Dlev, a, b)

  D2 = Dlev;
  sd = std (C(a,:), [], 1);
  sd(sd == 0 | ! isfinite (sd)) = 1;
  for j = 1:numel (G)
    D2 += G{j} / sd(j) ^ 2;
  endfor
  Dx = D2(a,a);
  s = median (Dx(triu (true (numel (a)), 1)));
  if (! (s > 0))
    s = 1;
  endif
  K = exp (-D2 / s);
  v = mean (K(a,a)(:)) - 2 * mean (K(a,b)(:)) + mean (K(b,b)(:));

endfunction

## Expected values from MATLAB R2026a
%!assert_equal (mmdtest ([0; 1], 0.5), 0.126338154442911, -1e-13)
%!assert_equal (mmdtest ([0; 1], 1.5), 0.799739712952452, -1e-13)
%!assert_equal (mmdtest ([0; 1], 3), 1.665500671892900, -1e-13)
%!assert_equal (mmdtest ([0; 1], 10), 1.683939720585721, -1e-13)
%!assert_equal (mmdtest ([0; 1], 1), 0.316060279414279, -1e-13)
%!assert_equal (mmdtest ([0; 2], 1), 0.126338154442911, -1e-13)
%!assert_equal (mmdtest ([0; 10], 5), 0.126338154442911, -1e-13)
%!assert_equal (mmdtest ([0; 1; 4], 2), 0.200029287193120, -1e-13)
%!assert_equal (mmdtest ([0; 1; 4], [2; 6]), 0.269711792340645, -1e-13)
%!assert_equal (mmdtest ([2; 6], [0; 1; 4]), 0.237458896831246, -1e-13)
%!assert_equal (mmdtest ([0; 1; NaN], 1.5), 0.799739712952452, -1e-13)
%!assert_equal (mmdtest ([0, 0; 1, 10], [1, 0]), 0.470878401160454, -1e-13)
%!assert_equal (mmdtest ([0, 0; 1, 10; 3, 2], [1, 0; 2, 5]), ...
%!              0.102827836885193, -1e-13)
%!assert_equal (mmdtest ([0, 0; 1, 10; 3, 2], [1, 0; 2, 5; 4, 4; 0, 9]), ...
%!              0.029108941526544, -1e-12)
%!assert_equal (mmdtest ([1; 2; 1], [2; 2], 'CategoricalVariables', 1), ...
%!              0.561884941180940, -1e-13)
%!assert_equal (mmdtest ([1; 2; 1; 3], [2; 2; 3], ...
%!                      'CategoricalVariables', 1), 0.272163018384518, -1e-13)
%!assert_equal (mmdtest ([1, 1; 2, 1; 1, 2], [2, 2; 2, 1], ...
%!                      'CategoricalVariables', 'all'), ...
%!              0.419413238496425, -1e-13)
%!assert_equal (mmdtest ([1, 0; 2, 1; 1, 4], [2, 2; 1, 6], ...
%!                      'CategoricalVariables', 1), 0.266631734445491, -1e-13)
%!test
%! ## Six pairs in X: the median of the squared distances, not its square
%! v = mmdtest ([1, 0; 2, 1; 1, 4; 2, 3], [2, 2; 1, 6; 3, 1], ...
%!              'CategoricalVariables', 1);
%! assert_equal (v, 0.149533632990419, -1e-13);
%!test
%! ## A table's categorical variable holds levels without being named
%! v = mmdtest (table (categorical ([1; 2; 1])), table (categorical ([2; 2])));
%! assert_equal (v, 0.561884941180940, -1e-13);
%!test
%! ## A standard deviation of 0 is taken as 1, and so is a scale of 0
%! assert_equal (mmdtest ([1, 0; 1, 1; 1, 2], [2, 0; 3, 1]), ...
%!               0.887991689349540, -1e-13);
%! assert_equal (mmdtest ([1; 1; 1], [2; 3]), 1.297744640525545, -1e-13);
%!test
%! load fisheriris
%! assert_equal (mmdtest (meas(1:25,:), meas(26:50,:)), ...
%!               0.028159364150432, -1e-12);
%! assert_equal (mmdtest (meas(51:100,:), meas(101:150,:)), ...
%!               0.540247447946648, -1e-12);
%! assert_equal (mmdtest (meas(101:150,:), meas(51:100,:)), ...
%!               0.604761188982786, -1e-12);
%!test
%! ## Each permutation is standardised afresh; the exact p-value over all 70
%! ## dealings is 0.2286, and MATLAB gave 0.2303 and 0.2273 over 20000
%! rand ('seed', 1);
%! [v, p] = mmdtest ([2.5; 4; 3; 2], [3; -8; 2.5; 0], 'NumPermutations', 4000);
%! assert_equal (v, 0.249999985933103, -1e-13);
%! assert_equal (abs (p - 0.2286) < 0.03, true);
%!test
%! rand ('seed', 1);
%! [~, p, h] = mmdtest ([0; 0.1; 0.2], [10; 10.1; 10.2], ...
%!                      'NumPermutations', 9, 'Alpha', 0.2);
%! assert_equal (p * 9, round (p * 9));
%! assert_equal (h, double (p <= 0.2));
%!assert_equal (mmdtest ([0; 1], 0.5, 'Options', ...
%!                      statset ('UseParallel', false)), ...
%!              0.126338154442911, -1e-13)

%!error<mmdtest: too few input arguments.> mmdtest (1)
%!error<mmdtest: invalid optional paired argument.> mmdtest (1, 2, 'Tail', 1)
%!error<mmdtest: 'Alpha' must be a scalar between 0 and 1.> ...
%! mmdtest ([0; 1], 2, 'Alpha', 0)
%!error<mmdtest: 'NumPermutations' must be a positive integer.> ...
%! mmdtest ([0; 1], 2, 'NumPermutations', 0)
%!error<mmdtest: 'Options' must be a structure.> ...
%! mmdtest ([0; 1], 2, 'Options', 1)
%!error<mmdtest: parallel computing is not implemented.> ...
%! mmdtest ([0; 1], 2, 'Options', struct ('UseParallel', true))
%!error<mmdtest: random streams are not implemented.> ...
%! mmdtest ([0; 1], 2, 'Options', struct ('Streams', 1))
%!error<mmdtest: X and Y must have the same number of columns.> ...
%! mmdtest ([0, 1; 1, 2], 2)
%!error<mmdtest: X and Y must each hold an observation with no missing value.> ...
%! mmdtest (NaN, 2)
