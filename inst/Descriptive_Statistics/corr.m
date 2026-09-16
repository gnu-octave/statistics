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
## @deftypefn  {statistics} {@var{rho} =} corr (@var{x})
## @deftypefnx {statistics} {@var{rho} =} corr (@var{x}, @var{y})
## @deftypefnx {statistics} {@var{rho} =} corr (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{rho}, @var{pval}] =} corr (@dots{})
##
## Linear or rank correlation coefficients.
##
## @code{@var{rho} = corr (@var{x})} returns the matrix of pairwise correlation
## coefficients between the columns of @var{x}, whose rows are observations and
## whose columns are variables.  For an @math{n}-by-@math{k} matrix @var{x},
## @var{rho} is @math{k}-by-@math{k} and @code{@var{rho}(i,j)} is the
## correlation between the @math{i}-th and the @math{j}-th column of @var{x}.
##
## @code{@var{rho} = corr (@var{x}, @var{y})} returns the correlations between
## the columns of @var{x} and the columns of @var{y}.  For an
## @math{n}-by-@math{k1} matrix @var{x} and an @math{n}-by-@math{k2} matrix
## @var{y}, @var{rho} is @math{k1}-by-@math{k2} and @code{@var{rho}(i,j)} is the
## correlation between the @math{i}-th column of @var{x} and the @math{j}-th
## column of @var{y}.  @var{x} and @var{y} must have the same number of rows.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.2 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Type'} @tab The coefficient to compute: @qcode{'Pearson'}
## (default) for the linear correlation coefficient, @qcode{'Kendall'} for
## Kendall's tau-b, or @qcode{'Spearman'} for Spearman's rho.
##
## @item @qcode{'Rows'} @tab How missing values are handled: @qcode{'all'}
## (default) uses every row, so an entry is @qcode{NaN} whenever either of its
## columns holds a @qcode{NaN}; @qcode{'complete'} first removes every row
## holding a @qcode{NaN} in any column; @qcode{'pairwise'} computes each entry
## from the rows where its own pair of columns is present.
##
## @item @qcode{'Tail'} @tab The alternative hypothesis of the test reported in
## @var{pval}: @qcode{'both'} (default) for a correlation different from zero,
## @qcode{'right'} for a positive one, or @qcode{'left'} for a negative one.
##
## @item @qcode{'Weights'} @tab A column vector of @math{n} nonnegative
## observation weights.  The default weights every observation equally.
## @end multitable
##
## Option values are matched without regard to case, and an unambiguous
## abbreviation is accepted, so @qcode{'spear'} selects @qcode{'Spearman'}.
##
## @code{[@var{rho}, @var{pval}] = corr (@dots{})} also returns @var{pval}, the
## p-value of a test of the null hypothesis that the corresponding correlation
## is zero, against the alternative named by @qcode{'Tail'}.  A @var{pval} entry
## is @qcode{NaN} wherever its @var{rho} entry is, and every entry is
## @qcode{NaN} when @qcode{'Weights'} is given.
##
## The p-value is computed as follows:
##
## @itemize
## @item @qcode{'Pearson'}: from a Student's t distribution with @math{n-2}
## degrees of freedom applied to
## @code{@var{t} = @var{rho} * sqrt ((n - 2) / (1 - @var{rho}^2))}, which is
## exact when the data are normally distributed.
##
## @item @qcode{'Kendall'}: from the exact permutation distribution when
## @math{n < 10}, and for a sample of fewer than 50 observations holding no
## tied values; otherwise from a normal approximation with a continuity
## correction and the usual correction for ties.
##
## @item @qcode{'Spearman'}: from the exact permutation distribution when
## @math{n < 10}; from the approximation of Best and Roberts (algorithm AS 89)
## for a larger sample holding no tied values; and from the same Student's t
## transformation as the Pearson coefficient for a larger sample that does
## hold them.
## @end itemize
##
## The exact p-value is the proportion of the @math{n!} orderings of one
## variable giving a coefficient as extreme as the observed one, so tied
## values need no special treatment there.
##
## A column of constant values has no defined correlation, so its entries,
## the diagonal one included, are @qcode{NaN}.  @qcode{Inf} makes a Pearson
## coefficient @qcode{NaN}, while the rank coefficients order it like any other
## value.
##
## Two deviations from MATLAB.  A character array is refused, where MATLAB
## correlates the character codes; its acceptance there is incidental,
## correlation having no meaning for text.  An array of more than two
## dimensions is refused, where MATLAB refuses most of them from inside a
## matrix multiplication but flattens a leading singleton dimension and
## answers @qcode{NaN}; core Octave's @code{corr} answers @qcode{NaN} there
## too.
##
## References:
## @enumerate
## @item
## M. G. Kendall.  A new measure of rank correlation.  @emph{Biometrika},
## 30(1-2):81--93, 1938.
## @item
## D. J. Best and D. E. Roberts.  Algorithm AS 89: the upper tail probabilities
## of Spearman's rho.  @emph{Applied Statistics}, 24(3):377--379, 1975.
## @end enumerate
##
## @seealso{corrcoef, cov, partialcorr, tiedrank, kendall, spearman}
## @end deftypefn

function [rho, pval] = corr (x, varargin)

  ## Input validation
  if (nargin < 1)
    print_usage ();
  endif

  ## Y is present unless the second argument opens the Name-Value pairs
  y = [];
  haveY = false;
  if (numel (varargin) > 0 ...
      && ! (ischar (varargin{1}) || isstring (varargin{1})))
    y = varargin{1};
    varargin(1) = [];
    haveY = true;
  endif

  if (! (isnumeric (x) || islogical (x)))
    error ("corr: X must be a numeric or logical matrix.");
  endif
  if (ndims (x) > 2)
    error ("corr: X must be a matrix or a vector.");
  endif
  issingle = isa (x, 'single');
  x = double (x);
  if (columns (x) == 0)
    error ("corr: X must have at least one column.");
  endif
  n = rows (x);
  if (haveY)
    if (! (isnumeric (y) || islogical (y)))
      error ("corr: Y must be a numeric or logical matrix.");
    endif
    if (ndims (y) > 2)
      error ("corr: Y must be a matrix or a vector.");
    endif
    issingle = issingle || isa (y, 'single');
    y = double (y);
    if (rows (y) != n)
      error ("corr: X and Y must have the same number of rows.");
    endif
  endif

  ## Parse optional Name-Value paired arguments
  optNames = {'Type', 'Rows', 'Tail', 'Weights'};
  dfValues = {'Pearson', 'all', 'both', []};
  [type, rowopt, tail, w, args] = parsePairedArguments (optNames, dfValues, ...
                                                        varargin(:));
  if (! isempty (args))
    error (strcat ("corr: optional arguments must be 'Type', 'Rows',", ...
                   " 'Tail' or 'Weights' Name-Value pairs."));
  endif
  type = matchopt_ (type, {'Pearson', 'Kendall', 'Spearman'}, 'Type');
  rowopt = matchopt_ (rowopt, {'all', 'complete', 'pairwise'}, 'Rows');
  tail = matchopt_ (tail, {'both', 'right', 'left'}, 'Tail');
  if (! isempty (w))
    if (! (isa (w, 'double') || isa (w, 'single')) || ! isreal (w) ...
        || ! iscolumn (w) || numel (w) != n)
      error ("corr: 'Weights' must be a column vector of %d real values.", n);
    endif
    if (any (w < 0))
      error ("corr: 'Weights' must be nonnegative.");
    endif
    w = double (w);
  endif

  ## Both forms are computed as X against Y, the one-argument form against
  ## itself, which is what makes its result symmetric with a unit diagonal.
  if (haveY)
    b = y;
  else
    b = x;
  endif
  if (columns (b) == 0)
    rho = zeros (columns (x), 0);
    pval = zeros (columns (x), 0);
    return;
  endif

  ## Drop the rows holding a missing value before anything else
  if (strcmp (rowopt, 'complete'))
    keep = ! any (isnan ([x, b]), 2);
    x = x(keep,:);
    b = b(keep,:);
    if (! isempty (w))
      w = w(keep);
    endif
  endif

  [rho, pval] = corr_compute_ (x, b, ! haveY, type, rowopt, tail, w, ...
                               nargout > 1);

  if (issingle)
    rho = single (rho);
    pval = single (pval);
  endif

endfunction

## Resolve an option value against the allowed spellings, without regard to
## case and accepting an unambiguous abbreviation.
function out = matchopt_ (val, allowed, name)
  if (isstring (val) && isscalar (val))
    val = char (val);
  endif
  if (! (ischar (val) && isrow (val)))
    val = '';
  endif
  idx = strncmpi (val, allowed, numel (val));
  if (isempty (val) || sum (idx) != 1)
    error ("corr: '%s' must be %s.", name, list_ (allowed));
  endif
  out = allowed{idx};
endfunction

## Render the allowed spellings as "'a', 'b' or 'c'"
function s = list_ (allowed)
  s = sprintf ("'%s', ", allowed{1:end-1});
  s = sprintf ("%s or '%s'", s(1:end-2), allowed{end});
endfunction

## The coefficients and their p-values, column pair by column pair
function [rho, pval] = corr_compute_ (x, b, symmetric, type, rowopt, tail, ...
                                      w, needp)
  k1 = columns (x);
  k2 = columns (b);
  rho = NaN (k1, k2);
  pval = NaN (k1, k2);
  pairwise = strcmp (rowopt, 'pairwise');
  for i = 1:k1
    if (symmetric)
      jj = i:k2;
    else
      jj = 1:k2;
    endif
    for j = jj
      a = x(:,i);
      c = b(:,j);
      wij = w;
      if (pairwise)
        keep = ! isnan (a) & ! isnan (c);
        a = a(keep);
        c = c(keep);
        if (! isempty (wij))
          wij = wij(keep);
        endif
      endif
      [r, p] = corr_pair_ (a, c, type, tail, wij, needp);
      rho(i,j) = r;
      if (needp)
        pval(i,j) = p;
      endif
    endfor
  endfor
  if (symmetric)
    rho = triu (rho) + triu (rho, 1)';
    ## A column correlates perfectly with itself unless it is degenerate,
    ## which the arithmetic above already reports as NaN.
    dd = find (! isnan (diag (rho)));
    di = sub2ind (size (rho), dd, dd);
    rho(di) = 1;
    if (needp)
      pval = triu (pval) + triu (pval, 1)';
      ## A weighted correlation has no p-value at all, the diagonal included
      if (isempty (w))
        pval(di) = 1;
      endif
    endif
  endif
endfunction

## One coefficient and its p-value
function [r, p] = corr_pair_ (a, c, type, tail, w, needp)
  n = numel (a);
  p = NaN;
  switch (type)
    case 'Pearson'
      r = wpearson_ (a, c, w);
      if (needp && isempty (w))
        p = pearson_pval_ (r, n, tail);
      endif
    case 'Spearman'
      if (isempty (w))
        ra = tiedrank (a);
        rc = tiedrank (c);
      else
        ra = wtiedrank_ (a, w);
        rc = wtiedrank_ (c, w);
      endif
      r = wpearson_ (ra, rc, w);
      if (needp && isempty (w))
        p = spearman_pval_ (a, c, ra, rc, r, tail);
      endif
    case 'Kendall'
      [r, S, V] = wkendall_ (a, c, w);
      if (needp && isempty (w))
        p = kendall_pval_ (a, c, r, S, V, tail);
      endif
  endswitch
endfunction

## Pearson coefficient, weighted when W is given
function r = wpearson_ (a, c, w)
  if (isempty (w))
    ac = a - mean (a);
    cc = c - mean (c);
    r = (ac' * cc) / sqrt ((ac' * ac) * (cc' * cc));
  else
    w = w / sum (w);
    ac = a - sum (w .* a);
    cc = c - sum (w .* c);
    r = sum (w .* ac .* cc) / sqrt (sum (w .* ac .^ 2) * sum (w .* cc .^ 2));
  endif
endfunction

## Midranks under observation weights: an observation of weight two occupies
## two consecutive ranks, so replicating it and ranking gives the same value.
function r = wtiedrank_ (x, w)
  [xs, idx] = sort (x(:));
  ws = w(idx);
  hi = cumsum (ws);
  lo = hi - ws;
  rs = (lo + hi) / 2;
  ## Tied values share the midrank of the block they form
  k = 1;
  n = numel (xs);
  while (k <= n)
    j = k;
    while (j < n && xs(j+1) == xs(k))
      j++;
    endwhile
    if (j > k)
      rs(k:j) = (lo(k) + hi(j)) / 2;
    endif
    k = j + 1;
  endwhile
  r = zeros (size (rs));
  r(idx) = rs;
endfunction

## Kendall's tau-b with its score S and the variance of S under the null
function [t, S, V] = wkendall_ (a, c, w)
  n = numel (a);
  sa = sign (a - a.');
  sc = sign (c - c.');
  if (isempty (w))
    S = sum (sum (sa .* sc)) / 2;
    na = sum (sum (sa != 0)) / 2;
    nc = sum (sum (sc != 0)) / 2;
  else
    W = w * w.';
    S = sum (sum (W .* sa .* sc)) / 2;
    na = sum (sum (W .* (sa != 0))) / 2;
    nc = sum (sum (W .* (sc != 0))) / 2;
  endif
  t = S / sqrt (na * nc);
  ta = tiecounts_ (a);
  tc = tiecounts_ (c);
  V = (n * (n - 1) * (2 * n + 5) - sum (ta .* (ta - 1) .* (2 * ta + 5)) ...
       - sum (tc .* (tc - 1) .* (2 * tc + 5))) / 18 ...
      + sum (ta .* (ta - 1) .* (ta - 2)) * sum (tc .* (tc - 1) .* (tc - 2)) ...
        / (9 * n * (n - 1) * (n - 2)) ...
      + sum (ta .* (ta - 1)) * sum (tc .* (tc - 1)) / (2 * n * (n - 1));
endfunction

## Sizes of the tied groups of X
function t = tiecounts_ (x)
  [~, ~, j] = unique (x(:));
  t = accumarray (j, 1);
endfunction

function tf = hasties_ (x)
  tf = numel (unique (x)) < numel (x);
endfunction

## Student's t p-value of a Pearson coefficient
function p = pearson_pval_ (r, n, tail)
  if (n < 3 || isnan (r))
    p = NaN;
    return;
  endif
  ## The two-sided probability comes from the incomplete beta function
  ## directly, which keeps its precision in the far tail where doubling a
  ## computed t quantile does not.
  df = n - 2;
  t = r * sqrt (df / (1 - r ^ 2));
  two = betainc (df / (df + t ^ 2), df / 2, 0.5);
  switch (tail)
    case 'both'
      p = two;
    case 'right'
      p = ifelse_ (t >= 0, two / 2, 1 - two / 2);
    case 'left'
      p = ifelse_ (t >= 0, 1 - two / 2, two / 2);
  endswitch
endfunction

function out = ifelse_ (cond, a, b)
  if (cond)
    out = a;
  else
    out = b;
  endif
endfunction

## Spearman p-value: exact below ten observations, else AS 89 without ties
## and the t transformation with them
function p = spearman_pval_ (a, c, ra, rc, r, tail)
  n = numel (a);
  if (n < 3 || isnan (r))
    p = NaN;
  elseif (n < 10)
    P = perms (1:n);
    p = perm_pval_ (rc(P) * ra, rc' * ra, tail);
  elseif (hasties_ (a) || hasties_ (c))
    p = pearson_pval_ (r, n, tail);
  else
    S = (1 - r) * (n ^ 3 - n) / 6;
    lower = as89_ (n, S);
    upper = 1 - as89_ (n, S + 2);
    switch (tail)
      case 'both'
        p = min (1, 2 * min (upper, lower));
      case 'right'
        p = upper;
      case 'left'
        p = lower;
    endswitch
  endif
endfunction

## Kendall p-value: exact below ten observations and for a smaller untied
## sample, else a normal approximation with a continuity correction
function p = kendall_pval_ (a, c, t, S, V, tail)
  n = numel (a);
  ties = hasties_ (a) || hasties_ (c);
  if (n < 3 || isnan (t))
    p = NaN;
  elseif (n < 10)
    p = perm_pval_ (kendall_scores_ (a, c), S, tail);
  elseif (! ties && n < 50)
    d = kendall_dist_ (n);
    p = disc_pval_ (d, t, tail);
  else
    switch (tail)
      case 'both'
        p = 2 * normcdf (-max ((abs (S) - 1) / sqrt (V), 0));
      case 'right'
        p = normcdf (-(S - 1) / sqrt (V));
      case 'left'
        p = normcdf ((S + 1) / sqrt (V));
    endswitch
    p = min (p, 1);
  endif
endfunction

## Kendall's score S under every ordering of C against A
function S = kendall_scores_ (a, c)
  n = numel (a);
  P = perms (1:n);
  A = sign (a - a.');
  C = sign (c - c.');
  S = zeros (rows (P), 1);
  for i = 1:n-1
    for j = i+1:n
      if (A(i,j) != 0)
        S += A(i,j) * C(sub2ind ([n, n], P(:,i), P(:,j)));
      endif
    endfor
  endfor
endfunction

## Tail probabilities of a statistic enumerated over every ordering
function p = perm_pval_ (S, s0, tail)
  tol = 1e-9 * max (1, abs (s0));
  upper = mean (S >= s0 - tol);
  lower = mean (S <= s0 + tol);
  switch (tail)
    case 'both'
      p = min (1, 2 * min (upper, lower));
    case 'right'
      p = upper;
    case 'left'
      p = lower;
  endswitch
endfunction

## Tail probabilities of a discrete null distribution given as [value, prob]
function p = disc_pval_ (d, r, tail)
  upper = sum (d(d(:,1) >= r - 1e-12, 2));
  lower = sum (d(d(:,1) <= r + 1e-12, 2));
  switch (tail)
    case 'both'
      p = min (1, 2 * min (upper, lower));
    case 'right'
      p = upper;
    case 'left'
      p = lower;
  endswitch
endfunction

## Null distribution of Kendall's tau for N untied observations, from the
## distribution of the number of inversions of a random permutation
function d = kendall_dist_ (n)
  c = 1;
  for k = 1:n-1
    c = conv (c, ones (1, k + 1));
  endfor
  K = n * (n - 1) / 2;
  d = [(2 * (0:K)' / K) - 1, (c / sum (c))'];
endfunction

## Upper tail probability of Spearman's S by algorithm AS 89
function q = as89_ (n, S)
  c = [0.2274, 0.2531, 0.1745, 0.0758, 0.1033, 0.3932, 0.0879, 0.0151, ...
       0.0072, 0.0831, 0.0131, 4.6e-4];
  b = 1 / n;
  x = (6 * (S - 1) * b / (n ^ 2 - 1) - 1) * sqrt (1 / b - 1);
  y = x * x;
  u = x * b * (c(1) + b * (c(2) + c(3) * b) ...
       + y * (-c(4) + b * (c(5) + c(6) * b) - y * b * (c(7) + c(8) * b ...
       - y * (c(9) - c(10) * b + y * b * (c(11) - c(12) * y)))));
  q = min (1, max (0, u / exp (y / 2) + normcdf (x, 0, 1, 'upper')));
endfunction

%!demo
%! ## Correlation between the columns of a matrix, with p-values
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! [rho, pval] = corr (x)

%!demo
%! ## Spearman's rank correlation, which a monotone relation makes exact
%! x = [1; 2; 3; 4; 5];
%! y = [1; 4; 9; 16; 25];
%! [rho, pval] = corr (x, y, 'Type', 'Spearman')

## Test output
## RHO and PVAL checked against corr in MATLAB R2024a Update 4
## (Statistics and Machine Learning Toolbox 24.1)
%!test  # Pearson, one matrix
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! [rho, pval] = corr (x);
%! assert_equal (rho, [1, 0.80516104831610558; 0.80516104831610558, 1], 1e-14);
%! assert_equal (pval, [1, 0.10016803410643096; 0.10016803410643096, 1], 1e-14);
%!test  # Pearson, two matrices
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! y = [2 1; 5 4; 4 3; 6 9; 8 7];
%! [rho, pval] = corr (x, y);
%! assert_equal (rho, [0.94518905671890663, 0.87745098039215708; ...
%!                     0.79999999999999993, 0.98019605881960692], 1e-14);
%! assert_equal (pval, [0.015276771734465051, 0.050541693400160549; ...
%!                      0.10408803866182803, 0.0033355462806318316], 1e-14);
%!test  # the one-sided tails of the Pearson test
%! x = [1; 3; 4; 7; 9];
%! y = [2; 5; 4; 8; 6];
%! [~, pboth] = corr (x, y);
%! [~, pright] = corr (x, y, 'Tail', 'right');
%! [~, pleft] = corr (x, y, 'Tail', 'left');
%! assert_equal (pboth, 0.10016803410643098, 1e-14);
%! assert_equal (pright, 0.050084017053215489, 1e-14);
%! assert_equal (pleft, 0.94991598294678448, 1e-14);
%!test  # Kendall's tau-b, exact below ten observations
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! [rho, pval] = corr (x, 'Type', 'Kendall');
%! assert_equal (rho, [1, 0.6; 0.6, 1], 1e-14);
%! assert_equal (pval, [1, 0.23333333333333334; 0.23333333333333334, 1], 1e-14);
%!test  # Spearman's rho, exact below ten observations
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! [rho, pval] = corr (x, 'Type', 'Spearman');
%! assert_equal (rho, [1, 0.79999999999999993; 0.79999999999999993, 1], 1e-14);
%! assert_equal (pval, [1, 0.13333333333333333; 0.13333333333333333, 1], 1e-14);
%!test  # the one-sided tails of the exact rank tests
%! x = [1; 3; 4; 7; 9];
%! y = [2; 5; 4; 8; 6];
%! [~, pk] = corr (x, y, 'Type', 'Kendall', 'Tail', 'right');
%! [~, ps] = corr (x, y, 'Type', 'Spearman', 'Tail', 'left');
%! assert_equal (pk, 0.11666666666666667, 1e-14);
%! assert_equal (ps, 0.95833333333333337, 1e-14);
%!test  # tied values, where the exact p-value needs no correction
%! x = [1 1; 2 1; 2 3; 4 3; 5 6; 5 6];
%! [rhok, pk] = corr (x, 'Type', 'Kendall');
%! [rhos, ps] = corr (x, 'Type', 'Spearman');
%! assert_equal (rhok(1,2), 0.88070484592797926, 1e-14);
%! assert_equal (pk(1,2), 0.044444444444444446, 1e-14);
%! assert_equal (rhos(1,2), 0.92318618234499539, 1e-14);
%! assert_equal (ps(1,2), 0.044444444444444446, 1e-14);
%!test  # nine observations, the largest sample the exact route covers
%! x = (1:9)';
%! y = [3; 1; 4; 6; 5; 9; 2; 8; 7];
%! [rk, pk] = corr (x, y, 'Type', 'Kendall');
%! [rs, ps] = corr (x, y, 'Type', 'Spearman');
%! assert_equal (pk, 0.11943893298059964, 1e-14);
%! assert_equal (ps, 0.096797839506172836, 1e-14);
%!test  # twenty untied observations: Kendall is still exact, Spearman is AS 89
%! x = (1:20)';
%! y = mod ((1:20) * 7, 101)';
%! [~, pk] = corr (x, y, 'Type', 'Kendall');
%! [~, ps] = corr (x, y, 'Type', 'Spearman');
%! assert_equal (pk, 0.098330218734756586, 1e-12);
%! assert_equal (ps, 0.65798513768943223, 1e-12);
%!test  # sixty untied observations: the normal approximation takes over
%! x = (1:60)';
%! y = mod ((1:60) * 7, 101)';
%! [~, pk] = corr (x, y, 'Type', 'Kendall');
%! [~, pkr] = corr (x, y, 'Type', 'Kendall', 'Tail', 'right');
%! assert_equal (pk, 0.15127746350184884, 1e-12);
%! assert_equal (pkr, 0.075638731750924421, 1e-12);
%!test  # AS 89 and its one-sided tails
%! x = (1:60)';
%! y = mod ((1:60) * 7, 101)';
%! [~, ps] = corr (x, y, 'Type', 'Spearman');
%! [~, psr] = corr (x, y, 'Type', 'Spearman', 'Tail', 'right');
%! [~, psl] = corr (x, y, 'Type', 'Spearman', 'Tail', 'left');
%! assert_equal (ps, 0.49214166036069795, 1e-12);
%! assert_equal (psr, 0.24607083018034898, 1e-12);
%! assert_equal (psl, 0.75406296187189426, 1e-12);
%!test  # tied values above the exact route: normal for tau, Student t for rho
%! x = (1:30)';
%! y = mod ((1:30) * 7, 5)';
%! [~, pk] = corr (x, y, 'Type', 'Kendall');
%! [~, ps] = corr (x, y, 'Type', 'Spearman');
%! assert_equal (pk, 0.67454343834388397, 1e-12);
%! assert_equal (ps, 0.66780132080922749, 1e-12);
%!test  # 'Rows' 'complete' drops every row holding a NaN
%! x = [1 2; NaN 5; 4 4; 7 NaN; 9 6];
%! [rho, pval] = corr (x, 'Rows', 'complete');
%! assert_equal (rho, [1, 0.98974331861078702; 0.98974331861078702, 1], 1e-14);
%! assert_equal (pval, ...
%!               [1, 0.091257896685979945; 0.091257896685979945, 1], 1e-14);
%!test  # 'Rows' 'pairwise' uses the rows each pair of columns shares
%! x = [1 2; NaN 5; 4 4; 7 NaN; 9 6];
%! y = [2 1; 5 4; 4 3; 6 9; 8 7];
%! [rho, pval] = corr (x, y, 'Rows', 'pairwise');
%! assert_equal (rho, [0.99591000331047852, 0.88678890262741183; ...
%!                     0.95638207148956245, 0.95638207148956245], 1e-14);
%! assert_equal (pval, [0.0040899966895214749, 0.11321109737258821; ...
%!                      0.043617928510437588, 0.043617928510437588], 1e-14);
%!test  # the default keeps every row, so a NaN reaches its whole row and column
%! x = [1 2; NaN 5; 4 4; 7 8; 9 6];
%! rho = corr (x);
%! assert_equal (rho, [NaN, NaN; NaN, 1], 1e-14);
%!test  # 'Weights' gives a weighted Pearson coefficient
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! w = [0.5; 1.25; 2; 0.75; 1];
%! rho = corr (x, 'Weights', w);
%! assert_equal (rho(1,2), 0.7498537417425275, 1e-14);
%!test  # 'Weights' applies to the rank coefficients too
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! w = [0.5; 1.25; 2; 0.75; 1];
%! rk = corr (x, 'Weights', w, 'Type', 'Kendall');
%! rs = corr (x, 'Weights', w, 'Type', 'Spearman');
%! assert_equal (rk(1,2), 0.43169398907103823, 1e-14);
%! assert_equal (rs(1,2), 0.63438256658595638, 1e-14);
%!test  # a weighted correlation has no p-value
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! [~, pval] = corr (x, 'Weights', [0.5; 1.25; 2; 0.75; 1]);
%! assert_equal (isnan (pval), true (2));
%!test  # a constant column has no correlation, its diagonal entry included
%! x = [1 4; 3 4; 4 4; 7 4; 9 4];
%! [rho, pval] = corr (x);
%! assert_equal (rho, [1, NaN; NaN, NaN], 1e-14);
%! assert_equal (pval, [1, NaN; NaN, NaN], 1e-14);
%!test  # two observations leave no degrees of freedom
%! [rho, pval] = corr ([1; 2], [1; 3]);
%! assert_equal (rho, 0.99999999999999989, 1e-14);
%! assert_equal (isnan (pval), true);
%!test  # single input gives a single result
%! x = single ([1 2; 3 5; 4 4; 7 8; 9 6]);
%! [rho, pval] = corr (x);
%! assert_equal (class (rho), 'single');
%! assert_equal (class (pval), 'single');
%! assert_equal (rho, single ([1, 0.80516105890274048; ...
%!                             0.80516105890274048, 1]), single (1e-7));
%!test  # logical input is correlated as the numbers it stands for
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! rho = corr (x > 4);
%! assert_equal (rho(1,2), 0.66666666666666641, 1e-14);
%!test  # an option value may be abbreviated and is matched without case
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! assert_equal (corr (x, 'type', 'spear'), corr (x, 'Type', 'Spearman'));
%! assert_equal (corr (x, 'Type', 'P'), corr (x));
%!test  # a Y without columns gives a result without columns
%! x = [1 2; 3 5; 4 4; 7 8; 9 6];
%! [rho, pval] = corr (x, zeros (5, 0));
%! assert_equal (size (rho), [2, 0]);
%! assert_equal (size (pval), [2, 0]);
%!test  # a row vector is one observation of several variables, as in MATLAB
%! assert_equal (corr ([1 2 3]), NaN (3));
%! assert_equal (corr (5), NaN);

## Test input validation
%!error <Invalid call to corr> corr ()
%!error <corr: X must be a numeric or logical matrix.> corr ({1, 2, 3})
%!error <corr: X must be a numeric or logical matrix.> corr ('abcde')
%!error <corr: X must be a matrix or a vector.> corr (ones (2, 2, 2))
%!error <corr: X must have at least one column.> corr (zeros (5, 0))
%!error <corr: X must have at least one column.> corr ([])
%!error <corr: Y must be a numeric or logical matrix.> corr ([1; 2], {1, 2})
%!error <corr: Y must be a matrix or a vector.> corr ([1; 2], ones (2, 2, 2))
%!error <corr: X and Y must have the same number of rows.> ...
%! corr (ones (5, 2), ones (4, 2))
%!error <corr: optional arguments must be 'Type', 'Rows', 'Tail' or 'Weights' Name-Value pairs.> ...
%! corr (ones (5, 2), 'Foo', 1)
%!error <corr: optional arguments must be 'Type', 'Rows', 'Tail' or 'Weights' Name-Value pairs.> ...
%! corr (ones (5, 2), 'Type')
%!error <corr: 'Type' must be 'Pearson', 'Kendall' or 'Spearman'.> ...
%! corr (ones (5, 2), 'Type', 'foo')
%!error <corr: 'Type' must be 'Pearson', 'Kendall' or 'Spearman'.> ...
%! corr (ones (5, 2), 'Type', 5)
%!error <corr: 'Rows' must be 'all', 'complete' or 'pairwise'.> ...
%! corr (ones (5, 2), 'Rows', 'foo')
%!error <corr: 'Tail' must be 'both', 'right' or 'left'.> ...
%! corr (ones (5, 2), 'Tail', 'foo')
%!error <corr: 'Weights' must be a column vector of 5 real values.> ...
%! corr (ones (5, 2), 'Weights', [1; 2])
%!error <corr: 'Weights' must be a column vector of 5 real values.> ...
%! corr (ones (5, 2), 'Weights', ones (1, 5))
%!error <corr: 'Weights' must be nonnegative.> ...
%! corr (ones (5, 2), 'Weights', [1; -1; 1; 1; 1])
