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

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{h} =} ansaribradley (@var{x}, @var{y})
## @deftypefnx {statistics} {@var{h} =} ansaribradley (@var{x}, @var{y}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{h}, @var{p}] =} ansaribradley (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{stats}] =} ansaribradley (@dots{})
##
## Ansari-Bradley two-sample test for equal dispersions.
##
## @code{@var{h} = ansaribradley (@var{x}, @var{y})} performs an Ansari-Bradley
## test of the hypothesis that the two independent samples in the vectors
## @var{x} and @var{y} come from distributions with the same dispersion
## parameter, against the alternative that they come from distributions with
## different dispersions.  The result is @var{h} = 0 if the null hypothesis of
## equal dispersions cannot be rejected at the 5% significance level, or
## @var{h} = 1 if it can.
##
## The Ansari-Bradley test is a nonparametric alternative to the two-sample
## @math{F} test (@code{vartest2}) that does not assume normality.  It assumes
## that the two samples are independent and that they come from distributions
## with the same median and shape, differing (under the alternative) only in
## dispersion.  If the medians differ, the data should be recentred (e.g.@: by
## subtracting the sample medians) before applying the test.
##
## @code{ansaribradley} treats NaNs in @var{x} or @var{y} as missing values and
## ignores them.
##
## @code{[@var{h}, @var{p}] = ansaribradley (@dots{})} returns the p-value of
## the test, that is the probability, under the null hypothesis, of observing a
## value of the test statistic as or more extreme than the one observed.
##
## @code{[@var{h}, @var{p}, @var{stats}] = ansaribradley (@dots{})} returns a
## structure with the following fields:
##
## @multitable @columnfractions 0.2 0.75
## @item @qcode{W} @tab the value of the Ansari-Bradley test statistic, the sum
## of the Ansari-Bradley scores of the sample @var{x}
## @item @qcode{Wstar} @tab the value of the approximate normal (z) statistic
## @end multitable
##
## @code{[@dots{}] = ansaribradley (@dots{}, @var{name}, @var{value})} specifies
## one or more of the following name/value pairs:
##
## @multitable @columnfractions 0.2 0.75
## @headitem Name @tab Value
## @item @qcode{'alpha'} @tab the significance level.  Default is 0.05.
##
## @item @qcode{'tail'} @tab a string specifying the alternative hypothesis
##
## @item @qcode{'method'} @tab a string selecting the p-value computation,
## either @qcode{'exact'} to use the exact permutation distribution of the
## statistic, or @qcode{'approximate'} to use the normal approximation.  The
## default is @qcode{'exact'} when the total sample size is 25 or less, and
## @qcode{'approximate'} otherwise.
## @end multitable
##
## The @qcode{'tail'} option can take one of the following values:
##
## @multitable @columnfractions 0.15 0.75
## @item @qcode{'both'} @tab dispersions are not equal (two-tailed, default)
## @item @qcode{'right'} @tab dispersion of @var{x} is greater than dispersion
## of @var{y} (right-tailed)
## @item @qcode{'left'} @tab dispersion of @var{x} is less than dispersion of
## @var{y} (left-tailed)
## @end multitable
##
## @seealso{vartest2, vartestn, kstest2, ranksum}
## @end deftypefn

function [h, p, stats] = ansaribradley (x, y, varargin)

  ## Validate input arguments
  if (nargin < 2)
    error ("ansaribradley: too few input arguments.");
  endif
  if (! isvector (x) || ! isvector (y))
    error ("ansaribradley: X and Y must be vectors.");
  endif
  ## Remove missing data and make column vectors
  x = x(! isnan (x))(:);
  y = y(! isnan (y))(:);
  if (isempty (x))
    error ("ansaribradley: not enough data in X.");
  endif
  if (isempty (y))
    error ("ansaribradley: not enough data in Y.");
  endif

  ## Add defaults and parse optional name/value pairs
  alpha = 0.05;
  tail = 'both';
  method = [];
  if (mod (numel (varargin), 2) != 0)
    error ("ansaribradley: optional arguments must be in name/value pairs.");
  endif
  for idx = 1:2:numel (varargin)
    name = varargin{idx};
    value = varargin{idx + 1};
    switch (lower (name))
      case 'alpha'
        alpha = value;
        if (! isscalar (alpha) || ! isnumeric (alpha) || ! isreal (alpha) ...
            || alpha <= 0 || alpha >= 1)
          error ("ansaribradley: invalid value for alpha.");
        endif
      case 'tail'
        tail = value;
        if (! (ischar (tail) && isrow (tail)) ...
            || ! any (strcmpi (tail, {'both', 'left', 'right'})))
          error ("ansaribradley: invalid value for tail.");
        endif
      case 'method'
        method = value;
        if (! (ischar (method) && isrow (method)) ...
            || ! any (strcmpi (method, {'exact', 'approximate'})))
          error ("ansaribradley: invalid value for method.");
        endif
      otherwise
        error ("ansaribradley: invalid name for optional arguments.");
    endswitch
  endfor

  nx = numel (x);
  ny = numel (y);
  N = nx + ny;

  ## Select the default computation method
  if (isempty (method))
    if (N <= 25)
      method = 'exact';
    else
      method = 'approximate';
    endif
  endif

  ## Compute the Ansari-Bradley scores of the pooled sample.  The raw score of
  ## the observation ranked I among the N pooled values is min (I, N + 1 - I),
  ## so that the extremes get the smallest scores and the centre the largest.
  ## Tied observations receive the average of their raw scores (mid-ranks).
  z = [x; y];
  [zs, ord] = sort (z);
  pos = (1:N)';
  rawscore = min (pos, N + 1 - pos);
  [~, ~, grp] = unique (zs);
  meanscore = accumarray (grp, rawscore, [], @mean);
  sortedscore = meanscore(grp);
  abscore = zeros (N, 1);
  abscore(ord) = sortedscore;

  ## Ansari-Bradley statistic: sum of the scores belonging to sample X
  W = sum (abscore(1:nx));

  ## Normal approximation statistic (tie-corrected), always reported in STATS
  abar = mean (abscore);
  EW = nx * abar;
  VarW = (nx * ny) / (N * (N - 1)) * sum ((abscore - abar) .^ 2);
  Wstar = (W - EW) / sqrt (VarW);

  ## Compute the p-value.  Larger dispersion in X pushes its observations to the
  ## extremes, lowering its scores and hence W; thus the 'right' alternative
  ## (dispersion of X greater) corresponds to the lower tail of W.
  if (strcmpi (method, 'exact'))
    ## Exact permutation distribution of W by dynamic programming over the
    ## scores doubled to integers (mid-ranks are integers or half-integers).
    s2 = round (2 * abscore);
    S = sum (s2);
    ## counts(k+1,t+1) = number of size-k subsets of scores summing to t
    counts = zeros (nx + 1, S + 1);
    counts(1, 1) = 1;
    for i = 1:N
      si = s2(i);
      for k = min (i, nx):-1:1
        counts(k + 1, (si + 1):(S + 1)) += counts(k, 1:(S + 1 - si));
      endfor
    endfor
    dist = counts(nx + 1, :);
    total = sum (dist);
    w2 = round (2 * W);
    p_le = sum (dist(1:(w2 + 1))) / total;
    p_ge = sum (dist((w2 + 1):end)) / total;
    switch (lower (tail))
      case 'both'
        p = min (1, 2 * min (p_le, p_ge));
      case 'right'
        p = p_le;
      case 'left'
        p = p_ge;
    endswitch
  else
    switch (lower (tail))
      case 'both'
        p = 2 * normcdf (- abs (Wstar));
      case 'right'
        p = normcdf (Wstar);
      case 'left'
        p = normcdf (- Wstar);
    endswitch
  endif

  ## Determine the test outcome and assemble the STATS structure
  h = double (p <= alpha);
  if (nargout > 2)
    stats = struct ('W', W, 'Wstar', Wstar);
  endif

endfunction

%!demo
%! ## Test whether two samples have the same dispersion.  The second sample is
%! ## drawn with twice the standard deviation, so the null hypothesis of equal
%! ## dispersions should be rejected.
%! x = [42, 44, 38, 52, 48, 46, 40, 50];
%! y = [30, 62, 25, 70, 33, 58, 20, 65];
%! [h, p, stats] = ansaribradley (x, y)

## Test input validation
%!error<ansaribradley: too few input arguments.> ansaribradley (1);
%!error<ansaribradley: X and Y must be vectors.> ...
%! ansaribradley (ones (3, 2), ones (3, 1));
%!error<ansaribradley: X and Y must be vectors.> ...
%! ansaribradley (ones (3, 1), ones (2, 2));
%!error<ansaribradley: not enough data in X.> ansaribradley ([NaN, NaN], [1, 2]);
%!error<ansaribradley: not enough data in Y.> ansaribradley ([1, 2], [NaN, NaN]);
%!error<ansaribradley: optional arguments must be in name/value pairs.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'alpha');
%!error<ansaribradley: invalid value for alpha.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'alpha', 0);
%!error<ansaribradley: invalid value for alpha.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'alpha', 1);
%!error<ansaribradley: invalid value for alpha.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'alpha', -0.2);
%!error<ansaribradley: invalid value for alpha.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'alpha', [0.01, 0.05]);
%!error<ansaribradley: invalid value for alpha.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'alpha', 'x');
%!error<ansaribradley: invalid value for tail.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'tail', 'other');
%!error<ansaribradley: invalid value for tail.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'tail', 5);
%!error<ansaribradley: invalid value for method.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'method', 'other');
%!error<ansaribradley: invalid value for method.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'method', 5);
%!error<ansaribradley: invalid name for optional arguments.> ...
%! ansaribradley ([1, 2, 3], [4, 5, 6], 'name', 'value');

## Test results
%!test
%! ## A concentrated sample versus a dispersed one (default exact, N = 16).
%! x = [42, 44, 38, 52, 48, 46, 40, 50];
%! y = [30, 62, 25, 70, 33, 58, 20, 65];
%! [h, p, stats] = ansaribradley (x, y);
%! assert_equal (h, 1);
%! assert_equal (p, 0.000155400155400155, 1e-15);
%! assert_equal (stats.W, 52);
%! assert_equal (stats.Wstar, 3.38061701891407, 1e-13);
%!test
%! ## Same data, left-tailed exact test (dispersion of X less than Y).
%! x = [42, 44, 38, 52, 48, 46, 40, 50];
%! y = [30, 62, 25, 70, 33, 58, 20, 65];
%! [h, p] = ansaribradley (x, y, 'tail', 'left');
%! assert_equal (h, 1);
%! assert_equal (p, 7.77000777000777e-05, 1e-16);
%!test
%! ## Same data, normal approximation.
%! x = [42, 44, 38, 52, 48, 46, 40, 50];
%! y = [30, 62, 25, 70, 33, 58, 20, 65];
%! [h, p, stats] = ansaribradley (x, y, 'method', 'approximate');
%! assert_equal (h, 1);
%! assert_equal (p, 0.000723232716430194, 1e-15);
%! assert_equal (stats.Wstar, 3.38061701891407, 1e-13);
%!test
%! ## Total sample size 30 (> 25) defaults to the approximate method.
%! u = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
%! v = [3, 3, 4, 5, 6, 7, 8, 8, 9, 9, 10, 11, 12, 12, 13];
%! [h, p, stats] = ansaribradley (u, v);
%! assert_equal (h, 0);
%! assert_equal (p, 0.184251738662491, 1e-13);
%! assert_equal (stats.Wstar, -1.32777714715389, 1e-13);
%!test
%! ## Tied observations receive mid-rank scores (exact method).
%! x = [1, 1, 2, 3, 3];
%! y = [0, 2, 2, 2, 4, 4];
%! [h, p, stats] = ansaribradley (x, y);
%! assert_equal (h, 0);
%! assert_equal (p, 0.904761904761905, 1e-14);
%! assert_equal (stats.W, 17);
%! assert_equal (stats.Wstar, 0.245274554572897, 1e-14);
%!test
%! ## NaNs are ignored as missing values.
%! x = [1, 2, 3, NaN, 5];
%! y = [2, NaN, 4, 6];
%! [h, p, stats] = ansaribradley (x, y);
%! assert_equal (h, 0);
%! assert_equal (p, 0.971428571428571, 1e-14);
%! assert_equal (stats.W, 9.5);
%! assert_equal (stats.Wstar, 0.253836541283405, 1e-14);
