## Copyright (C) 1996-2017 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{h} =} ztest2 (@var{x1}, @var{n1}, @var{x2}, @var{n2})
## @deftypefnx {statistics} {@var{h} =} ztest2 (@var{x1}, @var{n1}, @var{x2}, @var{n2}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{h}, @var{p}] =} ztest2 (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{ci}] =} ztest2 (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{ci}, @var{stats}] =} ztest2 (@dots{})
##
## Two proportions Z-test.
##
## If @var{x1} and @var{n1} are the counts of successes and trials in one
## sample, and @var{x2} and @var{n2} those in a second one, test the null
## hypothesis that the success probabilities @math{p1} and @math{p2} are the
## same.  The result is @var{h} = 0 if the null hypothesis cannot be rejected at
## the significance level @var{alpha}, or @var{h} = 1 if it can.
##
## Under the null, the test statistic approximately follows a standard normal
## distribution.
##
## The size of @var{h} and @var{p} is the common size of @var{x1}, @var{n1},
## @var{x2}, and @var{n2}, which must be scalars or of common size.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## @code{[@var{h}, @var{p}] = ztest2 (@dots{})} returns the p-value.  That
## is the probability of observing the given result, or one more extreme, by
## chance if the null hypothesis true.
##
## @code{[@var{h}, @var{p}, @var{ci}] = ztest2 (@dots{})} returns a
## @math{100 (1 - alpha)}% confidence interval for @math{p1 - p2}, Newcombe's
## hybrid score interval, which combines the Wilson score intervals of the two
## proportions.  It has one row per test, and is one-sided for a one-sided
## test: @code{[lower, 1]} for @qcode{'right'}, @code{[-1, upper]} for
## @qcode{'left'}.
##
## @code{[@var{h}, @var{p}, @var{ci}, @var{stats}] = ztest2 (@dots{})} also
## returns a structure with the following fields:
##
## @multitable @columnfractions 0.2 0.75
## @item @qcode{zval} @tab the value of the test statistic.
## @item @qcode{CohensH} @tab Cohen's @math{h}, the difference between the
## arcsine transforms @math{2 asin (sqrt (p1)) - 2 asin (sqrt (p2))} of the
## two estimated proportions.
## @item @qcode{CohensHCI} @tab its confidence interval,
## @math{h -/+ z sqrt (1/n1 + 1/n2)} with @math{z} the normal quantile at the
## level of the test, one row per test, within @math{[-pi, pi]} and one-sided
## as @var{ci} is.
## @end multitable
##
## @code{[@dots{}] = ztest2 (@dots{}, @var{Name}, @var{Value}, @dots{})}
## specifies one or more of the following @var{Name}/@var{Value} pairs:
##
## @multitable @columnfractions 0.2 0.75
## @headitem @var{Name} @tab @var{Value}
## @item @qcode{'alpha'} @tab the significance level of @var{h} and the level
## of both confidence intervals. Default is 0.05.
##
## @item @qcode{'tail'} @tab a string specifying the alternative hypothesis
## @end multitable
## @multitable @columnfractions 0.25 0.65
## @item @qcode{'both'} @tab @math{p1} is not @math{p2}
## (two-tailed, default)
## @item @qcode{'left'} @tab @math{p1} is less than @math{p2}
## (left-tailed)
## @item @qcode{'right'} @tab @math{p1} is greater than @math{p2}
## (right-tailed)
## @end multitable
##
## @seealso{chi2test, fishertest}
## @end deftypefn

function [h, p, ci, stats] = ztest2 (x1, n1, x2, n2, varargin)

  if (nargin < 4)
    print_usage ();
  endif

  if (! isscalar (x1) || ! isscalar (n1) || ! isscalar (x2) || ! isscalar (n2))
    [retval, x1, n1, x2, n2] = common_size (x1, n1, x2, n2);
    if (retval > 0)
      error ("ztest2: X1, N1, X2, and N2 must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x1) || iscomplex (n1) || iscomplex (x2) || iscomplex (n2))
    error ("ztest2: X1, N1, X2, and N2 must not be complex.");
  endif

  if (any (x1(:) > n1(:)) || any (x2(:) > n2(:)))
    error ("ztest2: X1 must be <= N1 and X2 must be <= N2.");
  endif

  ## Add defaults and parse optional arguments
  alpha = 0.05;
  tail = 'both';
  if (nargin > 4)
    params = numel (varargin);
    if ((params / 2) != fix (params / 2))
      error ("ztest2: optional arguments must be in NAME-VALUE pairs.")
    endif
    for idx = 1:2:params
      name = varargin{idx};
      value = varargin{idx+1};
      switch (lower (name))
        case 'alpha'
          alpha = value;
          if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                alpha <= 0 || alpha >= 1)
            error ("ztest2: invalid VALUE for alpha.");
          endif
        case 'tail'
          tail = lower (value);
          if (! any (strcmpi (tail, {'both', 'left', 'right'})))
            error ("ztest2: invalid VALUE for tail.");
          endif
        otherwise
          error ("ztest2: invalid NAME for optional arguments.");
      endswitch
    endfor
  endif

  p1 = x1 ./ n1;
  p2 = x2 ./ n2;
  pc = (x1 + x2) ./ (n1 + n2);

  zvalue  = (p1 - p2) ./ sqrt (pc .* (1 - pc) .* (1 ./ n1 + 1 ./ n2));

  switch (tail)
    case 'both'
      p = 2 * normcdf (-abs (zvalue));
      z = norminv (1 - alpha / 2);
    case 'right'
      p = normcdf (-zvalue);
      z = norminv (1 - alpha);
    case 'left'
      p = normcdf (zvalue);
      z = norminv (1 - alpha);
  endswitch

  ## Determine the test outcome
  h = double (p < alpha);
  h(isnan (p)) = NaN;

  ## Newcombe's hybrid score interval for p1 - p2, from the Wilson score
  ## interval of each proportion
  [l1, u1] = wilson (x1(:), n1(:), z);
  [l2, u2] = wilson (x2(:), n2(:), z);
  d = p1(:) - p2(:);
  lo = d - sqrt ((p1(:) - l1) .^ 2 + (u2 - p2(:)) .^ 2);
  hi = d + sqrt ((u1 - p1(:)) .^ 2 + (p2(:) - l2) .^ 2);
  ci = [lo, hi];

  ## Cohen's h, with its interval on the arcsine scale
  stats.zval = zvalue;
  stats.CohensH = 2 * asin (sqrt (p1)) - 2 * asin (sqrt (p2));
  se = sqrt (1 ./ n1(:) + 1 ./ n2(:));
  hCI = stats.CohensH(:) + [-1, 1] .* z .* se;
  hCI = min (max (hCI, -pi), pi);

  ## A one-sided test has a one-sided interval
  if (strcmp (tail, 'right'))
    ci(:,2) = 1;
    hCI(:,2) = pi;
  elseif (strcmp (tail, 'left'))
    ci(:,1) = -1;
    hCI(:,1) = -pi;
  endif
  stats.CohensHCI = hCI;

endfunction

## The Wilson score interval of the proportion X / N at the normal quantile Z
function [lo, hi] = wilson (x, n, z)

  phat = x ./ n;
  centre = (x + z ^ 2 / 2) ./ (n + z ^ 2);
  half = z * sqrt (n) ./ (n + z ^ 2) ...
         .* sqrt (phat .* (1 - phat) + z ^ 2 ./ (4 * n));
  lo = centre - half;
  hi = centre + half;

endfunction

%!test
%! ## Values from R's prop.test without continuity correction
%! [h, p] = ztest2 (30, 100, 18, 100);
%! assert_equal ([h, p], [1, 0.046944726978481718], -1e-12);
%!test
%! [~, ~, ~, st] = ztest2 (30, 100, 18, 100);
%! assert_equal (st.zval ^ 2, 3.9473684210526314, -1e-12);
%!test
%! ## Newcombe's interval, from the Wilson intervals of R's prop.test
%! [~, ~, ci] = ztest2 (30, 100, 18, 100);
%! assert_equal (ci, [0.0013340883914530755, 0.23469812426823411], -1e-12);
%!test
%! [h, p] = ztest2 (45, 60, 20, 70, 'tail', 'right');
%! assert_equal ([h, p], [1, 6.5305496467962847e-08], -1e-12);
%!test
%! [~, ~, ci] = ztest2 (45, 60, 20, 70, 'tail', 'right');
%! assert_equal (ci, [0.32502269598878675, 1], -1e-12);
%!test
%! [~, ~, ci] = ztest2 (45, 60, 20, 70, 'tail', 'left');
%! assert_equal (ci(1), -1);
%!test
%! [h, p] = ztest2 (30, 100, 18, 100, 'alpha', 0.01);
%! assert_equal (h, 0);
%!test
%! ## Below the resolution of 1 - normcdf
%! [~, p] = ztest2 (900, 1000, 300, 1000);
%! assert_equal (p, 4.01237554141706e-165, -1e-10);
%!test
%! [~, ~, ~, st] = ztest2 (30, 100, 18, 100);
%! assert_equal (st.CohensH, 2 * asin (sqrt (0.3)) - 2 * asin (sqrt (0.18)), ...
%!               -1e-14);
%!test
%! [~, ~, ~, st] = ztest2 (30, 100, 18, 100);
%! assert_equal (st.CohensHCI, ...
%!               st.CohensH + [-1, 1] * norminv (0.975) * sqrt (2 / 100), ...
%!               -1e-14);
%!test
%! [~, ~, ~, st] = ztest2 (45, 60, 20, 70, 'tail', 'right');
%! assert_equal (st.CohensHCI(2), pi);
%!test
%! ## One row of the interval per test
%! [~, ~, ci, st] = ztest2 ([30; 45], [100; 60], [18; 20], [100; 70]);
%! assert_equal (size (ci), [2, 2]);
%!test
%! [~, ~, ci] = ztest2 ([30; 45], [100; 60], [18; 20], [100; 70]);
%! [~, ~, c1] = ztest2 (45, 60, 20, 70);
%! assert_equal (ci(2,:), c1);

## Test input validation
%!error ztest2 ();
%!error ztest2 (1);
%!error ztest2 (1, 2);
%!error ztest2 (1, 2, 3);
%!error ztest2 (1, 2, 3, 2);
%!error<ztest2: optional arguments must be in NAME-VALUE pairs.> ...
%! ztest2 (1, 2, 3, 4, 'alpha')
%!error<ztest2: invalid VALUE for alpha.> ...
%! ztest2 (1, 2, 3, 4, 'alpha', 0);
%!error<ztest2: invalid VALUE for alpha.> ...
%! ztest2 (1, 2, 3, 4, 'alpha', 1.2);
%!error<ztest2: invalid VALUE for alpha.> ...
%! ztest2 (1, 2, 3, 4, 'alpha', 'val');
%!error<ztest2: invalid VALUE for tail.>  ...
%! ztest2 (1, 2, 3, 4, 'tail', 'val');
%!error<ztest2: invalid VALUE for tail.>  ...
%! ztest2 (1, 2, 3, 4, 'alpha', 0.01, 'tail', 'val');
%!error<ztest2: invalid NAME for optional arguments.> ...
%! ztest2 (1, 2, 3, 4, 'alpha', 0.01, 'tail', 'both', 'badoption', 3);
