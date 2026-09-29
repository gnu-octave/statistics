## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{h} =} chi2test (@var{x})
## @deftypefnx {statistics} {[@var{h}, @var{p}] =} chi2test (@var{x})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{stats}] =} chi2test (@var{x})
## @deftypefnx {statistics} {[@dots{}] =} chi2test (@var{x}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@dots{}] =} chi2test (@dots{}, @qcode{'Alpha'}, @var{alpha})
##
## Perform a chi-squared test (for independence or homogeneity).
##
## For 2-way contingency tables, @code{chi2test} performs a chi-squared test
## for independence or homogeneity, according to the sampling scheme and related
## question.  Independence means that the two variables forming the 2-way
## table are not associated, hence you cannot predict from one another.
## Homogeneity refers to the concept of similarity, hence they all come from the
## same distribution.
##
## Both tests are computationally identical and will produce the same result.
## Nevertheless, they answer to different questions.  Consider two variables,
## one for gender and another for smoking.  To test independence (whether gender
## and smoking is associated), we would randomly sample from the general
## population and break them down into categories in the table.  To test
## homogeneity (whether men and women share the same smoking habits), we would
## sample individuals from within each gender, and then measure their smoking
## habits (e.g. smokers vs non-smokers).
##
## When @code{chi2test} is called without any output arguments, it will print
## the result in the terminal including p-value, chi^2 statistic, and degrees of
## freedom.  Otherwise it returns @var{h}, 1 if the null hypothesis is rejected
## at the significance level @var{alpha} and 0 otherwise, the p-value @var{p},
## and a structure @var{stats} with the following fields:
##
## @multitable @columnfractions 0.2 0.75
## @item @qcode{chi2stat} @tab the chi^2 statistic of the test.
## @item @qcode{df} @tab its degrees of freedom.
## @item @qcode{O} @tab the table tested: @var{x}, or for @qcode{"marginal"} the
## two-way table left after collapsing @var{x}.
## @item @qcode{E} @tab the expected values under the tested model, in the
## layout of @qcode{O}.
## @item @qcode{CohensW} @tab Cohen's @math{w}, the effect size of the test.
## @item @qcode{CohensWCI} @tab its confidence interval.
## @item @qcode{CramersV} @tab Cramer's @math{V}, for a two-way table and for
## the @qcode{"joint"} and @qcode{"marginal"} models, which test one.
## @item @qcode{CramersVCI} @tab its confidence interval.
## @end multitable
##
## Both effect sizes are corrected for bias.  With @math{n} observations,
## @math{w = sqrt (max (chi^2 - df, 0) / n)}, and @math{V = w / sqrt (k)} with
## @math{k} one less than the smaller dimension of the two-way table; for
## @qcode{"joint"} that table has the variable it names as rows and every
## combination of the other two as columns.  The uncorrected
## @math{sqrt (chi^2 / n)} is biased upwards, by @math{df / n} in its square.
## The confidence intervals invert the noncentral chi^2 distribution: the
## bounds @math{lambda} of the noncentrality consistent with the observed
## chi^2 give @math{sqrt (lambda / n)} for @math{w}, divided by
## @math{sqrt (k)} for @math{V}, whose upper bound is at most 1.
##
## @qcode{'Alpha'} sets the significance level @var{alpha} of @var{h} and the
## level @math{100 (1 - alpha)}% of both confidence intervals.  It is 0.05 by
## default.
##
## Unlike MATLAB, in GNU Octave @code{chi2test} also supports 3-way tables,
## which involve three categorical variables (each in a different dimension of
## @var{x}.  In its simplest form, @code{[@dots{}] = chi2test (@var{x})} will
## test for mutual independence among the three variables.  Alternatively,
## when called in the form @code{[@dots{}] = chi2test (@var{x}, @var{name},
## @var{value})}, it can perform the following tests:
##
## @multitable @columnfractions 0.2 0.1 0.7
## @headitem @var{name} @tab @var{value} @tab Description
## @item "mutual" @tab [] @tab Mutual independence.  All variables are
## independent from each other, (A, B, C).  It takes no value; an empty
## matrix is accepted.
## @item "joint" @tab scalar @tab Joint independence.  Two variables are jointly
## independent of the third, (AB, C). The scalar value corresponds to the
## dimension of the independent variable (i.e. 3 for C).
## @item "marginal" @tab scalar @tab Marginal independence.  Two variables are
## independent if you ignore the third, (A, C).  The scalar value corresponds
## to the dimension of the variable to be ignored (i.e. 2 for B); the table is
## collapsed over it and independence is tested in the two-way table that
## remains.
## @item "conditional" @tab scalar @tab Conditional independence.  Two variables
## are independent given the third, (AC, BC).  The scalar value corresponds to
## the dimension of the variable that forms the conditional dependence
## (i.e. 3 for C).
## @item "homogeneous" @tab [] @tab Homogeneous associations.  Conditional
## (partial) odds-ratios are not related on the value of the third,
## (AB, AC, BC).  It takes no value; an empty matrix is accepted.
## @end multitable
##
## When testing for homogeneous associations in 3-way tables, the iterative
## proportional fitting procedure is used.  For small samples it is better to
## use the Cochran-Mantel-Haenszel Test.  K-way tables for k > 3 are supported
## only for testing mutual independence.  Similar to 2-way tables, no optional
## parameters are required for k > 3 multi-way tables.
##
## @code{chi2test} produces a warning if any cell of a 2x2 table has an expected
## frequency less than 5 or if more than 20% of the cells in larger 2-way tables
## have expected frequencies less than 5 or any cell with expected frequency
## less than 1.  In such cases, use @code{fishertest}.
##
## @seealso{crosstab, fishertest, mcnemar_test}
## @end deftypefn

function [h, p, stats] = chi2test (x, varargin)
  ## Check input arguments
  if (nargin < 1)
    print_usage ();
  endif
  if (isvector (x))
    error ("chi2test: X must be a matrix.");
  endif
  if (! isreal (x))
    error ("chi2test: values in X must be real numbers.");
  endif
  if (any (isnan (x(:))))
    error ("chi2test: X must not have missing values (NaN).");
  endif
  ## Get size and dimensions of contingency table
  sz = size (x);
  dim = length (sz);
  ## Parse the optional arguments: a model for a 3-way table with the
  ## dimension it names ('mutual' and 'homogeneous' take none, though an
  ## empty value is accepted), and 'Alpha'
  alpha = 0.05;
  model = 'mutual';
  c_dim = [];
  models = {'mutual', 'joint', 'marginal', 'conditional', 'homogeneous'};
  i = 1;
  while (i <= numel (varargin))
    name = varargin{i};
    if (ischar (name) && strcmpi (name, 'alpha'))
      if (i == numel (varargin))
        error ("chi2test: optional arguments must be in pairs.");
      endif
      alpha = varargin{i+1};
      if (! (isnumeric (alpha) && isscalar (alpha) && isreal (alpha)
             && alpha > 0 && alpha < 1))
        error ("chi2test: invalid value for alpha.");
      endif
      i += 2;
    elseif (ischar (name) && any (strcmpi (name, models)))
      if (dim != 3)
        error ("chi2test: a model applies only to 3-way tables.");
      endif
      model = lower (name);
      novalue = any (strcmp (model, {'mutual', 'homogeneous'}));
      if (novalue && (i == numel (varargin) || ischar (varargin{i+1})))
        i += 1;
        continue;
      endif
      if (i == numel (varargin))
        error ("chi2test: optional arguments must be in pairs.");
      endif
      value = varargin{i+1};
      if (! isnumeric (value))
        error (strcat ("chi2test: value must be numeric in optional", ...
                       " argument name/value pair, for 3-way tables."));
      endif
      if (numel (value) > 1)
        error (strcat ("chi2test: value must be empty or scalar in", ...
                       " optional argument name/value pair, for 3-way", ...
                       " tables."));
      endif
      if (! novalue && ! (isscalar (value) && any (value == 1:3)))
        error ("chi2test: the dimension must be 1, 2, or 3.");
      endif
      c_dim = value;
      i += 2;
    elseif (dim == 3)
      error ("chi2test: invalid model name for testing a 3-way table.");
    else
      error ("chi2test: invalid optional argument.");
    endif
  endwhile

  ## The smaller dimension of a two-way table less one, for Cramer's V;
  ## empty where the test has no two-way table
  vk = [];
  ## Calculate total sample size
  n = sum (x(:));
  ## For 2-way contingency table
  if (length (sz) == 2)
    ## Calculate degrees of freedom
    df = prod (sz - 1);
    vk = min (sz) - 1;
    ## Calculate expected values
    E = sum (x')' * sum (x) / n;
  ## For 3-way contingency table
  elseif (length (sz) == 3)
    if (strcmp (model, 'mutual'))
      ## Calculate degrees of freedom
      df = prod (sz) - sum (sz) + 2;
      ## Calculate marginal table sums
      q1 = sum (sum (x, 2), 3);
      q2 = sum (sum (x, 1), 3);
      q3 = sum (sum (x, 1), 2);
      ns = sum (x(:)) ^ 2;
      for d1 = 1:size (x, 1)
        for d2 = 1:size (x, 2)
          for d3 = 1:size (x, 3)
            E(d1,d2,d3) = q1(d1,:,:) * q2(:,d2,:) * q3(:,:,d3) / ns;
          endfor
        endfor
      endfor
    elseif (strcmp (model, 'joint'))
      ## The test is of the table with the named variable as rows and every
      ## combination of the other two as columns
      vk = min (sz(c_dim), prod (sz) / sz(c_dim)) - 1;
      ## Calculate degrees of freedom
      c_sz = sz;
      c_sz(c_dim) = [];
      df = (sz(c_dim) - 1) * (prod (c_sz) - 1);
      ## Rearrange dimensions so that independent variable goes in dim 1
      dm = [1, 2, 3];
      dm(c_dim) = [];
      x = permute (x, [c_dim, dm]);
      ## Calculate partial table sums
      q1 = sum (sum (x, 1), 1);
      q2 = sum (sum (x, 2), 3);
      n = sum (x(:));
      for d1 = 1:size (x, 1)
        for d2 = 1:size (x, 2)
          for d3 = 1:size (x, 3)
            E(d1,d2,d3) = q1(:,d2,d3) * q2(d1) / n;
          endfor
        endfor
      endfor
      ## Rearrange OBSERVED and EXPECTED matrices in original dimensions
      x = ipermute (x, [c_dim, dm]);
      E = ipermute (E, [c_dim, dm]);
    elseif (strcmp (model, 'marginal'))
      ## Collapse over the ignored variable and test independence in the
      ## two-way table that remains
      x = reshape (sum (x, c_dim), sz(setdiff (1:3, c_dim)));
      sz = size (x);
      vk = min (sz) - 1;
      dim = 2;
      df = prod (sz - 1);
      E = sum (x, 2) * sum (x, 1) / n;
    elseif (strcmp (model, 'conditional'))
      ## Calculate degrees of freedom
      c_sz = sz;
      c_sz(c_dim) = [];
      df = prod (c_sz - 1) * sz(c_dim);
      ## Rearrange dimensions so that conditional variable goes in dim 1
      dm = [1, 2, 3];
      dm(c_dim) = [];
      x = permute (x, [c_dim, dm]);
      ## Calculate partial table sums
      q1 = sum (sum (x, 3), 3);
      q2 = sum (sum (x, 2), 2);
      q3 = sum (sum (x, 2), 3);
      ## Calculate expected values
      for d1 = 1:size (x, 1)
        for d2 = 1:size (x, 2)
          for d3 = 1:size (x, 3)
            E(d1,d2,d3) = q1(d1,d2) * q2(d1,:,d3) / q3(d1);
          endfor
        endfor
      endfor
      ## Rearrange OBSERVED and EXPECTED matrices in original dimensions
      x = ipermute (x, [c_dim, dm]);
      E = ipermute (E, [c_dim, dm]);
    else  # homogeneous
      ## Calculate degrees of freedom
      df = prod (sz - 1);
      ## Fit the three two-way margins by iterative proportional fitting,
      ## until those it last left behind match to within rounding
      ratio = @(a, b) a ./ (b + (b == 0));
      E = ones (sz);
      converged = false;
      for iter = 1:1000
        E = E .* ratio (sum (x, 3), sum (E, 3));
        E = E .* ratio (sum (x, 2), sum (E, 2));
        E = E .* ratio (sum (x, 1), sum (E, 1));
        gap = max ([abs(sum (E, 3) - sum (x, 3))(:); ...
                    abs(sum (E, 2) - sum (x, 2))(:)]);
        if (gap <= 1e-12 * n)
          converged = true;
          break;
        endif
      endfor
      if (! converged)
        warning ("chi2test: the homogeneous model did not converge.");
      endif
    endif
  ## For k-way contingency table, where k > 3
  else
    ## Calculate degrees of freedom
    df = prod (sz) - sum (sz) + 2;
    ## Calculate squared sample size
    ns = sum (x(:)) ^ (dim - 1);
    ## Calculate marginal table sums for each available dimension
    for i = 1:dim
      qi(i) = {x};
      remdim = [1:dim];
      remdim(remdim == i) = [];
      for j = 1:length (remdim)
        qi(i) = sum (qi{i}, remdim(j));
      endfor
      qi(i) = squeeze (qi{i});
    endfor
    ## Iterate through all cells
    cn = numel (x);
    for i = 1:cn
      E(i) = 1;
      cid = i;
      ## Keep track of indexing
      for d = dim - 1:-1:1
        idx(d+1) = ceil (cid / prod (sz(1:d)));
        if (idx(d+1) > 1)
          cid -= (idx(d+1) - 1) * prod (sz(1:d));
        endif
      endfor
      idx(1) = cid;
      ## Calculate the expected value
      for j = 1:dim
        E(i) = E(i) * qi{j}(idx(j));
      endfor
      E(i) = E(i) / ns;
    endfor
    ## Reshape to original dimensions
    E = reshape (E, sz);
  endif
  ## Check expected values and display warnings
  if ((dim == 2 && isequal (sz, [2, 2]) && any (E(:) < 5)) || ...
      (dim == 2 && any (sz > 2) && sum (E(:) < 5) > 0.2 * numel (E)) || ...
      (dim > 2 && sum (E(:) < 5) > 0.2 * numel (E)))
    warning ("chi2test: Expected values less than 5.");
  endif
  if (any (E(:) < 1))
    warning ("chi2test: Expected values less than 1.");
  endif
  ## Calculate chi-squared and p-value
  cells = ((x - E) .^2) ./ E;
  chisq = sum (cells(:));
  p = chi2cdf (chisq, df, 'upper');
  h = double (p < alpha);

  ## Statistics, and the effect sizes with their confidence intervals
  stats.chi2stat = chisq;
  stats.df = df;
  stats.O = x;
  stats.E = E;
  if (isfinite (chisq))
    n = sum (x(:));
    stats.CohensW = sqrt (max (chisq - df, 0) / n);
    stats.CohensWCI = sqrt (ncx2bounds (chisq, df, alpha) / n);
  else
    stats.CohensW = NaN;
    stats.CohensWCI = [NaN, NaN];
  endif
  if (! isempty (vk) && vk > 0)
    stats.CramersV = stats.CohensW / sqrt (vk);
    stats.CramersVCI = min (stats.CohensWCI / sqrt (vk), 1);
  endif

  ## Print results if no output requested
  if (nargout == 0)
    printf ("p-val = %f with chi^2 statistic = %f and d.f. = %d.\n", ...
            p, chisq, df);
  endif

endfunction

## The noncentralities at which the observed statistic falls at the upper and
## at the lower ALPHA/2 point of the noncentral chi^2 distribution; zero
## where the central distribution already puts it beyond that point.
function lambda = ncx2bounds (chisq, df, alpha)

  target = [1 - alpha / 2, alpha / 2];
  lambda = [0, 0];
  for j = 1:2
    if (chi2cdf (chisq, df) > target(j))
      hi = max (chisq, 1);
      while (ncx2cdf (chisq, df, hi) > target(j))
        hi *= 2;
      endwhile
      lambda(j) = fzero (@(l) ncx2cdf (chisq, df, l) - target(j), [0, hi]);
    endif
  endfor

endfunction

%!function p = chi2p (varargin)
%!  [~, p] = chi2test (varargin{:});
%!endfunction

%!test
%! ## Below the resolution of 1 - chi2cdf
%! [~, p, st] = chi2test ([300, 100; 100, 300]);
%! assert_equal (p, erfc (sqrt (st.chi2stat / 2)), -1e-12);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [~, p] = chi2test (x);
%! assert_equal (p, 0.017787, 1e-6);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [~, ~, st] = chi2test (x);
%! assert_equal (st.chi2stat, 11.9421, 1e-4);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [~, ~, st] = chi2test (x);
%! assert_equal (st.df, 4);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [~, ~, st] = chi2test (x);
%! assert_equal (st.O, x);
%!assert_equal (chi2test ([11, 3, 8; 2, 9, 14; 12, 13, 28]), 1)
%!assert_equal (chi2test ([11, 3, 8; 2, 9, 14; 12, 13, 28], 'Alpha', 0.01), 0)
%!test
%! ## Effect sizes corrected for bias; the intervals are R's inversion of
%! ## pchisq with ncp by uniroot
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [~, ~, st] = chi2test (x);
%! assert_equal (st.CohensW, sqrt ((st.chi2stat - 4) / 100), -1e-14);
%!test
%! [~, ~, st] = chi2test ([11, 3, 8; 2, 9, 14; 12, 13, 28]);
%! assert_equal (st.CohensWCI, [0.054364192439641336, 0.50550722873705722], ...
%!               -1e-10);
%!test
%! [~, ~, st] = chi2test ([11, 3, 8; 2, 9, 14; 12, 13, 28]);
%! assert_equal (st.CramersV, 0.19927484317340177, -1e-12);
%!test
%! [~, ~, st] = chi2test ([11, 3, 8; 2, 9, 14; 12, 13, 28]);
%! assert_equal (st.CramersVCI, ...
%!               [0.054364192439641336, 0.50550722873705722] / sqrt (2), ...
%!               -1e-10);
%!test
%! ## A higher confidence level widens the interval on both sides
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [~, ~, s95] = chi2test (x);
%! [~, ~, s99] = chi2test (x, 'Alpha', 0.01);
%! assert_equal ([s99.CohensWCI(1) < s95.CohensWCI(1), ...
%!                s99.CohensWCI(2) > s95.CohensWCI(2)], [true, true]);
%!test
%!shared x
%! x(:,:,1) = [59, 32; 9,16];
%! x(:,:,2) = [55, 24;12,33];
%! x(:,:,3) = [107,80;17,56];
%!assert_equal (chi2p (x), 2.282063427117009e-11, 1e-14);
%!assert_equal (chi2p (x, 'mutual', []), 2.282063427117009e-11, 1e-14);
%!assert_equal (chi2p (x, 'joint', 1), 1.164834895206468e-11, 1e-14);
%!assert_equal (chi2p (x, 'joint', 2), 7.771350230001417e-11, 1e-14);
%!assert_equal (chi2p (x, 'joint', 3), 0.07151361728026107, 1e-14);
%!assert_equal (chi2p (x, 'marginal', 1), 0.12455768155123595, -1e-12);
%!assert_equal (chi2p (x, 'marginal', 2), 0.039793350279010681, -1e-12);
%!assert_equal (chi2p (x, 'marginal', 3), 9.0141038839122684e-13, -1e-12);
%!assert_equal (chi2p (x, 'conditional', 1), 0.2303114201312508, 1e-14);
%!assert_equal (chi2p (x, 'conditional', 2), 0.0958810684407079, 1e-14);
%!assert_equal (chi2p (x, 'conditional', 3), 2.648037344954446e-11, 1e-14);
%!assert_equal (chi2p (x, 'homogeneous', []), 0.57357539370887889, -1e-10);
%!assert_equal (chi2p (x, 'homogeneous'), chi2p (x, 'homogeneous', []));
%!assert_equal (chi2p (x, 'mutual'), chi2p (x));
%!assert_equal (chi2p (x, 'joint', 3, 'Alpha', 0.01), chi2p (x, 'joint', 3));
%!test
%! [~, ~, st] = chi2test (x);
%! assert_equal (st.chi2stat, 64.0982, 1e-4);
%! assert_equal (st.df, 7);
%! assert_equal (st.E(:,:,1), [42.903, 39.921; 17.185, 15.991], 1e-3);
%!test
%! [~, ~, st] = chi2test (x, 'joint', 2);
%! assert_equal (st.chi2stat, 56.0943, 1e-4);
%! assert_equal (st.df, 5);
%! assert_equal (st.E(:,:,2), [40.922, 38.078; 23.310, 21.690], 1e-3);
%!test
%! [~, ~, st] = chi2test (x, 'marginal', 3);
%! assert_equal (st.chi2stat, 51.0479, 1e-4);
%! assert_equal (st.df, 1);
%! assert_equal (st.O, sum (x, 3));
%! assert_equal (st.E, [184.926, 172.074; 74.074, 68.926], 1e-3);
%!test
%! [~, ~, st] = chi2test (x, 'conditional', 3);
%! assert_equal (st.chi2stat, 52.2509, 1e-4);
%! assert_equal (st.df, 3);
%! assert_equal (st.E(:,:,1), [53.345, 37.655; 14.655, 10.345], 1e-3);
%!test
%! [~, ~, st] = chi2test (x, 'homogeneous', []);
%! assert_equal (st.chi2stat, 1.1117, 1e-4);
%! assert_equal (st.df, 2);
%! assert_equal (st.E(:,:,1), [60.469, 30.531; 7.531, 17.469], 1e-3);
%!test
%! ## The homogeneous model reproduces every two-way margin
%! [~, ~, st] = chi2test (x, 'homogeneous');
%! E = st.E;
%! assert_equal ([sum(E, 1)(:); sum(E, 2)(:); sum(E, 3)(:)], ...
%!               [sum(x, 1)(:); sum(x, 2)(:); sum(x, 3)(:)], -1e-10);
%!test
%! ## 'joint' measures V on the table flattened to its variable against the
%! ## other two; bounds from R as above
%! [h, ~, st] = chi2test (x, 'joint', 3, 'Alpha', 0.01);
%! assert_equal (st.CohensWCI, [0, 0.24137302257386881], -1e-10);
%!test
%! [~, ~, st] = chi2test (x, 'marginal', 3);
%! assert_equal (st.CramersV, 0.31637908166758733, -1e-12);
%!test
%! [~, ~, st] = chi2test (x, 'marginal', 3);
%! assert_equal (st.CohensWCI, [0.23187195991809992, 0.40717646803341606], ...
%!               -1e-10);
%!test
%! ## No two-way table, no Cramer's V
%! [~, ~, st] = chi2test (x, 'conditional', 3);
%! assert_equal (isfield (st, 'CramersV'), false);
%!test
%! [~, ~, st] = chi2test (x);
%! assert_equal (isfield (st, 'CramersV'), false);
%!test
%! ## E keeps the layout of a table whose dimensions differ
%! y = reshape ([12 7 3 9 15 4 6 8 11 5 14 2 9 10 3 7 6 13], [3, 2, 3]);
%! [~, ~, st] = chi2test (y, 'joint', 2);
%! assert_equal (sum (st.E, 2), sum (y, 2), -1e-12);
%!test
%! y = reshape ([12 7 3 9 15 4 6 8 11 5 14 2 9 10 3 7 6 13], [3, 2, 3]);
%! [~, ~, st] = chi2test (y, 'conditional', 2);
%! assert_equal (size (st.E), [3, 2, 3]);
%!test
%! ## Chi-squares of R's loglin on a 3-by-2-by-3 table
%! y = reshape ([12 7 3 9 15 4 6 8 11 5 14 2 9 10 3 7 6 13], [3, 2, 3]);
%! [~, ~, st] = chi2test (y, 'homogeneous');
%! assert_equal (st.chi2stat, 15.05772742, -1e-9);
%!test
%! ## Marginal independence is tested in the collapsed table, as R's
%! ## chisq.test does it
%! y = reshape ([12 7 3 9 15 4 6 8 11 5 14 2 9 10 3 7 6 13], [3, 2, 3]);
%! [~, ~, st] = chi2test (y, 'marginal', 2);
%! assert_equal ([st.chi2stat, st.df], [7.584460, 4], -1e-6);

## Check warnings
%!warning<chi2test: Expected values less than 5.> chi2test (ones (2));
%!warning<chi2test: Expected values less than 5.> chi2test (ones (3, 2));
%!warning<chi2test: Expected values less than 1.> chi2test (0.4 * ones (3));

## Test input validation
%!error chi2test ();
%!error<chi2test: X must be a matrix.> chi2test ([1, 2, 3, 4, 5]);
%!error<chi2test: values in X must be real numbers.> ...
%! chi2test ([1, 2; 2, 1+3i]);
%!error<chi2test: X must not have missing values \(NaN\).> ...
%! chi2test ([NaN, 6; 34, 12]);
%!error<chi2test: a model applies only to 3-way tables.> ...
%! chi2test (ones (3, 3), 'mutual', []);
%!error<chi2test: a model applies only to 3-way tables.> ...
%! chi2test (ones (3, 3, 3, 4), 'mutual');
%!error<chi2test: invalid model name for testing a 3-way table.> ...
%! chi2test (ones (3, 3, 3), 'testtype', 2);
%!error<chi2test: invalid optional argument.> ...
%! chi2test (ones (3, 3), 'testtype', 2);
%!error<chi2test: optional arguments must be in pairs.> ...
%! chi2test (ones (3, 3, 3), 'joint');
%!error<chi2test: optional arguments must be in pairs.> ...
%! chi2test (ones (3, 3), 'Alpha');
%!error<chi2test: value must be numeric in optional argument name/value pair, for 3-way tables.> ...
%! chi2test (ones (3, 3, 3), 'joint', 'a');
%!error<chi2test: value must be empty or scalar in optional argument name/value pair, for 3-way tables.> ...
%! chi2test (ones (3, 3, 3), 'joint', [2, 3]);
%!error<chi2test: the dimension must be 1, 2, or 3.> ...
%! chi2test (ones (3, 3, 3), 'joint', 4);
%!error<chi2test: the dimension must be 1, 2, or 3.> ...
%! chi2test (ones (3, 3, 3), 'marginal', []);
%!error<chi2test: invalid value for alpha.> chi2test (ones (3, 3), 'Alpha', 0);
%!error<chi2test: invalid value for alpha.> chi2test (ones (3, 3), 'Alpha', 1.5);
%!error<chi2test: invalid value for alpha.> ...
%! chi2test (ones (3, 3), 'Alpha', [0.1, 0.2]);
