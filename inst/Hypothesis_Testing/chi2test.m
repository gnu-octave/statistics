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
## @deftypefn  {statistics} {@var{pval} =} chi2test (@var{x})
## @deftypefnx {statistics} {[@var{pval}, @var{chisq}] =} chi2test (@var{x})
## @deftypefnx {statistics} {[@var{pval}, @var{chisq}, @var{dF}] =} chi2test (@var{x})
## @deftypefnx {statistics} {[@var{pval}, @var{chisq}, @var{dF}, @var{E}] =} chi2test (@var{x})
## @deftypefnx {statistics} {[@dots{}] =} chi2test (@var{x}, @var{name}, @var{value})
##
## Perform a chi-squared test (for independence or homogeneity).
##
## For 2-way contingency tables, @code{chi2test} performs and a chi-squared test
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
## freedom.  Otherwise it can return the following output arguments:
##
## @multitable @columnfractions 0.1 0.85
## @item @var{pval} @tab the p-value of the relevant test.
## @item @var{chisq} @tab the chi^2 statistic of the relevant test.
## @item @var{dF} @tab the degrees of freedom of the relevant test.
## @item @var{E} @tab the expected values under the tested model, in the
## layout of @var{x}; for @qcode{"marginal"}, of the two-way table left after
## collapsing @var{x}.
## @end multitable
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

function [pval, chisq, df, E] = chi2test (x, varargin)
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
  ## Check optional arguments
  if (dim == 2 && nargin > 1)
    error ("chi2test: optional arguments are not supported for 2-way tables.");
  endif
  ## 'mutual' and 'homogeneous' take no value, though an empty one is accepted
  if (dim == 3 && numel (varargin) == 1 && ischar (varargin{1}) ...
      && any (strcmpi (varargin{1}, {'mutual', 'homogeneous'})))
    varargin{2} = [];
  endif
  if (dim == 3 && mod (numel (varargin(:)), 2) != 0)
    error ("chi2test: optional arguments must be in pairs.");
  endif
  if (dim == 3 && nargin > 1 && ! isnumeric (varargin{2}))
    error (strcat ("chi2test: value must be numeric in optional argument", ...
                   " name/value pair, for 3-way tables."));
  endif
  if (dim == 3 && nargin > 1 && numel (varargin{2}) > 1)
    error (strcat ("chi2test: value must be empty or scalar in optional", ...
                   " argument name/value pair, for 3-way tables."));
  endif
  if (dim >= 4 && nargin > 1)
    error ("chi2test: optional arguments are not supported for k>3.");
  endif
  ## Calculate total sample size
  n = sum (x(:));
  ## For 2-way contingency table
  if (length (sz) == 2)
    ## Calculate degrees of freedom
    df = prod (sz - 1);
    ## Calculate expected values
    E = sum (x')' * sum (x) / n;
  ## For 3-way contingency table
  elseif (length (sz) == 3)
    ## Check optional arguments
    if (nargin == 1 || strcmpi (varargin{1}, 'mutual'))
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
    elseif (strcmpi (varargin{1}, 'joint'))
      ## Get dimension of independent variable (dim)
      c_dim = varargin{2};
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
    elseif (strcmpi (varargin{1}, 'marginal'))
      ## Collapse over the ignored variable and test independence in the
      ## two-way table that remains
      c_dim = varargin{2};
      x = reshape (sum (x, c_dim), sz(setdiff (1:3, c_dim)));
      sz = size (x);
      dim = 2;
      df = prod (sz - 1);
      E = sum (x, 2) * sum (x, 1) / n;
    elseif (strcmpi (varargin{1}, 'conditional'))
      ## Get dimension of conditional variable (dim)
      c_dim = varargin{2};
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
    elseif (strcmpi (varargin{1}, 'homogeneous'))
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
    else
      error ("chi2test: invalid model name for testing a 3-way table.");
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
  pval = chi2cdf (chisq, df, 'upper');
  ## Print results if no output requested
  if (nargout == 0)
    printf ("p-val = %f with chi^2 statistic = %f and d.f. = %d.\n", ...
            pval, chisq, df);
  endif
endfunction

%!test
%! ## Below the resolution of 1 - chi2cdf
%! [p, c] = chi2test ([300, 100; 100, 300]);
%! assert_equal (p, erfc (sqrt (c / 2)), -1e-12);

## Input validation tests
%!error chi2test ();
%!error chi2test ([1, 2, 3, 4, 5]);
%!error chi2test ([1, 2; 2, 1+3i]);
%!error chi2test ([NaN, 6; 34, 12]);
%!error<chi2test: optional arguments are not supported for 2-way> ...
%! p = chi2test (ones (3, 3), 'mutual', []);
%!error<chi2test: invalid model name for testing a 3-way table.> ...
%! p = chi2test (ones (3, 3, 3), 'testtype', 2);
%!error<chi2test: optional arguments must be in pairs.> ...
%! p = chi2test (ones (3, 3, 3), 'joint');
%!error<chi2test: value must be numeric in optional argument> ...
%! p = chi2test (ones (3, 3, 3), 'joint', ['a']);
%!error<chi2test: value must be empty or scalar in optional argument> ...
%! p = chi2test (ones (3, 3, 3), 'joint', [2, 3]);
%!error<chi2test: optional arguments are not supported for k> ...
%! p = chi2test (ones (3, 3, 3, 4), 'mutual', [])

## Check warning
%!warning<chi2test: Expected values less than 5.> p = chi2test (ones (2));
%!warning<chi2test: Expected values less than 5.> p = chi2test (ones (3, 2));
%!warning<chi2test: Expected values less than 1.> p = chi2test (0.4 * ones (3));
## Output validation tests
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! p = chi2test (x);
%! assert_equal (p, 0.017787, 1e-6);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [p, chisq] = chi2test (x);
%! assert_equal (chisq, 11.9421, 1e-4);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [p, chisq, df] = chi2test (x);
%! assert_equal (df, 4);
%!test
%!shared x
%! x(:,:,1) = [59, 32; 9,16];
%! x(:,:,2) = [55, 24;12,33];
%! x(:,:,3) = [107,80;17,56];%!
%!assert_equal (chi2test (x), 2.282063427117009e-11, 1e-14);
%!assert_equal (chi2test (x, 'mutual', []), 2.282063427117009e-11, 1e-14);
%!assert_equal (chi2test (x, 'joint', 1), 1.164834895206468e-11, 1e-14);
%!assert_equal (chi2test (x, 'joint', 2), 7.771350230001417e-11, 1e-14);
%!assert_equal (chi2test (x, 'joint', 3), 0.07151361728026107, 1e-14);
%!assert_equal (chi2test (x, 'marginal', 1), 0.12455768155123595, -1e-12);
%!assert_equal (chi2test (x, 'marginal', 2), 0.039793350279010681, -1e-12);
%!assert_equal (chi2test (x, 'marginal', 3), 9.0141038839122684e-13, -1e-12);
%!assert_equal (chi2test (x, 'conditional', 1), 0.2303114201312508, 1e-14);
%!assert_equal (chi2test (x, 'conditional', 2), 0.0958810684407079, 1e-14);
%!assert_equal (chi2test (x, 'conditional', 3), 2.648037344954446e-11, 1e-14);
%!assert_equal (chi2test (x, 'homogeneous', []), 0.57357539370887889, -1e-10);
%!assert_equal (chi2test (x, 'homogeneous'), chi2test (x, 'homogeneous', []));
%!assert_equal (chi2test (x, 'mutual'), chi2test (x));
%!test
%! [pval, chisq, df, E] = chi2test (x);
%! assert_equal (chisq, 64.0982, 1e-4);
%! assert_equal (df, 7);
%! assert_equal (E(:,:,1), [42.903, 39.921; 17.185, 15.991], ones (2, 2) * 1e-3);
%!test
%! [pval, chisq, df, E] = chi2test (x, 'joint', 2);
%! assert_equal (chisq, 56.0943, 1e-4);
%! assert_equal (df, 5);
%! assert_equal (E(:,:,2), [40.922, 38.078; 23.310, 21.690], ones (2, 2) * 1e-3);
%!test
%! [pval, chisq, df, E] = chi2test (x, 'marginal', 3);
%! assert_equal (chisq, 51.0479, 1e-4);
%! assert_equal (df, 1);
%! assert_equal (E, [184.926, 172.074; 74.074, 68.926], ones (2, 2) * 1e-3);
%!test
%! [pval, chisq, df, E] = chi2test (x, 'conditional', 3);
%! assert_equal (chisq, 52.2509, 1e-4);
%! assert_equal (df, 3);
%! assert_equal (E(:,:,1), [53.345, 37.655; 14.655, 10.345], ones (2, 2) * 1e-3);
%!test
%! [pval, chisq, df, E] = chi2test (x, 'homogeneous', []);
%! assert_equal (chisq, 1.1117, 1e-4);
%! assert_equal (df, 2);
%! assert_equal (E(:,:,1), [60.469, 30.531; 7.531, 17.469], ones (2, 2) * 1e-3);
%!test
%! ## The homogeneous model reproduces every two-way margin
%! [~, ~, ~, E] = chi2test (x, 'homogeneous');
%! assert_equal ([sum(E, 1)(:); sum(E, 2)(:); sum(E, 3)(:)], ...
%!               [sum(x, 1)(:); sum(x, 2)(:); sum(x, 3)(:)], -1e-10);
%!test
%! ## E keeps the layout of a table whose dimensions differ
%! y = reshape ([12 7 3 9 15 4 6 8 11 5 14 2 9 10 3 7 6 13], [3, 2, 3]);
%! [~, ~, ~, E] = chi2test (y, 'joint', 2);
%! assert_equal (sum (E, 2), sum (y, 2), -1e-12);
%!test
%! y = reshape ([12 7 3 9 15 4 6 8 11 5 14 2 9 10 3 7 6 13], [3, 2, 3]);
%! [~, ~, ~, E] = chi2test (y, 'conditional', 2);
%! assert_equal (size (E), [3, 2, 3]);
%!test
%! ## Chi-squares of R's loglin on a 3-by-2-by-3 table
%! y = reshape ([12 7 3 9 15 4 6 8 11 5 14 2 9 10 3 7 6 13], [3, 2, 3]);
%! [~, c] = chi2test (y, 'homogeneous');
%! assert_equal (c, 15.05772742, -1e-9);
%!test
%! ## Marginal independence is tested in the collapsed table, as R's
%! ## chisq.test does it
%! y = reshape ([12 7 3 9 15 4 6 8 11 5 14 2 9 10 3 7 6 13], [3, 2, 3]);
%! [~, c, df] = chi2test (y, 'marginal', 2);
%! assert_equal ([c, df], [7.584460, 4], -1e-6);
