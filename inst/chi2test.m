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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} @var{pval} = chi2test (@var{x})
## @deftypefnx {Function File} [@var{pval}, @var{chisq}] = chi2test (@var{x})
## @deftypefnx {Function File} [@var{pval}, @var{chisq}, @var{dF}] = chi2test (@var{x})
## @deftypefnx {Function File} [@var{pval}, @var{chisq}, @var{dF}, @var{E}] = chi2test (@var{x})
## @deftypefnx {Function File} [@dots{}] = chi2test (@var{x}, "name", value)
##
## Perform a chi-squared test (for independence or homogeneity).
##
## For 2-way contingency tables, @code{chi2test} performs and a chi-squared test
## for independence or homogeneity, according to the sampling scheme and related
## question.  Independence means that the the two variables forming the 2-way
## table are not associated, hence you cannot predict from one another.
## Homogeneity refers to the concept of similarity, hence they all come from the
## same distribution.
##
## Both tests are computationally identical and will produce the same result.
## Nevertheless, they anwser to different questions.  Consider two variables,
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
## @multitable @columnfractions 0.05 0.1 0.85
## @item @tab @var{pval} @tab the p-value of the relevant test.
## @item @tab @var{chisq} @tab the chi^2 statistic of the relevant test.
## @item @tab @var{dF} @tab the degrees of freedom of the relevant test.
## @item @tab @var{E} @tab the EXPECTED values of the original contigency table.
## @end multitable
##
## Unlike MATLAB, in GNU Octave @code{chi2test} also supports 3-way tables,
## which involve three categorical variables (each in a different dimension of
## @var{x}.  In its simplest form, @code{[@dots{}] = chi2test (@var{x})} will
## will test for mutual independence among the three variables.  Alternatively,
## when called in the form @code{[@dots{}] = chi2test (@var{x}, "name", value)},
## it can perform the following tests:
##
## @multitable @columnfractions 0.2 0.1 0.7
## @headitem Name @tab Value @tab Description
## @item "mutual" @tab [] @tab Mutual independence.  All variables are
## independent from each other, (A, B, C).  Value must be an empty matrix.
## @item "joint" @tab scalar @tab Joint independence.  Two variables are jointly
## independent of the third, (AB, C). The scalar value corresponds to the
## dimension of the independent variable (i.e. 3 for C).
## @item "marginal" @tab scalar @tab Marginal independence.  Two variables are
## independent if you ignore the third, (A, C).  The scalar value corresponds
## to the dimension of the variable to be ignored (i.e. 2 for B).
## @item "conditional" @tab scalar @tab Conditional independence.  Two variables
## are independent given the third, (AC, BC).  The scalar value corresponds to
## the dimension of the variable that forms the conditional dependence
## (i.e. 3 for C).
## @item "homogeneous" @tab [] @tab Homogeneous associations.  Conditional
## (partial) odds-ratios are not related on the value of the third,
## (AB, AC, BC).  Value must be an empty matrix.
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
## @seealso{fishertest}
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
  ## Get size and dimensions of contigency table
  sz = size (x);
  dim = length (sz);
  ## Check optional arguments
  if (dim == 2 && nargin > 1)
    error ("chi2test: optional arguments are not supported for 2-way tables.");
  endif
  if (dim == 3 && mod (numel (varargin(:)), 2) != 0)
    error ("chi2test: optional arguments must be in pairs.");
  endif
  if (dim == 3 && nargin > 1 && ! isnumeric (varargin{2}))
    error (strcat (["chi2test: value must be numeric in optional argument"], ...
                   [" name/value pair, for 3-way tables."]));
  endif
  if (dim == 3 && nargin > 1 && numel (varargin{2}) > 1)
    error (strcat (["chi2test: value must be empty or scalar in optional"], ...
                   [" argument name/value pair, for 3-way tables."]));
  endif
  if (dim >= 4 && nargin > 1)
    error ("chi2test: optional arguments are not supported for k>3.");
  endif
  ## Calculate total sample size
  n = sum (x(:));
  ## For 2-way contigency table
  if (length (sz) == 2)
    ## Calculate degrees of freedom
    df = prod (sz - 1);
    ## Calculate expected values
    E = sum (x')' * sum (x) / n;
  ## For 3-way contigency table
  elseif (length (sz) == 3)
    ## Check optional arguments
    if (nargin == 1 || strcmpi (varargin{1}, "mutual"))
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
    elseif (strcmpi (varargin{1}, "joint"))
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
      x = permute (x, [c_dim, dm]);
      x = permute (x, [c_dim, dm]);
      E = permute (E, [c_dim, dm]);
      E = permute (E, [c_dim, dm]);
    elseif (strcmpi (varargin{1}, "marginal"))
      ## Get dimension of marginal variable (dim)
      c_dim = varargin{2};
      ## Calculate degrees of freedom
      c_sz = sz;
      c_sz(c_dim) = [];
      df = prod (sz) - sum (c_sz) + 1;
      ## Rearrange dimensions so that marginal variable goes in dim 1
      dm = [1, 2, 3];
      dm(c_dim) = [];
      x = permute (x, [c_dim, dm]);
      ## Calculate partial table sums
      q1 = sum (sum (x, 1), 3);
      q2 = sum (sum (x, 1), 2);
      n2 = sz(c_dim) * sum (x(:));
      ## Calculate expected values
      for d1 = 1:size (x, 1)
        for d2 = 1:size (x, 2)
          for d3 = 1:size (x, 3)
            E(d1,d2,d3) = q1(:,d2) * q2(:,:,d3) / n2;
          endfor
        endfor
      endfor
      ## Rearrange OBSERVED and EXPECTED matrices in original dimensions
      x = permute (x, [c_dim, dm]);
      E = permute (E, [c_dim, dm]);
    elseif (strcmpi (varargin{1}, "conditional"))
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
      x = permute (x, [c_dim, dm]);
      x = permute (x, [c_dim, dm]);
      E = permute (E, [c_dim, dm]);
      E = permute (E, [c_dim, dm]);
    elseif (strcmpi (varargin{1}, "homogeneous"))
      ## Calculate degrees of freedom
      df = prod(sz - 1);
      ## Compute observed marginal totals for any two dimensions
      omt12 = sum (sum (x, 3), 3);
      omt13 = sum (sum (x, 2), 2);
      omt23 = sum (sum (x, 1), 1);
      ## Produce initial seed 3-way table
      S = ones (sz);
      ## Calculate initial expected marginal totals
      emt12 = sum (sum (S, 3), 3);
      emt13 = sum (sum (S, 2), 2);
      emt23 = sum (sum (S, 1), 1);
      ## Compute difference to converge within certain tolerance or iterations
      OEdiff = sum (omt12(:) - emt12(:)) + sum (omt13(:) - emt13(:)) + ...
               sum (omt23(:) - emt23(:));
      iter = 1;
      tol = 1e-6;
      ## Start Iterative Proportional Fitting Procedure
      while (OEdiff > tol || iter > 50)
        ## Rows x Columns
        for d1 = 1:size (x, 1)
          for d2 = 1:size (x, 2)
            for d3 = 1:size (x, 3)
              E(d1,d2,d3) = S(d1,d2,d3) * omt12(d1,d2) / emt12(d1,d2);
            endfor
          endfor
        endfor
        ## Update seed and recalculate Rows x Layers expected marginal totals
        S = E;
        emt13 = sum (sum (S, 2), 2);
        ## Rows x Layers
        for d1 = 1:size (x, 1)
          for d2 = 1:size (x, 2)
            for d3 = 1:size (x, 3)
              E(d1,d2,d3) = S(d1,d2,d3) * omt13(d1,:,d3) / emt13(d1,:,d3);
            endfor
          endfor
        endfor
        ## Update seed and recalculate Columns x Layers expected marginal totals
        S = E;
        emt23 = sum (sum (S, 1), 1);
        ## Columns x Layers
        for d1 = 1:size (x, 1)
          for d2 = 1:size (x, 2)
            for d3 = 1:size (x, 3)
              E(d1,d2,d3) = S(d1,d2,d3) * omt23(:,d2,d3) / emt23(:,d2,d3);
            endfor
          endfor
        endfor
        ## Update seed and recalculate Rows x Layers expected marginal totals
        S = E;
        emt12 = sum (sum (S, 3), 3);
        ## Update difference between OBSERVED and EXPECTED tables
        OEdiff = sum (omt12(:) - emt12(:)) + sum (omt13(:) - emt13(:)) + ...
                 sum (omt23(:) - emt23(:));
        iter += 1;
      endwhile
    else
      error ("chi2test: invalid model name for testing a 3-way table.");
    endif
  ## For k-way contigency table, where k > 3
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
  pval = 1 - chi2cdf (chisq, df);
  ## Print results if no output requested
  if (nargout == 0)
    printf ("p-val = %f with chi^2 statistic = %f and d.f. = %d.\n", ...
            pval, chisq, df);
  endif
endfunction

## Input validation tests
%!error chi2test ();
%!error chi2test ([1, 2, 3, 4, 5]);
%!error chi2test ([1, 2; 2, 1+3i]);
%!error chi2test ([NaN, 6; 34, 12]);
%!error<chi2test: optional arguments are not supported for 2-way> ...
%! p = chi2test (ones (3, 3), "mutual", []);
%!error<chi2test: invalid model name for testing a 3-way table.> ...
%! p = chi2test (ones (3, 3, 3), "testtype", 2);
%!error<chi2test: optional arguments must be in pairs.> ...
%! p = chi2test (ones (3, 3, 3), "mutual");
%!error<chi2test: value must be numeric in optional argument> ...
%! p = chi2test (ones (3, 3, 3), "joint", ["a"]);
%!error<chi2test: value must be empty or scalar in optional argument> ...
%! p = chi2test (ones (3, 3, 3), "joint", [2, 3]);
%!error<chi2test: optional arguments are not supported for k> ...
%! p = chi2test (ones (3, 3, 3, 4), "mutual", [])

## Check warning
%!warning<chi2test: Expected values less than 5.> p = chi2test (ones (2));
%!warning<chi2test: Expected values less than 5.> p = chi2test (ones (3, 2));
%!warning<chi2test: Expected values less than 1.> p = chi2test (0.4 * ones (3));
## Output validation tests
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! p = chi2test (x);
%! assert (p, 0.017787, 1e-6);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [p, chisq] = chi2test (x);
%! assert (chisq, 11.9421, 1e-4);
%!test
%! x = [11, 3, 8; 2, 9, 14; 12, 13, 28];
%! [p, chisq, df] = chi2test (x);
%! assert (df, 4);
%!test
%!shared x
%! x(:,:,1) = [59, 32; 9,16];
%! x(:,:,2) = [55, 24;12,33];
%! x(:,:,3) = [107,80;17,56];%!
%!assert (chi2test (x), 2.282063427117009e-11, 1e-14);
%!assert (chi2test (x, "mutual", []), 2.282063427117009e-11, 1e-14);
%!assert (chi2test (x, "joint", 1), 1.164834895206468e-11, 1e-14);
%!assert (chi2test (x, "joint", 2), 7.771350230001417e-11, 1e-14);
%!assert (chi2test (x, "joint", 3), 0.07151361728026107, 1e-14);
%!assert (chi2test (x, "marginal", 1), 0, 1e-14);
%!assert (chi2test (x, "marginal", 2), 6.347555814301131e-11, 1e-14);
%!assert (chi2test (x, "marginal", 3), 0, 1e-14);
%!assert (chi2test (x, "conditional", 1), 0.2303114201312508, 1e-14);
%!assert (chi2test (x, "conditional", 2), 0.0958810684407079, 1e-14);
%!assert (chi2test (x, "conditional", 3), 2.648037344954446e-11, 1e-14);
%!assert (chi2test (x, "homogeneous", []), 0.4485579470993741, 1e-14);
%!test
%! [pval, chisq, df, E] = chi2test (x);
%! assert (chisq, 64.0982, 1e-4);
%! assert (df, 7);
%! assert (E(:,:,1), [42.903, 39.921; 17.185, 15.991], ones (2, 2) * 1e-3);
%!test
%! [pval, chisq, df, E] = chi2test (x, "joint", 2);
%! assert (chisq, 56.0943, 1e-4);
%! assert (df, 5);
%! assert (E(:,:,2), [40.922, 23.310; 38.078, 21.690], ones (2, 2) * 1e-3);
%!test
%! [pval, chisq, df, E] = chi2test (x, "marginal", 3);
%! assert (chisq, 146.6058, 1e-4);
%! assert (df, 9);
%! assert (E(:,1,1), [61.642; 57.358], ones (2, 1) * 1e-3);
%!test
%! [pval, chisq, df, E] = chi2test (x, "conditional", 3);
%! assert (chisq, 52.2509, 1e-4);
%! assert (df, 3);
%! assert (E(:,:,1), [53.345, 37.655; 14.655, 10.345], ones (2, 2) * 1e-3);
%!test
%! [pval, chisq, df, E] = chi2test (x, "homogeneous", []);
%! assert (chisq, 1.6034, 1e-4);
%! assert (df, 2);
%! assert (E(:,:,1), [60.827, 31.382; 7.173, 16.618], ones (2, 2) * 1e-3);
