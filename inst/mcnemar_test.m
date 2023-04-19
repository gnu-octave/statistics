## Copyright (C) 1996-2017 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{h}, @var{pval}, @var{chisq}] =} mcnemar_test (@var{x})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{chisq}] =} mcnemar_test (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{chisq}] =} mcnemar_test (@var{x}, @var{testtype})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{chisq}] =} mcnemar_test (@var{x}, @var{alpha}, @var{testtype})
##
## Perform a McNemar's test on paired nominal data.
##
## @nospell{McNemar's} test is applied to a @math{2x2} contingency table @var{x}
## with a dichotomous trait, with matched pairs of subjects, of data
## cross-classified on the row and column variables to testing the null
## hypothesis of symmetry of the classification probabilities.  More formally,
## the null hypothesis of marginal homogeneity states that the two marginal
## probabilities for each outcome are the same.
##
## Under the null, with a sufficiently large number of discordants
## (@qcode{@var{x}(1,2) + @var{x}(2,1) >= 25}), the test statistic, @var{chisq},
## follows a chi-squared distribution with 1 degree of freedom.  When the number
## of discordants is less than 25, then the mid-P exact McNemar test is used.
##
## @var{testtype} will force @code{mcnemar_test} to apply a particular method
## for testing the null hypothesis independently of the number of discordants.
## Valid options for @var{testtype}:
## @itemize
## @item @qcode{"asymptotic"} Original McNemar test statistic
## @item @qcode{"corrected"} Edwards' version with continuity correction
## @item @qcode{"exact"} An exact binomial test
## @item @qcode{"mid-p"} The mid-P McNemar test (mid-p binomial test)
## @end itemize
##
## The test decision is returned in @var{h}, which is 1 when the null hypothesis
## is rejected (@qcode{@var{pval} < @var{alpha}}) or 0 otherwise.  @var{alpha}
## defines the critical value of statistical significance for the test.
##
## Further information about the McNemar's test can be found at
## @url{https://en.wikipedia.org/wiki/McNemar%27s_test}
##
## @seealso{crosstab, chi2test, fishertest}
## @end deftypefn

function [h, pval, chisq] = mcnemar_test (x, varargin)

  ## Check for valid number of input arguments
  if (nargin > 3)
    error ("mcnemar_test: too many input arguments.");
  endif

  ## Check contigency table
  if (! isequal (size (x), [2, 2]))
    error ("mcnemar_test: X must be a 2x2 matrix.");
  elseif (! (all ((x(:) >= 0)) && all (x(:) == fix (x(:)))))
    error ("mcnemar_test: all entries of X must be non-negative integers.");
  endif

  ## Add defaults
  alpha = 0.05;
  b = x(1,2);
  c = x(2,1);
  if (b + c < 25)
    testtype = "mid-p";
  else
    testtype = "asymptotic";
  endif

  ## Parse optional arguments
  if (nargin == 2)
    if (isnumeric (varargin{1}))
      alpha = varargin{1};
    elseif (ischar (varargin{1}))
      testtype = varargin{1};
    else
      error ("mcnemar_test: invalid 2nd input argument.");
    endif
  elseif (nargin == 3)
    alpha = varargin{1};
    testtype = varargin{2};
  endif

  ## Check optional arguments
  if (! isscalar (alpha) || alpha <= 0 || alpha >= 1)
    error ("mcnemar_test: invalid value for ALPHA.");
  endif
  types = {"exact", "asymptotic", "mid-p", "corrected"};
  if (! any (strcmpi (testtype, types)))
    error ("mcnemar_test: invalid value for TESTTYPE.");
  endif

  ## Calculate test
  switch (lower (testtype))
    case "asymptotic"
      chisq = (b - c) .^2 / (b + c);
      pval = 1 - chi2cdf (chisq, 1);
    case "corrected"
      chisq = (abs (b - c) - 1) .^2 / (b + c);
      pval = 1 - chi2cdf (chisq, 1);
    case "exact"
      chisq = [];
      pval = 2 * (binocdf (b, b + c, 0.5));
    case "mid-p"
      chisq = [];
      pval = 2 * (binocdf (b, b + c, 0.5)) - binopdf (b, b + c, 0.5);
  endswitch

  ## Get null hypothesis test result
  if (pval < alpha)
    h = 1;
  else
    h = 0;
  endif

endfunction

%!test
%! [h, pval, chisq] = mcnemar_test ([101,121;59,33]);
%! assert (h, 1);
%! assert (pval, 3.8151e-06, 1e-10);
%! assert (chisq, 21.356, 1e-3);

%!test
%! [h, pval, chisq] = mcnemar_test ([59,6;16,80]);
%! assert (h, 1);
%! assert (pval, 0.034690, 1e-6);
%! assert (isempty (chisq), true);

%!test
%! [h, pval, chisq] = mcnemar_test ([59,6;16,80], 0.01);
%! assert (h, 0);
%! assert (pval, 0.034690, 1e-6);
%! assert (isempty (chisq), true);

%!test
%! [h, pval, chisq] = mcnemar_test ([59,6;16,80], "mid-p");
%! assert (h, 1);
%! assert (pval, 0.034690, 1e-6);
%! assert (isempty (chisq), true);

%!test
%! [h, pval, chisq] = mcnemar_test ([59,6;16,80], "asymptotic");
%! assert (h, 1);
%! assert (pval, 0.033006, 1e-6);
%! assert (chisq, 4.5455, 1e-4);

%!test
%! [h, pval, chisq] = mcnemar_test ([59,6;16,80], "exact");
%! assert (h, 0);
%! assert (pval, 0.052479, 1e-6);
%! assert (isempty (chisq), true);

%!test
%! [h, pval, chisq] = mcnemar_test ([59,6;16,80], "corrected");
%! assert (h, 0);
%! assert (pval, 0.055009, 1e-6);
%! assert (chisq, 3.6818, 1e-4);

%!test
%! [h, pval, chisq] = mcnemar_test ([59,6;16,80], 0.1, "corrected");
%! assert (h, 1);
%! assert (pval, 0.055009, 1e-6);
%! assert (chisq, 3.6818, 1e-4);

%!error<mcnemar_test: too many input arguments.> mcnemar_test (59, 6, 16, 80)
%!error<mcnemar_test: X must be a 2x2 matrix.> mcnemar_test (ones (3, 3))
%!error<mcnemar_test: all entries of X must be non-negative integers.> ...
%! mcnemar_test ([59,6;16,-80])
%!error<mcnemar_test: all entries of X must be non-negative integers.> ...
%! mcnemar_test ([59,6;16,4.5])
%!error<mcnemar_test: invalid 2nd input argument.> ...
%! mcnemar_test ([59,6;16,80], {""})
%!error<mcnemar_test: invalid value for ALPHA.> ...
%! mcnemar_test ([59,6;16,80], -0.2)
%!error<mcnemar_test: invalid value for ALPHA.> ...
%! mcnemar_test ([59,6;16,80], [0.05, 0.1])
%!error<mcnemar_test: invalid value for ALPHA.> ...
%! mcnemar_test ([59,6;16,80], 1)
%!error<mcnemar_test: invalid value for TESTTYPE.> ...
%! mcnemar_test ([59,6;16,80], "")
