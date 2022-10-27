## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
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
## @deftypefn  {Function File} @var{n} = sampsizepwr (@var{testtype}, @var{effsz})
## @deftypefnx {Function File} @var{n} = sampsizepwr (@var{testtype}, @var{effsz}, @var{pow})
## @deftypefnx {Function File} @var{n} = sampsizepwr (@var{testtype}, @var{effsz}, @var{pow}, @var{alpha})
## @deftypefnx {Function File} @var{n} = sampsizepwr (@var{testtype}, @var{effsz}, @var{pow}, @var{alpha}, @var{tails})
##
## Perform sample size calculations by power analysis. 
##
## @code{@var{n} = sampsizepwr (@var{testtype}, @var{effsz})} returns the
## required sample size to reach the significance level (alpha) of 0.05 in a
## two-tailed version of the test specified in @var{testtype} for the specified
## effect size, @var{effsz}, with a power of 0.8 (i.e. a type II error rate of
## 1 - 0.8 = 0.2)
##
## @var{testtype} can be:
##
## @itemize
## @item
## "t2" (default) : two-sample unpaired t-test
##
## @item
## "t" : paired t-test or one-sample t-test
##
## @item
## "r" : significance test for no correlation
##
## @item
## "z2" and "z" are also accepted values and represent z-test versions of "t2" 
## and "t" (using a Normal approximation instead of Student's t-distribution)
##
## @end itemize
##
## @var{effsz} can be numeric value corresponding to the standardized effect
## size, Cohen's d (when @var{testtype} is "t2" or "t") or Pearson's correlation
## coefficient (when @var{testtype} is "r"). For convenience, @var{effsz} can
## also be one of the following strings:
##
## @itemize
## @item
## "small" : which is 0.2 for Cohen's d or 0.1 for Pearson's r.
##
## @item
## "medium" : which is 0.5 for Cohen's d or 0.3 for Pearson's r.
##
## @item
## "large" : which is 0.8 for Cohen's d or 0.5 for Pearson's r.
##
## @end itemize
##
## @code{@var{n} = sampsizepwr (@var{testtype}, @var{effsz}, @var{pow})} also
## sets the desired power of the test. The power corresponds to 1 - beta, where
## beta is the type II error rate (i.e. the probability of not rejecting the
## null hypothesis when it is actually false). (Default is 0.8)
##
##
## @code{@var{n} = sampsizepwr (@var{testtype}, @var{effsz}, @var{pow}, @var{alpha})}
## also sets the desired significance level, @var{alpha}, of the test. @var{alpha}
## corresponds to the type I error ratepower corresponds to 1 - beta, where
## beta is the type II error rate (i.e. the probability of rejecting the
## null hypothesis when it is actually true). (Default is 0.05)
##
## HINT: If the test is expected to be among a family of tests, divide @var{alpha}
## by the number of tests so that the sample size calculations will maintain the
## desired power after correction for multiple comparisons.
##
## @code{@var{n} = sampsizepwr (@var{testtype}, @var{effsz}, @var{pow}, @var{alpha}, @var{tails})}
## also sets whether the test is one-sided or two-sided (Default is 2)
##
## @seealso{ztest, ttest, ttest2, corr}
## @end deftypefn

function n = sampsizepwr (testtype, effsz, power, alpha, tails, ncomp)

  ## Set default values
  if (nargin < 2)
    error ("sampsizepwr: at least two input arguments required")
  endif
  if (nargin < 3)
    power = 0.8; 
  end
  if (nargin < 4)
    alpha = 0.05;
  endif
  if (nargin < 5)
    tails = 2;
  endif

  ## Perform sample size calculation
  switch (lower (testtype))
    case {"z", "z2", "t", "t2"}
      ## Sample size calculations for the difference between means
      ## Assume effect size is Cohen's d
      if (isnumeric (effsz))
        d = abs (effsz);
      elseif (ischar (effsz))
        switch (lower (effsz))
          case "small"
            d = 0.2;
          case "medium"
            d = 0.5;
          case "large"
            d = 0.8;
          otherwise
            error ("sampsizepwr: string description for EFFSIZE not recognised")
        endswitch
      endif
      k = numel (testtype);
      ## Calculate group sample size based on Normal distribution
      n0 = (k * ((norminv (power) + norminv (1 - alpha / tails)) / d)^2);
      switch ( lower (testtype) )
        case {"z","z2"}
          n = ceil (n0);
        case {"t","t2"}
          ## Create function to optimize sample size based on Student-t distribution
          ## and n * type - type degrees of freedom
          func = @(n) n - k * ...
                  ((tinv (power, n * k - k) + ...
                    tinv (1 - alpha / tails, n * k - k)) / d)^2;
          n = ceil (fzero (func, n0)); # Find the root using fzero
      endswitch
    case "r"
      ## Sample size calculation for testing a correlation between 2 variables
      ## Assume effect size is Pearson's correlation coefficient (r)
      if (isnumeric (effsz))
        r = abs (effsz);
      elseif (ischar (effsz))
        switch (lower (effsz))
          case "small"
            r = 0.1;
          case "medium"
            r = 0.3;
          case "large"
            r = 0.5;
          otherwise
            error ("sampsizepwr: string description for EFFSIZE not recognised")
        endswitch
      endif
      z = 0.5 * log ((1 + r) / (1 - r));
      n = ceil (((norminv (power) + norminv (1 - alpha / tails)) / z)^2 + 3);
    otherwise
      error ("sampsizepwr: TESTTYPE not supported")
  endswitch

end

%!demo
%!
%! # The difference between a sample mean from a zero constant (one sample test)
%! # or the difference between two dependent means (matched pair)
%!
%! n = sampsizepwr ('t', "large")

%!demo
%!
%! # The difference between two independent means (two groups)
%!
%! n = sampsizepwr ('t2', "large")

%!demo
%!
%! # The test for Pearson's correlation coefficient equal to 0
%!
%! n = sampsizepwr ('r', "large")

%!test
%! # The difference between a sample mean from a zero constant (one sample test)
%! # or the difference between two dependent means (matched pair)
%! # Required sample sizes for small, medium and large effects with power, alpha 
%! # and the number of tails at 0.8, 0.05 and 2 respectively (defaults)
%! # The standardized effect size corresponds to Cohen's d
%! # Results compared to G*Power 3.1
%! ns = sampsizepwr ('t', 0.20, 0.80, 0.05, 2);
%! assert (ns, 199, 1);
%! nm = sampsizepwr ('t', 0.50, 0.80, 0.05, 2);
%! assert (nm, 34, 1);
%! nl = sampsizepwr ('t', 0.80, 0.80, 0.05, 2);
%! assert (nl, 15, 1);
%! ns = sampsizepwr ('t', "small");
%! assert (ns, 199, 1);
%! nm = sampsizepwr ('t', "medium");
%! assert (nm, 34, 1);
%! nl = sampsizepwr ('t', "large");
%! assert (nl, 15, 1);

%!test
%! # The difference between two independent means (two groups)
%! # Required sample sizes for small, medium and large effects with power, alpha 
%! # and the number of tails at 0.8, 0.05 and 2 respectively (defaults)
%! # The standardized effect size corresponds to Cohen's d
%! # Results compared to G*Power 3.1
%! ns = sampsizepwr ('t2', 0.20, 0.80, 0.05, 2);
%! assert (ns, 394, 1);
%! nm = sampsizepwr ('t2', 0.50, 0.80, 0.05, 2);
%! assert (nm, 64, 1);
%! nl = sampsizepwr ('t2', 0.80, 0.80, 0.05, 2);
%! assert (nl, 26, 1);
%! ns = sampsizepwr ('t2', "small");
%! assert (ns, 394, 1);
%! nm = sampsizepwr ('t2', "medium");
%! assert (nm, 64, 1);
%! nl = sampsizepwr ('t2', "large");
%! assert (nl, 26, 1);

%!test
%! # The test for Pearson's correlation coefficient equal to 0
%! # Required sample sizes for small, medium and large effects with power, alpha 
%! # and the number of tails at 0.8, 0.05 and 2 respectively (defaults)
%! # Results compared to G*Power 3.1
%! ns = sampsizepwr ('r', 0.10, 0.80, 0.05, 2);
%! assert (ns, 783, 1);
%! nm = sampsizepwr ('r', 0.30, 0.80, 0.05, 2);
%! assert (nm, 85, 1);
%! nl = sampsizepwr ('r', 0.50, 0.80, 0.05, 2);
%! assert (nl, 30, 1);
%! ns = sampsizepwr ('r', "small");
%! assert (ns, 783, 1);
%! nm = sampsizepwr ('r', "medium");
%! assert (nm, 85, 1);
%! nl = sampsizepwr ('r', "large");
%! assert (nl, 30, 1);