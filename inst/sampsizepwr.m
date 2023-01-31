## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
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
## @deftypefn  {statistics} @var{n} = sampsizepwr (@var{testtype}, @var{params}, @var{p1})
## @deftypefnx {statistics} @var{n} = sampsizepwr (@var{testtype}, @var{params}, @var{p1}, @var{power})
## @deftypefnx {statistics} @var{power} = sampsizepwr (@var{testtype}, @var{params}, @var{p1}, [], @var{n})
## @deftypefnx {statistics} @var{p1} = sampsizepwr (@var{testtype}, @var{params}, [], @var{power}, @var{n})
## @deftypefnx {statistics} [@var{n1}, @var{n2}] = sampsizepwr ("t2", @var{params}, @var{p1}, @var{power})
## @deftypefnx {statistics} [@dots{}] = sampsizepwr (@var{testtype}, @var{params}, @var{p1}, @var{power}, @var{n}, @var{name}, @var{value})
##
## Sample size and power calculation for hypothesis test.
##
## @code{sampsizepwr} computes the sample size, power, or alternative parameter
## value for a hypothesis test, given the other two values. For example, you can
## compute the sample size required to obtain a particular power for a
## hypothesis test, given the parameter value of the alternative hypothesis.
##
## @code{@var{n} = sampsizepwr (@var{testtype}, @var{params}, @var{p1})} returns
## the sample size N required for a two-sided test of the specified type to have
## a power (probability of rejecting the null hypothesis when the alternative is
## true) of 0.90 when the significance level (probability of rejecting the null
## hypothesis when the null hypothesis is true) is 0.05.  @var{params} specifies
## the parameter values under the null hypothesis.   P1 specifies the value of
## the single parameter being tested under the alternative hypothesis.  For the
## two-sample t-test, N is the value of the equal sample size for both samples,
## @var{params} specifies the parameter values of the first sample under the
## null and alternative hypotheses, and P1 specifies the value of the single
## parameter from the other sample under the alternative hypothesis.
##
## The following TESTTYPE values are available:
##
## @multitable @columnfractions 0.05 0.1 0.85
## @item @tab "z" @tab one-sample z-test for normally distributed data with
## known standard deviation.  @var{params} is a two-element vector [MU0 SIGMA0]
## of the mean and standard deviation, respectively, under the null hypothesis.
## P1 is the value of the mean under the alternative hypothesis.
## @item @tab "t" @tab one-sample t-test or paired t-test for normally
## distributed data with unknown standard deviation.  @var{params} is a
## two-element vector [MU0 SIGMA0] of the mean and standard deviation,
## respectively, under the null hypothesis.  P1 is the value of the mean under
## the alternative hypothesis.
## @item @tab "t2" @tab two-sample pooled t-test (test for equal means) for
## normally distributed data with equal unknown standard deviations.
## @var{params} is a two-element vector [MU0 SIGMA0] of the mean and standard
## deviation of the first sample under the null and alternative hypotheses.  P1
## is the the mean of the second sample under the alternative hypothesis.
## @item @tab "var" @tab chi-square test of variance for normally distributed
## data.  @var{params} is the variance under the null hypothesis.  P1 is the
## variance under the alternative hypothesis.
## @item @tab "p" @tab test of the P parameter (success probability) for a
## binomial distribution.  @var{params} is the value of P under the null
## hypothesis.  P1 is the value of P under the alternative hypothesis.
## @item @tab "r" @tab test of the correlation coefficient parameter for
## significance.  @var{params} is the value of r under the null hypothesis.
## P1 is the value of r under the alternative hypothesis.
## @end multitable
##
## The "p" test for the binomial distribution is a discrete test for which
## increasing the sample size does not always increase the power.  For N values
## larger than 200, there may be values smaller than the returned N value that
## also produce the desired power.
##
## @code{@var{n} = sampsizepwr (@var{testtype}, @var{params}, @var{p1},
## @var{power})} returns the sample size N such that the power is @var{power}
## for the parameter value P1.  For the two-sample t-test, N is the equal sample
## size of both samples.
##
## @code{[@var{n1}, @var{n2}] = sampsizepwr ("t2", @var{params}, @var{p1},
## @var{power})} returns the sample sizes @var{n1} and @var{n2} for the two
## samples.  These values are the same unless the "ratio" parameter,
## @code{@var{ratio} = @var{n2} / @var{n2}}, is set to a value other than
## the default (See the name/value pair definition of ratio below).
##
## @code{@var{power} = sampsizepwr (@var{testtype}, @var{params}, @var{p1}, [],
## @var{n})} returns the power achieved for a sample size of @var{n} when the
## true parameter value is @var{p1}.  For the two-sample t-test, @var{n} is the
## smaller one of the two sample sizes.
##
## @code{@var{p1} = sampsizepwr (@var{testtype}, @var{params}, [], @var{power},
## @var{n})} returns the parameter value detectable with the specified sample
## size @var{n} and power @var{power}.  For the two-sample t-test, @var{n} is
## the smaller one of the two sample sizes.  When computing @var{p1} for the "p"
## test, if no alternative can be rejected for a given @var{params}, @var{n} and
## @var{power} value, the function displays a warning message and returns NaN.
##
## @code{[@dots{}] = sampsizepwr (@dots{}, @var{n}, @var{name}, @var{value})}
## specifies one or more of the following @var{name} / @var{value} pairs:
##
## @multitable @columnfractions 0.05 0.1 0.85
## @item @tab "alpha" @tab significance level of the test (default is 0.05)
## @item @tab "tail" @tab the type of test which can be:
## @end multitable
## @multitable @columnfractions 0.1 0.15 0.75
## @item @tab "both" @tab two-sided test for an alternative @var{p1} not equal
## to @var{params}
## @item @tab "right" @tab one-sided test for an alternative @var{p1} larger
## than @var{params}
## @item @tab "left" @tab one-sided test for an alternative @var{p1} smaller
## than @var{params}
## @end multitable
## @multitable @columnfractions 0.05 0.1 0.85
## @item @tab "ratio" @tab desired ratio @var{n2} / @var{n2} of the larger
## sample size @var{n2} to the smaller sample size @var{n1}.  Used only for the
## two-sample t-test.  The value of @code{@var{ratio}} is greater than or equal
## to 1 (default is 1).
## @end multitable
##
## @code{sampsizepwr} computes the sample size, power, or alternative hypothesis
## value given values for the other two.  Specify one of these as [] to compute
## it.  The remaining parameters (and ALPHA, RATIO) can be scalars or arrays of
## the same size.
##
## @seealso{vartest, ttest, ttest2, ztest, binocdf}
## @end deftypefn

function [out, N2] = sampsizepwr (TestType, params, p1, power, n, varargin)
  ## Check for valid number of input arguments
  narginchk (3, Inf);
  ## Add defaults for 3 or 4 input arguments
  if (nargin == 3)
    power = 0.90;
    n = [];
  elseif (nargin == 4)
    n = [];
  endif
  ## Add defaults for extra arguments
  alpha = 0.05;
  tail = "both";
  ratio = 1;
  ## Check for valid test type and corresponding size of parameters
  t_types = {"z", "t", "t2", "var", "p", "r"};
  nparams = [2, 2, 2, 1, 1, 1];
  if (isempty (TestType) || ! ischar (TestType) || size (TestType, 1) != 1)
    error ("sampsizepwr: test type must be a non-empty string.");
  endif
  if (sum (strcmpi (TestType, t_types)) != 1)
    error ("sampsizepwr: invalid test type.");
  endif
  if (! isnumeric (params))
    error ("sampsizepwr: parameters must be numeric.");
  endif
  if (length (params) != nparams(strcmpi (TestType, t_types)))
    error ("sampsizepwr: invalid size of parameters for this test type.");
  endif
  ## Check for correct number of output arguments
  if ((nargout > 1) && ! strcmpi (TestType, "t2"))
    error ("sampsizepwr: wrong number of output arguments for this test type.");
  endif
  Lbound = [-Inf, -Inf, -Inf,   0, 0, -1];
  Ubound = [ Inf,  Inf,  Inf, Inf, 1,  1];
  Lbound = Lbound(strcmpi (TestType, t_types));
  Ubound = Ubound(strcmpi (TestType, t_types));
  ## Check for invalid parameters specific to each test type
  switch (lower (TestType))
    case "z"
      if (params(2) <= 0)
        error ("sampsizepwr: negative or zero variance.");
      endif
      PowerFunction = @PowerFunction_N;
    case "t"
      if (params(2) <= 0)
        error ("sampsizepwr: negative or zero variance.");
      endif
      PowerFunction = @PowerFunction_T;
    case "t2"
      if (params(2) <= 0)
        error ("sampsizepwr: negative or zero variance.");
      endif
      PowerFunction = @PowerFunction_T2;
    case "var"
      if (params(1) <= 0)
        error ("sampsizepwr: negative or zero variance.");
      endif
      PowerFunction = @PowerFunction_V;
    case "p"
      if (params(1) <= 0 || params(1) >= 1)
        error ("sampsizepwr: out of range probability.");
      endif
      PowerFunction = @PowerFunction_P;
    case "r"
      if (params(1) <= -1 || params(1) >= 1)
        error ("sampsizepwr: out of range regression coefficient.");
      elseif (params(1) == 0)
        error ("sampsizepwr: regression coefficient must not be 0.");
      endif
      PowerFunction = @PowerFunction_R;
  endswitch
  ## Get and validate optional parameters
  numarg = numel (varargin);
  argpos = 1;
  while (numarg)
    argname = varargin{argpos};
    switch (lower (argname))
      case "alpha"
        alpha = varargin{argpos + 1};
        if (! isnumeric (alpha) || any (alpha(:) <= 0) || any( alpha(:) >= 1))
          error ("sampsizepwr: invalid value for 'alpha' parameter.");
        endif
      case "tail"
        tail = varargin{argpos + 1};
        if (! ischar (tail) || (size (tail, 1) != 1))
          error ("sampsizepwr: 'tail' parameter must be a non-empty string.");
        endif
        if (sum (strcmpi (tail, {"left","both","right"})) != 1)
          error ("sampsizepwr: invalid value for 'tail' parameter.");
        endif
      case "ratio"
        ratio = varargin{argpos + 1};
        if (! isnumeric (ratio) || any (ratio(:) < 1))
          error ("sampsizepwr: invalid value for 'ratio' parameter.");
        endif
    endswitch
    numarg -= 2;
    argpos += 2;
  endwhile
  ## Check that only one of either p1, power, or n are missing
  if (isempty(p1) + isempty(power) + isempty(n) != 1)
    error ("sampsizepwr: only one of either p1, power, or n must be missing.");
  endif
  ## Check for valid P1
  if (! isempty (p1))
    if (! isnumeric (p1))
      error ("sampsizepwr: alternative hypothesis parameter must be numeric.");
    elseif ((! strcmpi (tail, "right") && any (p1(:) <= Lbound)) || ...
           (! strcmpi (tail, "left") && any (p1(:) >= Ubound)))
        error ("sampsizepwr: alternative hypothesis parameter out of range.");
    endif
  endif
  ## Check for valid POWER
  if (! isempty (power) && ...
     (! isnumeric (power) || any (power(:) <= 0) || any (power(:) >= 1)))
    error ("sampsizepwr: invalid value for POWER.");
  endif
  if (! isempty (power) && any (power(:) <= alpha(:)))
    error ("sampsizepwr: Cannot compute N or P1 unless POWER > 'alpha'.");
  endif
  ## Expand non-empty P1/POWER/N so they are all the same size
  if (isempty (p1))
    [err, power, n, alpha, ratio] = common_size (power, n, alpha, ratio);
    outclass = getclass (power, n, alpha, ratio);
  elseif (isempty (power))
    [err, p1, n, alpha, ratio] = common_size (p1, n, alpha, ratio);
    outclass = getclass (p1, n, alpha, ratio);
  else # n is empty
    [err, p1, power, alpha, ratio] = common_size (p1, power, alpha, ratio);
    outclass = getclass (power, p1, alpha, ratio);
  endif
  if (err > 0)
    error ("sampsizepwr: input arguments size mismatch.");
  end
  ## Check for valid options when computing N
  if (isempty (n))
    if (any (p1(:) == params(1)))
      error ("sampsizepwr: Same value for null and alternative hypothesis.");
    elseif (strcmpi (tail, "left") && any (p1(:) >= params(1)))
      error ("sampsizepwr: Invalid P1 for testing left tail.");
    elseif (strcmpi (tail, "right") && any(p1(:) <= params(1)))
      error ("sampsizepwr: Invalid P1 for testing right tail.");
    endif
  endif
  ## Allocate output of proper size and class
  out = zeros (size (alpha), outclass);
  ## Compute whichever one of P1/POWER/N that is now empty
  if (isempty (p1))
    ## Compute effect size given power and sample size
    switch (lower (TestType))
      case "z"
        ## z (normal) test
        out(:) = findP1z (params(1), params(2), power, n, alpha, tail);
      case "t"
        ## t-test
        out(:) = findP1t (params(1), params(2), power, n, alpha, tail);
      case "t2"
        ## two-sample t-test
        out(:) = findP1t2 (params(1), params(2), power, n, alpha, tail, ratio);
      case "var"
        ## chi-square (variance) test
        out(:) = findP1v (params(1), power, n, alpha, tail);
      case "p"
        ## binomial (p) test
        out(:) = findP1p (params(1), power, n, alpha, tail);
      case "r"
        ## regression coefficient (r) test
        out(:) = findP1r (params(1), power, n, alpha, tail);
    endswitch
  elseif (isempty (power))
    ## Compute power given effect size and sample size
    switch (lower (TestType))
      case {"z", "t"}
        out(:) = PowerFunction (params(1), p1, params(2), alpha, tail, n);
      case "t2"
        out(:) = PowerFunction (params(1), p1, params(2), alpha, tail, n, ratio);
      case {"var", "p", "r"}
        out(:) = PowerFunction (params(1), p1, alpha, tail, n);
    endswitch
  else
    ## Compute sample size given power and effect size
    switch (lower (TestType))
      case {"z", "t"}
        ## Calculate one-sided Z value directly
        out(:) = z1testN (params(1), p1, params(2), power, alpha, tail);
        ## Iterate upward from there for the other cases
        if (strcmpi (TestType, "t") || strcmp (tail, "both"))
          if (strcmpi (TestType, "t"))
            out = max (out, 2);
          endif
          ## Count upward until we get the value we need
          elem = 1:numel (alpha);
          while (! isempty (elem))
            actualpower = PowerFunction_T (params(1), p1(elem), params(2), ...
                                           alpha(elem), tail, out(elem));
            elem = elem(actualpower < power(elem));
            out(elem) = out(elem) + 1;
          endwhile
        endif
      case "t2"
        ## Initialize second output argument
        N2 = zeros (size (alpha), outclass);
        ## Caculate one-sided two-sample t-test iteratively
        [out(:), N2(:)] = t1testN (params(1), p1, params(2), power, ...
                                   alpha, tail, ratio);
      case "var"
        ## Use a binary search method
        out(:) = searchbinaryN (PowerFunction, [1, 100], params(1), ...
                                p1, power, alpha, tail);
      case "p"
        ## Use a binary search method
        out(:) = searchbinaryN (PowerFunction, [0, 100], params(1), ...
                                p1, power, alpha, tail);
        ## Adjust for discrete distribution
        t = out <= 200;
        if (any (t(:)))
          ## Try values from 1 up to N (out) and pick the smallest value
          out(t) = adjdiscreteN (out(t), PowerFunction, params(1), ...
                                 p1(t), alpha(t), tail, power(t));
        endif
        if (any (! t(:)))
          warning ("sampsizepwr: approximate N.");
        endif
      case "r"
        ## Calculate sample size using Student's t distribution
        out(:) = r1testN (params(1), p1, power, alpha, tail);
    endswitch
  endif
endfunction

## Define class for output
function out = getclass (varargin)
  if (class (varargin{1}) == "single" || class (varargin{2}) == "single" || ...
      class (varargin{3}) == "single" || class (varargin{4}) == "single")
    out = "single";
  else
    out = "double";
  endif
endfunction

## Sample size calculation for the one-sided Z test
function N = z1testN (mu0, mu1, sig, desiredpower, alpha, tail)
  ## Compute the one-sided normal value directly
  if (strcmp (tail, "both"))
    alpha = alpha ./ 2;
  endif
  z1 = -norminv (alpha);
  z2 = norminv (1 - desiredpower);
  mudiff = abs (mu0 - mu1) / sig;
  N = ceil (((z1 - z2) ./ mudiff) .^ 2);
endfunction

## Sample size calculation for R test
function N = rtestN (r0, r1, desiredpower, alpha, tail)
  ## Compute only for 2-tailed test
  if (strcmp (tail, "both"))
    alpha = alpha ./ 2;
  else
    error ("sampsizepwr: only 2-tailed testing for regression coefficient.");
  endif
  ## Get quantiles of the standard normal deviates for alpha and power
  Za = norminv (alpha);
  Zb = norminv (1 - desiredpower);
  ## Compute difference in regression coefficients
  rdiff = abs (r0 - r1)
  C = 0.5 * log ((1 + rdiff) / (1 - rdiff));
  ## Compute sample size
  N = ((Za + Zb) / C) .^ 2 + 3;
endfunction

## Find alternative hypothesis parameter value P1 for Z test
function mu1 = findP1z (mu0, sig, desiredpower, N, alpha, tail)
  if (strcmp (tail, "both"))
    alpha = alpha ./ 2;
  end
  sig = sig ./ sqrt (N);
  ## Get quantiles of the normal or t distribution
  if (strcmp (tail, "left"))
    z1 = norminv (alpha);
    z2 = norminv (desiredpower);
  else               # upper or two-tailed test
    z1 = norminv (1 - alpha);
    z2 = norminv (1 - desiredpower);
  endif
  mu1 = mu0 + sig .* (z1 - z2);
  ## For 2-sided test, refine by taking the other tail into account
  if (strcmp (tail, "both"))
    elem = 1:numel (alpha);
    desiredbeta = 1 - desiredpower;
    betahi = desiredbeta;
    betalo = zeros (size (desiredbeta));
    while (true)
      ## Compute probability of being below the lower critical value under H1
      betalo(elem) = normcdf (-z1(elem) + (mu0 - mu1(elem)) ./ sig(elem));
      ## See if the upper and lower probabilities are close enough
      elem = elem(abs ((betahi(elem) - betalo(elem)) - desiredbeta(elem)) > ...
                  1e-6 * desiredbeta(elem));
      if (isempty (elem))
        break
      endif
      ## Find a new mu1 by adjusting beta to take lower tail into account
      betahi(elem) = desiredbeta(elem) + betalo(elem);
      mu1(elem) = mu0 + sig(elem) .* (z1(elem) - norminv (betahi(elem)));
    endwhile
  endif
endfunction

## Find alternative hypothesis parameter value P1 for t-test
function mu1 = findP1t (mu0, sig, desiredpower, N, alpha, tail)
  if (strcmp (tail, "both"))
    a2 = alpha ./ 2;
  else
    a2 = alpha;
  endif
  ## Get quantiles of the normal or t distribution
  if (strcmp (tail, "left"))
    z1 = norminv(alpha);
    z2 = norminv(desiredpower);
  else               # upper or two-tailed test
    z1 = norminv(1-a2);
    z2 = norminv(1-desiredpower);
  endif
  mu1 = mu0 + sig .* (z1-z2) ./ sqrt (N);
  ## Refine using fzero
  for j=1:numel (mu1)
    if (mu1(j) > mu0)
      F0 = @(mu1arg) PowerFunction_T (mu0, max (mu0, mu1arg), sig, alpha(j), ...
                                      tail, N(j)) - desiredpower(j);
    else
      F0 = @(mu1arg) desiredpower(j) - PowerFunction_T (mu0, min (mu0, ...
                                       mu1arg), sig, alpha(j), tail, N(j));
    endif
    mu1(j) = fzero (F0, mu1(j));
  endfor
endfunction

## Sample size calculation for the one-sided two-sample t-test
function [N1, N2] = t1testN (mu0, mu1, sig, desiredpower, alpha, tail, ratio)
  if (strcmp (tail, "both"))
    alpha = alpha ./ 2;
  endif
  ## Compute the initial value of N, approximated by normal distribution
  z1 = -norminv (alpha);
  z2 = norminv (1 - desiredpower);
  n_0 = ceil ((z1 - z2) .^2 .* (sig ./ abs ((mu0 - mu1))) .^ 2 * 2);
  ## n need to be > 1, otherwise the degree of freedom of t < 0
  n_0(n_0 <= 1) = 2;
  N = ones (size (n_0));
  ## iteratively update the sample size
  if (strcmp (tail, "both"))
    for j = 1:numel (n_0)
      F = @(n) nctcdf (tinv (alpha(j), n + ratio(j) .* n - 2), ...
                       n + ratio(j) .* n - 2, abs (mu1(j) - mu0) ./ ...
                       (sig .* sqrt (1 ./ n + 1 ./ (ratio(j) .* n)))) + ...
                       (1 - nctcdf (- tinv (alpha(j), n + ratio(j) .* n - 2), ...
                       n + ratio(j) .* n - 2, abs (mu1(j) - mu0) ./ ...
                       (sig .* sqrt (1 ./ n + 1 ./ (ratio(j) .* n)))))- ...
                       desiredpower(j);
      N(j) = localfzero (F, n_0(j), ratio);
    endfor
  else
    for j = 1:numel (n_0)
      F = @(n) (1 - nctcdf (- tinv (alpha(j), n + ratio(j) .* n - 2), ...
                n + ratio(j) .* n - 2, abs (mu1(j) - mu0) ./ (sig .* ...
                sqrt (1 ./ n + 1 ./ (ratio(j) .* n))))) - desiredpower(j);
      N(j) = localfzero (F, n_0(j), ratio);
    endfor
  endif
  N1 = ceil (N);
  N2 = ceil (ratio .* N);
endfunction

## Find alternative hypothesis parameter value P1 for two-sample t-test
function mu1 = findP1t2 (mu0, sig, desiredpower, N,  alpha, tail, ratio)
  if (strcmp (tail, "both"))
    a2 = alpha ./ 2;
  else
    a2 = alpha;
  endif
  ## Get quantiles of the normal or t distribution
  if (strcmp (tail, "left"))
    t1 = tinv (alpha, N + ratio .* N - 2);
    t2 = tinv (desiredpower, N + ratio .* N - 2);
  else               # upper or two-tailed test
    t1 = tinv (1 - a2, N + ratio .* N - 2);            # upper tail under H0
    t2 = tinv (1 - desiredpower, N + ratio .* N - 2);  # lower tail under H1
  endif
  mu1 = mu0 + sig .* (t1 - t2) .* sqrt (1 ./ N + 1 ./ (ratio .* N));
  ## Refine using fzero
  for j = 1:numel (mu1)
    if (mu1(j) > mu0)
        F0 = @(mu1arg) PowerFunction_T2 (mu0, max (mu0, mu1arg), sig, ...
               alpha(j), tail, N(j), ratio(j)) - desiredpower(j);
    else
        F0 = @(mu1arg) desiredpower(j) - PowerFunction_T2 (mu0, min (mu0, ...
               mu1arg), sig, alpha(j), tail, N(j), ratio(j));
    endif
    mu1(j) = fzero (F0, mu1(j));
  endfor
endfunction

## Find alternative hypothesis parameter value P1 for variance test
function p1 = findP1v (p0, desiredpower, N, alpha, tail)
  ## F and Finv are the cdf and inverse cdf
  F = @(x,n,p1) chi2cdf (x .* (n - 1) ./ p1, n - 1);     # cdf for s^2
  Finv = @(p,n,p1) p1 .* chi2inv (p, n - 1) ./ (n - 1);  # inverse

  if (strcmp (tail, "both"))
    alpha = alpha ./ 2;
  endif
  desiredbeta = 1 - desiredpower;

  ## Calculate critical values and p1 for one-sided test
  if (! strcmp (tail, "left"))
    critU = Finv (1 - alpha, N, p0);
    p1 = 1 ./ Finv (desiredbeta, N, 1 ./ critU);
  endif
  if (! strcmp (tail, "right"))
    critL = Finv (alpha, N, p0);
  endif
  if (strcmp (tail, "left"))
    p1 = 1 ./ Finv (desiredpower, N, 1 ./ critL);
  endif

  if (strcmp (tail, "both"))
    ## For 2-sided test, we have the upper tail probability under H1.
    ## Refine by taking the other tail into account.
    elem = 1:numel (alpha);
    betahi = desiredbeta;
    betalo = zeros (size (desiredbeta));
    while (true)
      ## Compute probability of being in the lower tail under H1
      betalo(elem) = F (critL(elem), N(elem), p1(elem));
      ## See if the upper and lower probabilities are close enough
      obsbeta = betahi(elem) - betalo(elem);
      elem = elem(abs (obsbeta - desiredbeta(elem)) > 1e-6 * desiredbeta(elem));
      if (isempty (elem))
        break
      endif
      ## Find a new mu1 by adjusting beta to take lower tail into account
      betahi(elem) = desiredbeta(elem) + betalo(elem);
      p1(elem) = 1 ./ Finv (betahi(elem), N(elem), 1 ./ critU(elem));
    endwhile
  endif
endfunction

## Find alternative hypothesis parameter value P1 for p test
function p1 = findP1p (p0, desiredpower, N, alpha, tail)
  ## Get critical values
  [critL, critU] = getcritP (p0, N, alpha, tail);
  ## Use a normal approximation to find P1 values
  sigma = sqrt (p0 .* (1 - p0) ./ N);
  p1 = findP1z (p0, sigma, desiredpower, N, alpha, tail);
  ## Problem if we have no critical region left
  if (strcmp (tail, "both"))
    t = (critL == 0 & critU == N);
  elseif (strcmp (tail, "right"))
    t = (critU == N);
  else
    t = (critL == 0);
  endif
  if (any (t))
    warning ("sampsizepwr: No Valid Parameter");
    p1(t) = NaN;
  endif
  ## Force in bounds
  t = p1 <= 0;
  if (any (t(:)))
    p1(t) = p0 / 2;
  end
  t = p1 >= 1;
  if (any (t(:)))
    p1(t) = 1 - p0 / 2;
  end
  ## Refine using fzero
  for j=1:numel(p1)
    if (! isnan (p1(j)));
      if (p1(j) > p0)
        F0 = @(p1arg) PowerFunction_P (p0, max (p0, min (1, p1arg)), ...
               alpha(j), tail, N(j), critL(j), critU(j)) - desiredpower(j);
      else
        F0 = @(p1arg) desiredpower(j) - PowerFunction_P (p0, max (0, ...
               min (p0, p1arg)), alpha(j), tail, N(j), critL(j), critU(j));
      endif
      p1(j) = fzero (F0, p1(j));
    endif
  endfor
endfunction

## Find alternative hypothesis parameter value P1 for r test
function p1 = findP1r (p0, desiredpower, N, alpha, tail)
  ## Compute only for 2-tailed test
  if (! strcmp (tail, "both"))
    error ("sampsizepwr: only 2-tailed testing for regression coefficient.");
  endif
  ## Set initial search boundaries for p1
  p1_lo = eps;
  p1_hi = 1 - eps;
  ## Compute initial sample size N0 according to P0, POWER and ALPHA
  N0 = rtestN (p0, 0, desiredpower, alpha, tail);
  ## Find P0 for N0 == N
  while (N != N0)
    if (N0 < N)
      p1_hi = p0;
      p1 = (p0 + p1_lo) / 2;
      p0 = p1;
      N0 = rtestN (p1, 0, desiredpower, alpha, tail);
    else
      p1_lo = p0;
      p1 = (p0 + p1_hi) / 2;
      p0 = p1;
      N0 = rtestN (p1, 0, desiredpower, alpha, tail);
    endif
  endwhile
endfunction

## Get upper and lower critical values for binomial (p) test.
function [critL, critU] = getcritP (p0, N, alpha, tail)
  ## For two-sided tests, this function tries to compute critical values
  ## favorable for p0<.5.  It does this by allocating alpha/2 to the lower
  ## tail where the probabilities come in larger chunks, then using any
  ## left-over alpha, probably more than alpha/2, for the upper tail.

  ## Get part of alpha available for lower tail
  if (strcmp (tail, "both"))
    Alo = alpha ./ 2;
  elseif (strcmp (tail, "left"))
    Alo = alpha;
  else
    Alo = 0;
  endif
  ## Calculate critical values
  critU = N;
  critL = zeros(size(N));
  if (! strcmp (tail, "right"))
    critL = binoinv (Alo, N, p0);
    Alo = binocdf (critL, N, p0);
    t = (critL < N) & (Alo <= alpha / 2);
    critL(t) = critL(t) + 1;
    Alo(! t) = Alo(! t) - binopdf (critL(! t), N(! t), p0);
  endif
  if (! strcmp (tail, "left"))
    Aup = max(0, alpha - Alo);
    critU = binoinv(1 - Aup, N, p0);
  endif
endfunction

## Sample size calculation via binary search
function N = searchbinaryN (F, lohi, p0, p1, desiredpower, alpha, tail)
  ## Find uper and lower bounds
  nlo = repmat(lohi(1),size(alpha));
  nhi = repmat(lohi(2),size(alpha));
  obspower = F(p0,p1,alpha,tail,nhi);
  ## Iterate on n until we achieve the desired power
  elem = 1:numel (alpha);
  while (! isempty (elem))
    elem = elem(obspower(elem) < desiredpower(elem));
    nhi(elem) = nhi(elem) * 2;
    obspower(elem) = F (p0, p1(elem), alpha(elem), tail, nhi(elem));
  endwhile
  ## Binary search between these bounds for required sample size
  elem = find(nhi > nlo+1);
  while (! isempty (elem))
    n = floor ((nhi(elem) + nlo(elem)) / 2);
    obspower = F (p0, p1(elem), alpha(elem), tail, n);
    toohigh = (obspower > desiredpower(elem));
    nhi(elem(toohigh)) = n(toohigh);
    nlo(elem(! toohigh)) = n(! toohigh);
    elem = elem(nhi(elem) > nlo(elem) + 1);
  endwhile
  N = nhi;
endfunction

## Adjust sample size to take discreteness into account
function N = adjdiscreteN (N, PowerFunction, p0, p1, alpha, tail, power)
  for j=1:numel(N)
    allN = 1:N(j);
    obspower = PowerFunction (p0, p1(j), alpha(j), tail, allN);
    N(j) = allN(find (obspower >= power(j), 1, "first"));
  endfor
endfunction

## Normal power calculation
function power = PowerFunction_N (mu0, mu1, sig, alpha, tail, n)
  S = sig ./ sqrt (n);
  if (strcmp (tail, "both"))
    critL = norminv (alpha / 2, mu0, S);
    critU = mu0 + (mu0 - critL);
    power = normcdf (critL, mu1, S) + normcdf (-critU, -mu1, S);
  elseif (strcmp (tail, "right"))
    crit = mu0 + (mu0 - norminv (alpha, mu0, S));
    power = normcdf (-crit, -mu1, S);
  else
    crit = norminv (alpha, mu0, S);
    power = normcdf (crit, mu1, S);
  endif
endfunction

## T power calculation
function power = PowerFunction_T (mu0, mu1, sig, alpha, tail, n)
  S = sig ./ sqrt (n);
  ncp = (mu1 - mu0) ./ S;
  if (strcmp (tail, "both"))
    critL = tinv (alpha / 2, n - 1);
    critU = -critL;
    power = nctcdf (critL, n - 1, ncp) + nctcdf (-critU, n - 1, -ncp);
  elseif (strcmp (tail, "right"))
    crit = tinv (1 - alpha, n - 1);
    power = nctcdf (-crit, n - 1, -ncp);
  else
    crit = tinv (alpha, n - 1);
    power = nctcdf (crit, n - 1, ncp);
  endif
endfunction

## Two-sample T power calculation
function power = PowerFunction_T2 (mu0, mu1, sig, alpha, tail, n, ratio)
  ncp = (mu1 - mu0) ./ (sig .* sqrt (1 ./  n + 1 ./ (ratio .* n)));
  if (strcmp (tail, "both"))
    critL = tinv (alpha / 2, n + ratio .* n - 2);
    critU = -critL;
    power = nctcdf (critL, n + ratio .* n - 2, ncp) + ...
            nctcdf (-critU, n + ratio .* n - 2, -ncp);
  elseif (strcmp (tail, "right"))
    crit = tinv (1 - alpha, n + ratio .* n - 2);
    power = nctcdf (-crit, n + ratio .* n - 2, -ncp);
  else
    crit = tinv (alpha, n  + ratio .* n - 2);
    power = nctcdf (crit, n + ratio .* n - 2, ncp);
  endif
endfunction

## Chi-square power calculation
function power = PowerFunction_V (v0, v1, alpha, tail, n)
  if (strcmp (tail, "both"))
   critU = v0 .* chi2inv (1 - alpha / 2, n - 1);
   critL = v0 .* chi2inv (alpha / 2, n - 1);
   power = chi2cdf (critL ./ v1, n - 1) + chi2cdf (critU ./ v1, n - 1);
  elseif (strcmp (tail, "right"))
    crit = v0 .* chi2inv (1 - alpha, n - 1);
    power = chi2cdf (crit ./ v1, n - 1);
  else
    crit = v0 .* chi2inv (alpha, n - 1);
    power = chi2cdf (crit ./ v1, n - 1);
  endif
endfunction

## Binomial power calculation
function [power, critL, critU] = PowerFunction_P (p0, p1, alpha, ...
                                                  tail, n, critL, critU)
  if (nargin < 6)
    [critL, critU] = getcritP (p0, n, alpha, tail);
  endif
  if (strcmp (tail, "both"))
    power = binocdf (critL - 1, n, p1) + 1 - binocdf (critU, n, p1);
  elseif (strcmp (tail, "right"))
    power = 1 - binocdf (critU , n, p1);
  else
    power = binocdf (critL - 1, n, p1);
  endif
endfunction

## Regression power calculation
function power = PowerFunction_R (r0, r1, alpha, tail, n)
  ## Compute only for 2-tailed test
  if (! strcmp (tail, "both"))
    error ("sampsizepwr: only 2-tailed testing for regression coefficient.");
  endif
  ## Set initial search boundaries for power
  dp_lo = eps;
  dp_hi = 1 - eps;
  power = 0.5;
  ## Compute initial sample size N0 according to P0, POWER and ALPHA
  N0 = rtestN (r0, r1, power, alpha, tail);
  ## Find POWER for N0 == N
  while (N != N0)
    if (N0 < N)
      dp_hi = power;
      power = (power + dp_lo) / 2;
      N0 = rtestN (r0, r1, power, alpha, tail);
    else
      dp_lo = power;
      power = (power + pd_hi) / 2;
      N0 = rtestN (r0, r1, power, alpha, tail);
    endif
  endwhile
endfunction

## Local zero function for "t2" test
function N = localfzero (F, N0, ratio)
  ## Set minN according to ratio
  if (ratio >= 2)
    minN = 1;
  else
    minN = 2;
  endif
  ## Return minN if function gives a value above zero
  if (F(minN) > 0)
    N = minN;
    return;
  endif
  ## Make sure that fzero does not try values below minN
  if (N0 == minN)
    N0 = N0 + 1;
  endif
  ## Find solution
  if (F(N0) > 0)
    N = fzero (F, [minN, N0], optimset ('TolX',1e-6)); # N0 is an upper bound
  else
    N = fzero (F, N0, optimset ('TolX',1e-6));         # N0 is a starting value
  endif
endfunction

## Demos
%!demo
%! ## Compute the mean closest to 100 that can be determined to be
%! ## significantly different from 100 using a t-test with a sample size
%! ## of 60 and a power of 0.8.
%! mu1 = sampsizepwr ("t", [100, 10], [], 0.8, 60);
%! disp (mu1);

%!demo
%! ## Compute the sample sizes required to distinguish mu0 = 100 from
%! ## mu1 = 110 by a two-sample t-test with a ratio of the larger and the
%! ## smaller sample sizes of 1.5 and a power of 0.6.
%! [N1,N2] = sampsizepwr ("t2", [100, 10], 110, 0.6, [], "ratio", 1.5)

%!demo
%! ## Compute the sample size N required to distinguish p=.26 from p=.2
%! ## with a binomial test.  The result is approximate, so make a plot to
%! ## see if any smaller N values also have the required power of 0.6.
%! Napprox = sampsizepwr ("p", 0.2, 0.26, 0.6);
%! nn = 1:250;
%! pwr = sampsizepwr ("p", 0.2, 0.26, [], nn);
%! Nexact = min (nn(pwr >= 0.6));
%! plot(nn,pwr,'b-', [Napprox Nexact],pwr([Napprox Nexact]),'ro');
%! grid on

%!demo
%! ## The company must test 52 bottles to detect the difference between a mean
%! ## volume of 100 mL and 102 mL with a power of 0.80.  Generate a power curve
%! ## to visualize how the sample size affects the power of the test.
%!
%! nout = sampsizepwr('t',[100 5],102,0.80);
%! nn = 1:100;
%! pwrout = sampsizepwr('t',[100 5],102,[],nn);
%!
%! figure;
%! plot (nn, pwrout, "b-", nout, 0.8, "ro")
%! title ("Power versus Sample Size")
%! xlabel ("Sample Size")
%! ylabel ("Power")

## Input validation
%!error<sampsizepwr: test type must be a non-empty string.> ...
%! out = sampsizepwr ([], [100, 10], [], 0.8, 60);
%!error<sampsizepwr: test type must be a non-empty string.> ...
%! out = sampsizepwr (3, [100, 10], [], 0.8, 60);
%!error<sampsizepwr: test type must be a non-empty string.> ...
%! out = sampsizepwr ({"t", "t2"}, [100, 10], [], 0.8, 60);
%!error<sampsizepwr: invalid test type.> ...
%! out = sampsizepwr ("reg", [100, 10], [], 0.8, 60);
%!error<sampsizepwr: parameters must be numeric.> ...
%! out = sampsizepwr ("t", ["a", "e"], [], 0.8, 60);
%!error<sampsizepwr: invalid size of parameters for this test type.> ...
%! out = sampsizepwr ("z", 100, [], 0.8, 60);
%!error<sampsizepwr: invalid size of parameters for this test type.> ...
%! out = sampsizepwr ("t", 100, [], 0.8, 60);
%!error<sampsizepwr: invalid size of parameters for this test type.> ...
%! out = sampsizepwr ("t2", 60, [], 0.8, 60);
%!error<sampsizepwr: invalid size of parameters for this test type.> ...
%! out = sampsizepwr ("var", [100, 10], [], 0.8, 60);
%!error<sampsizepwr: invalid size of parameters for this test type.> ...
%! out = sampsizepwr ("p", [100, 10], [], 0.8, 60);
%!error<sampsizepwr: invalid size of parameters for this test type.> ...
%! out = sampsizepwr ("r", [100, 10], [], 0.8, 60);
%!error<sampsizepwr: wrong number of output arguments for this test type.> ...
%! [out, N1] = sampsizepwr ("z", [100, 10], [], 0.8, 60);
%!error<sampsizepwr: wrong number of output arguments for this test type.> ...
%! [out, N1] = sampsizepwr ("t", [100, 10], [], 0.8, 60);
%!error<sampsizepwr: wrong number of output arguments for this test type.> ...
%! [out, N1] = sampsizepwr ("var", 2, [], 0.8, 60);
%!error<sampsizepwr: wrong number of output arguments for this test type.> ...
%! [out, N1] = sampsizepwr ("p", 0.1, [], 0.8, 60);
%!error<sampsizepwr: wrong number of output arguments for this test type.> ...
%! [out, N1] = sampsizepwr ("r", 0.5, [], 0.8, 60);
%!error<sampsizepwr: negative or zero variance.> ...
%! out = sampsizepwr ("z", [100, 0], [], 0.8, 60);
%!error<sampsizepwr: negative or zero variance.> ...
%! out = sampsizepwr ("z", [100, -5], [], 0.8, 60);
%!error<sampsizepwr: negative or zero variance.> ...
%! out = sampsizepwr ("t", [100, 0], [], 0.8, 60);
%!error<sampsizepwr: negative or zero variance.> ...
%! out = sampsizepwr ("t", [100, -5], [], 0.8, 60);
%!error<sampsizepwr: negative or zero variance.> ...
%! [out, N1] = sampsizepwr ("t2", [100, 0], [], 0.8, 60);
%!error<sampsizepwr: negative or zero variance.> ...
%! [out, N1] = sampsizepwr ("t2", [100, -5], [], 0.8, 60);
%!error<sampsizepwr: negative or zero variance.> ...
%! out = sampsizepwr ("var", 0, [], 0.8, 60);
%!error<sampsizepwr: negative or zero variance.> ...
%! out = sampsizepwr ("var", -5, [], 0.8, 60);
%!error<sampsizepwr: out of range probability.> ...
%! out = sampsizepwr ("p", 0, [], 0.8, 60);
%!error<sampsizepwr: out of range probability.> ...
%! out = sampsizepwr ("p", 1.2, [], 0.8, 60);
%!error<sampsizepwr: out of range regression coefficient.> ...
%! out = sampsizepwr ("r", -1.5, [], 0.8, 60);
%!error<sampsizepwr: out of range regression coefficient.> ...
%! out = sampsizepwr ("r", -1, [], 0.8, 60);
%!error<sampsizepwr: out of range regression coefficient.> ...
%! out = sampsizepwr ("r", 1.2, [], 0.8, 60);
%!error<sampsizepwr: regression coefficient must not be 0.> ...
%! out = sampsizepwr ("r", 0, [], 0.8, 60);
%!error<sampsizepwr: invalid value for 'alpha' parameter.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "alpha", -0.2);
%!error<sampsizepwr: invalid value for 'alpha' parameter.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "alpha", 0);
%!error<sampsizepwr: invalid value for 'alpha' parameter.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "alpha", 1.5);
%!error<sampsizepwr: invalid value for 'alpha' parameter.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "alpha", "zero");
%!error<sampsizepwr: 'tail' parameter must be a non-empty string.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "tail", 1.5);
%!error<sampsizepwr: 'tail' parameter must be a non-empty string.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "tail", {"both", "left"});
%!error<sampsizepwr: invalid value for 'tail' parameter.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "tail", "other");
%!error<sampsizepwr: invalid value for 'ratio' parameter.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "ratio", "some");
%!error<sampsizepwr: invalid value for 'ratio' parameter.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "ratio", 0.5);
%!error<sampsizepwr: invalid value for 'ratio' parameter.> ...
%! out = sampsizepwr ("r", 0.2, [], 0.8, 60, "ratio", [2, 1.3, 0.3]);
%!error<sampsizepwr: only one of either p1, power, or n must be missing.> ...
%! out = sampsizepwr ("z", [100, 5], [], [], 60);
%!error<sampsizepwr: only one of either p1, power, or n must be missing.> ...
%! out = sampsizepwr ("z", [100, 5], 110, [], []);
%!error<sampsizepwr: only one of either p1, power, or n must be missing.> ...
%! out = sampsizepwr ("z", [100, 5], [], 0.8, []);
%!error<sampsizepwr: only one of either p1, power, or n must be missing.> ...
%! out = sampsizepwr ("z", [100, 5], 110, 0.8, 60);
%!error<sampsizepwr: alternative hypothesis parameter must be numeric.> ...
%! out = sampsizepwr ("z", [100, 5], "mu", [], 60);
%!error<sampsizepwr: alternative hypothesis parameter out of range.> ...
%! out = sampsizepwr ("var", 5, -1, [], 60);
%!error<sampsizepwr: alternative hypothesis parameter out of range.> ...
%! out = sampsizepwr ("p", 0.8, 1.2, [], 60, "tail", "right");
%!error<sampsizepwr: alternative hypothesis parameter out of range.> ...
%! out = sampsizepwr ("r", 0.8, 1.2, [], 60);
%!error<sampsizepwr: alternative hypothesis parameter out of range.> ...
%! out = sampsizepwr ("r", 0.8, -1.2, [], 60);
%!error<sampsizepwr: invalid value for POWER.> ...
%! out = sampsizepwr ("z", [100, 5], 110, 1.2);
%!error<sampsizepwr: invalid value for POWER.> ...
%! out = sampsizepwr ("z", [100, 5], 110, 0);
%!error<sampsizepwr: Cannot compute N or P1 unless POWER> ...
%! out = sampsizepwr ("z", [100, 5], 110, 0.05, [], "alpha", 0.1);
%!error<sampsizepwr: input arguments size mismatch.> ...
%! out = sampsizepwr ("z", [100, 5], [], [0.8, 0.7], [60, 80, 100]);
%!error<sampsizepwr: Same value for null and alternative hypothesis.> ...
%! out = sampsizepwr ("t", [100, 5], 100, 0.8, []);
%!error<sampsizepwr: Invalid P1 for testing left tail.> ...
%! out = sampsizepwr ("t", [100, 5], 110, 0.8, [], "tail", "left");
%!error<sampsizepwr: Invalid P1 for testing right tail.> ...
%! out = sampsizepwr ("t", [100, 5], 90, 0.8, [], "tail", "right");

## Warning test
%!warning<sampsizepwr: approximate N.> ...
%! Napprox = sampsizepwr ("p", 0.2, 0.26, 0.6);
%!warning<sampsizepwr: approximate N.> ...
%! Napprox = sampsizepwr ("p", 0.30, 0.36, 0.8);

## Results validation
%!test
%! mu1 = sampsizepwr ("t", [100, 10], [], 0.8, 60);
%! assert (mu1, 103.67704316, 1e-8);
%!test
%! [N1,N2] = sampsizepwr ("t2", [100, 10], 110, 0.6, [], "ratio", 1.5);
%! assert (N1, 9);
%! assert (N2, 14);
%!test
%! nn = 1:250;
%! pwr = sampsizepwr ("p", 0.2, 0.26, [], nn);
%! pwr_out = [0, 0.0676, 0.0176, 0.0566, 0.0181, 0.0431, 0.0802, 0.0322];
%! assert (pwr([1:8]), pwr_out, 1e-4 * ones (1,8));
%! pwr_out = [0.59275, 0.6073, 0.62166, 0.6358, 0.6497, 0.6087, 0.6229, 0.6369];
%! assert (pwr([243:end]), pwr_out, 1e-4 * ones (1,8));
%!test
%! nout = sampsizepwr ("t", [100, 5], 102, 0.80);
%! assert (nout, 52);
%!test
%! power = sampsizepwr ("t", [20, 5], 25, [], 5, "Tail", "right");
%! assert (power, 0.5797373588621888, 1e-14);
%!test
%! nout = sampsizepwr ("t", [20, 5], 25, 0.99, [], "Tail", "right");
%! assert (nout, 18);
%!test
%! p1out = sampsizepwr ("t", [20, 5], [], 0.95, 10, "Tail", "right");
%! assert (p1out, 25.65317979360237, 1e-14);
%!test
%! pwr = sampsizepwr ("t2", [1.4, 0.2], 1.7, [], 5, "Ratio", 2);
%! assert (pwr, 0.716504004686586, 1e-14);
%!test
%! n = sampsizepwr ("t2", [1.4, 0.2], 1.7, 0.9, []);
%! assert (n, 11);
%!test
%! [n1, n2] = sampsizepwr ("t2", [1.4, 0.2], 1.7, 0.9, [], "Ratio", 2);
%! assert ([n1, n2], [8, 16]);
