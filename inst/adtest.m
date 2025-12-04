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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{h} =} adtest (@var{x})
## @deftypefnx {statistics} {@var{h} =} adtest (@var{x}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{h}, @var{pval}] =} adtest (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{adstat}, @var{cv}] =} adtest (@dots{})
##
## Anderson-Darling goodness-of-fit hypothesis test.
##
## @code{@var{h} = adtest (@var{x})} returns a test decision for the null
## hypothesis that the data in vector @var{x} is from a population with a normal
## distribution, using the Anderson-Darling test.  The alternative hypothesis is
## that x is not from a population with a normal distribution. The result
## @var{h} is 1 if the test rejects the null hypothesis at the 5% significance
## level, or 0 otherwise.
##
## @code{@var{h} = adtest (@var{x}, @var{Name}, @var{Value})} returns a test
## decision for the Anderson-Darling test with additional options specified by
## one or more Name-Value pair arguments.  For example, you can specify a null
## distribution other than normal, or select an alternative method for
## calculating the p-value, such as a Monte Carlo simulation.
##
## The following parameters can be parsed as Name-Value pair arguments.
##
## @multitable @columnfractions 0.25 0.75
## @headitem Name @tab Description
## @item "Distribution" @tab The distribution being tested for.  It tests
## whether @var{x} could have come from the specified distribution.  There are
## two choices available for parsing distribution parameters:
## @end multitable
##
## @itemize
## @item
## One of the following char strings: "norm", "exp", "ev", "logn", "weibull",
## for defining either the 'normal', 'exponential', 'extreme value', lognormal,
## or 'Weibull' distribution family, respectively.  In this case, @var{x} is
## tested against a composite hypothesis for the specified distribution family
## and the required distribution parameters are estimated from the data in
## @var{x}.  The default is "norm".
##
## @item
## A cell array defining a distribution in which the first cell contains a char
## string with the distribution name, as mentioned above, and the consecutive
## cells containing all specified parameters of the null distribution.  In this
## case, @var{x} is tested against a simple hypothesis.
## @end itemize
##
## @multitable @columnfractions 0.25 0.75
## @item "Alpha" @tab Significance level alpha for the test.  Any scalar
## numeric value between 0 and 1.  The default is 0.05 corresponding to the 5%
## significance level.
##
## @item "MCTol" @tab Monte-Carlo standard error for the p-value,
## @var{pval}, value. which must be a positive scalar value.  In this case, an
## approximation for the p-value is computed directly, using Monte-Carlo
## simulations.
##
## @item "Asymptotic" @tab Method for calculating the p-value of the
## Anderson-Darling test, which can be either true or false logical value.  If
## you specify 'true', adtest estimates the p-value using the limiting
## distribution of the Anderson-Darling test statistic.  If you specify 'false',
## adtest calculates the p-value based on an analytical formula.  For sample
## sizes greater than 120, the limiting distribution estimate is likely to be
## more accurate than the small sample size approximation method.
## @end multitable
##
## @itemize
## @item
## If you specify a distribution family with unknown parameters for the
## distribution Name-Value pair (i.e. composite distribution hypothesis test),
## the "Asymptotic" option must be false.
## @item
##
## If you use MCTol to calculate the p-value using a Monte Carlo simulation,
## the "Asymptotic" option must be false.
## @end itemize
##
## @code{[@var{h}, @var{pval}] = adtest (@dots{})} also returns the p-value,
## @var{pval}, of the Anderson-Darling test, using any of the input arguments
## from the previous syntaxes.
##
## @code{[@var{h}, @var{pval}, @var{adstat}, @var{cv}] = adtest (@dots{})} also
## returns the test statistic, @var{adstat}, and the critical value, @var{cv},
## for the Anderson-Darling test.
##
## The Anderson-Darling test statistic belongs to the family of Quadratic
## Empirical Distribution Function statistics, which are based on the weighted
## sum of the difference @math{[Fn(x)-F(x)]^2} over the ordered sample values
## @math{X1 < X2 < ... < Xn}, where @math{F} is the hypothesized continuous
## distribution and @math{Fn} is the empirical CDF based on the data sample with
## @math{n} sample points.
##
## @seealso{kstest}
## @end deftypefn

function [H, pVal, ADStat, CV] = adtest (x, varargin)

  ## Check for valid input data
  if (nargin < 1)
    print_usage;
  endif
  ## Ensure the sample data is a real vector.
  if (! isvector (x) || ! isreal(x))
    error ("adtest: X must be a vector of real numbers.");
  endif

  ## Add defaults
  distribution = "norm";
  isCompositeTest = true;
  alpha = 0.05;
  MCTol = [];
  asymptotic = false;

  ## Parse arguments and check if parameter/value pairs are valid
  i = 1;
  while (i <= length (varargin))
    switch lower (varargin{i})
      case "distribution"
        i = i + 1;
        distribution = varargin{i};
        ## Check for char string or cell array
        ## If distribution is a char string, then X is tested against a
        ## composite hypothesis for the specified distribution family and
        ## distribution parameters are estimated from the sample data.
        valid_dists = {"norm", "exp", "ev", "logn", "weibull"};
        if (ischar (distribution))
          if (! any (strcmpi (distribution, valid_dists)))
            error ("adtest: invalid distribution family in char string.");
          endif
        ## If distribution is a cell array, then X is tested against a
        ## simple hypothesis for the specified distribution and
        ## distribution parameters given in the cell array.
        elseif (iscell (distribution))
          if (! any (strcmpi (distribution(1), valid_dists)))
            error ("adtest: invalid distribution family in cell array.");
          endif
          ## Check for valid distribution parameter in cell array
          err_msg = "adtest: invalid distribution parameters in cell array.";
          if (strcmpi (distribution(1), "norm"))
            if (numel (distribution) != 3)
              error (err_msg);
            endif
            z = normcdf (x, distribution{2}, distribution{3});
          elseif (strcmpi (distribution(1), "exp"))
            if (numel (distribution) != 2)
              error (err_msg);
            endif
            z = expcdf (x, distribution{2});
          elseif (strcmpi (distribution(1), "ev"))
            if (numel (distribution) != 3)
              error (err_msg);
            endif
            z = evcdf (x, distribution{2}, distribution{3});
          elseif (strcmpi (distribution(1), "logn"))
            if (numel (distribution) != 3)
              error (err_msg);
            endif
            z = logncdf (x, distribution{2}, distribution{3});
          elseif (strcmpi (distribution(1), "weibull"))
            if (numel (distribution) != 3)
              error (err_msg);
            endif
            z = wblcdf (x, distribution{2}, distribution{3});
          endif
          isCompositeTest = false;
        else
          error ("adtest: invalid distribution option.");
        endif
      case "alpha"
        i = i + 1;
        alpha = varargin{i};
        ## Check for valid alpha
        if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                    alpha <= 0 || alpha >= 1)
          error ("adtest: invalid value for alpha.");
        endif
      case "mctol"
        i = i + 1;
        MCTol = varargin{i};
        ## Check Monte Carlo tolerance is a numeric scalar greater than 0
        if (! isempty (MCTol) && (! isscalar (MCTol) || MCTol <= 0))
          error ("adtest: invalid Monte Carlo Tolerance.");
        endif
      case "asymptotic"
        i = i + 1;
        asymptotic = varargin{i};
        ## Check that it is either true or false
        if (! isbool (asymptotic))
          error ("adtest: asymptotic option must be boolean.");
        endif
      otherwise
        error ("adtest: invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  ## Check conflicts with asymptotic option
  if (asymptotic && isCompositeTest)
    error (strcat ("adtest: asymptotic option is not valid", ...
                   " for the composite distribution test."));
  elseif (asymptotic && ! isempty (MCTol))
    error (strcat ("adtest: asymptotic option is not valid", ...
                   " for the Monte Carlo simulation test."));
  endif

  ## Remove missing values.
  x = x(! isnan (x));
  ## Compute sample size n.
  n = length (x);

  ## For composite tests
  if (isCompositeTest)
    ## Abort for sample size less than 4
    if (n < 4)
      error ("adtest: not enough data for composite testing.");
    endif

    ## If data follow a lognormal distribution, log(x) is normally distributed
    if (strcmpi (distribution, "logn"))
      x = log (x);
      distribution = "norm";
    ## If data follow a Weibull distribution, log(x) has a type I extreme-value
    ## distribution
    elseif (strcmpi (distribution, "weibull"))
      x = log (x);
      distribution = "ev";
    endif

    ## Compute ADStat
    switch distribution
      case "norm"
        ## Check for complex numbers due to log (x)
        if (any (! isreal (x)))
          ## Data is not compatible with logn distribution test
          warning ("adtest: bad data for lognormal distribution.");
          ADStat = NaN;
        else
          z = normcdf (x, mean (x), std (x));
          ADStat = ComputeADStat (z, n);
        endif
      case "exp"
        z = expcdf (x, mean (x));
        ADStat = ComputeADStat (z, n);
      case "ev"
        ## Check for complex numbers due to log (x)
        if (any (! isreal (x)))
          ## Data is not compatible with Weibull distribution test
          warning ("adtest: bad data for Weibull distribution.");
          ADStat = NaN;
        else
          params = evfit (x);
          z = evcdf (x, params (1), params (2));
          ADStat = ComputeADStat (z, n);
        endif
    endswitch

    ## Compute p-value and critical values without Monte Carlo simulation
    if (isempty (MCTol))
      alphas = [0.0005, 0.0010, 0.0015, 0.0020, 0.0050, 0.0100, 0.0250, ...
                0.0500, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.3500, ...
                0.4000, 0.4500, 0.5000, 0.5500, 0.6000, 0.6500, 0.7000, ...
                0.7500, 0.8000, 0.8500, 0.9000, 0.9500, 0.9900];
      switch distribution
        case "norm"
          CVs = computeCriticalValues_norm (n);
        case "exp"
          CVs = computeCriticalValues_exp (n);
        case "ev"
          CVs = computeCriticalValues_ev (n);
      endswitch

      ## 1-D interpolation into the tabulated results
      pp = pchip (log (alphas), CVs);
      CV = ppval (pp, log (alpha));

      ## If alpha is not within the lookup table, throw a warning
      ## Hypothesis result is computed by comparing the p-value with
      ## alpha, rather than CV with ADStat
      if alpha < alphas(1)
        CV = CVs(1);
        warning ("adtest: alpha not within the lookup table.");
      elseif  alpha > alphas(end)
        CV = CVs(end);
        warning ("adtest: alpha not within the lookup table.");
      endif

      if (ADStat > CVs(1))
        ## P value is smaller than smallest tabulated value
        warning (strcat (["adtest: out of range min p-value:"], ...
                         sprintf (" %g", alphas(1))));
        pVal = alphas(1);
      elseif (ADStat < CVs(end))
        ## P value is larger than largest tabulated value
        warning (strcat (["adtest: out of range max p-value:"], ...
                         sprintf (" %g", alphas(end))));
        pVal = alphas(end);
      elseif (isnan (ADStat))
        ## Handle certain cases of negative data (ADStat == NaN)
        pVal = 0;
      else
        ## Find p-value by inverse interpolation
        i = find (ADStat > CVs, 1, "first");
        logPVal = fzero (@(x)ppval(pp,x) - ADStat, log(alphas([i-1,i])));
        pVal = exp (logPVal);
      endif

    ## Compute p-value and critical values without Monte Carlo simulation
    else
        [CV, pVal] = adtestMC (ADStat, n, alpha, distribution, mctol);
    endif

    ## Calculate H
    if (isnan (ADStat))
      H = true;
    else
      if (isempty (MCTol))
        if (alpha < alphas(1) || alpha > alphas(end))
          H = (pVal < alpha);
        else
          H = (ADStat > CV);
        endif
      else
          H = (ADStat > CV);
      endif
    endif

  ## For simple tests
  else
    ## Compute the Anderson-Darling statistic
    ADStat = ComputeADStat (z, n);

    ## Compute p-value and critical values without Monte Carlo simulation
    if (isempty (MCTol))
      alphas = [0.0005, 0.0010, 0.0015, 0.0020, 0.0050, 0.0100, 0.0250, ...
                0.0500, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.3500, ...
                0.4000, 0.4500, 0.5000, 0.5500, 0.6000, 0.6500, 0.7000, ...
                0.7500, 0.8000, 0.8500, 0.9000, 0.9500, 0.9900];
      if (asymptotic)
        if (n <= 120)
          warning ("adtest: asymptotic distribution with small sample size.");
        endif
        pVal = 1 - ADInf (ADStat);
        ## For extra output arguments
        if (nargout > 3)
          ## Make sure alpha is within the lookup table
          validateAlpha (alpha, alphas);
          ## Find critical values
          critVals = findAsymptoticDistributionCriticalValues;
          i = find (alphas > alpha, 1, "first");
          startVal = critVals(i-1);
          CV = fzero(@(ad)1-ADInf(ad)-alpha, startVal);
        endif
      else
        if (n == 1)
          pVal = 1 - sqrt (1 - 4 * exp (-1 - ADStat));
        else
          if (n < 4)
            warning ("adtest: small sample size.");
          endif
          pVal = 1 - ADn (n, ADStat);
        endif
        ## For extra output arguments
        if (nargout > 3)
          ## Make sure alpha is within the lookup table
          validateAlpha (alpha, alphas);
          ## Find critical values
          [CVs, sampleSizes] = findCriticalValues;
          [OneOverSampleSizes, LogAlphas] = meshgrid (1 ./ sampleSizes, ...
                                                      log (alphas));
          CV = interp2 (OneOverSampleSizes, LogAlphas, CVs', 1./n, log (alpha));
        endif
      endif

    ## Compute p-value and critical values with Monte Carlo simulation
    else
      [CV, pVal] = adtestMC (ADStat, n, alpha, "unif", mctol);
    endif

    ## Calculate H
    H  =  (pVal < alpha);

  endif

endfunction

## Compute Anderson-Darling Statistic
function ADStat = ComputeADStat (z, n)
  ## Sort the data and compute the statistic
  z = reshape (z, n, 1);
  z = sort (z);
  w = 2 * (1:n) - 1;
  ADStat = - w * (log (z)+ log (1-z(end:-1:1))) / n - n;
endfunction

## Anderson-Darling distribution Pr(An<z) for arbitrary n.
function p = ADn (n, ad)
  ## Check for invalid ADStat
  if (any (ad < 0))
    error ("adtest: invalid Anderson Darling statistic.");
  endif
  x = zeros (size (ad));
  adl = ad(ad < 2);
  x(ad < 2) = adl .^ (-1 / 2) .* exp (-1.2337 ./ adl) .* ...
              (2.00012 + (0.247105 - (0.0649821 - (0.0347962 - ...
              (0.0116720 - 0.00168691 * adl) .* adl) .* adl) .* adl) .* adl);
  adh = ad(ad >= 2);
  x(ad >= 2) = exp (-exp (1.0776 - (2.30695 - (0.43424 - (0.082433 - ...
               (0.008056 - 0.0003146 .* adh) .* adh) .* adh) .* adh) .* adh));
  ## Compute error function defined between 0 and 1
  if (any (x < 0 | x > 1))
    error ("adtest: invalid values for error function.");
  endif
  e = zeros (size (x));
  c = 0.01265 + 0.1757/n;

  ## Define function by intervals using 3 fixed functions g1, g2, g3.
  xc1 = x(x < c) / c;
  g1 = sqrt (xc1) .* (1 - xc1) .* (49 * xc1 - 102);
  e(x < c) = (0.0037 / n ^ 3 + 0.00078 / n ^ 2 + 0.00006 / n) * g1;

  xc2 = (x(x >= c & x < 0.8) - c) ./ (0.8 - c);
  g2 = -0.00022633 + (6.54034 - (14.6538 - (14.458 - (8.259 -...
                      1.91864 .* xc2) .* xc2) .* xc2) .* xc2) .* xc2;
  e(x >= c & x < 0.8) = (0.04213 / n + 0.01365 / n ^ 2) * g2;

  xc3 = x(x >= 0.8);
  e(x >= 0.8) = 1 / n * (-130.2137 + (745.2337 - (1705.091 - (1950.646 -...
                        (1116.360 - 255.7844 .* ...
                         xc3) .* xc3) .* xc3) .* xc3) .* xc3) .* xc3;
  p = x + e;
endfunction

## Evaluate the Anderson-Darling limit distribution
function ad = ADInf (z)
  ## Distribution is invalid for negative values
  if (z < 0)
    error ("adtest: invalid X for asymptotic distribution.");
  endif
  ## Due to floating point precision:
  ## Return 0 below a certain threshold
  if (z < 0.02)
    ad = 0;
    return;
  ## Return 1 above the following threshold
  elseif (z >= 32.4)
    ad = 1;
    return;
  endif
  n = 1:500;
  K = 1/z*[1, ((4*n + 1).*cumprod((1/2 - n)./n))];
  ADTerms = arrayfun(@(j)ADf(z,j),0:500);
  ad = ADTerms*K';
endfunction

## Series expansion for f(z,j) called by ADInf
function f = ADf(z,j)
  ## Compute t=tj=(4j+1)^2*pi^2/(8z)
  t = (4 * j + 1) ^ 2 * 1.233700550136170 / z;
  ## First 2 terms in recursive series
  ## c0=pi*exp(-t)/(sqrt(2t))
  ## c1=pi*sqrt(pi/2)*erfc(sqrt(t))
  c0 = 2.221441469079183 * exp (-t) / sqrt (t);
  c1 = 3.937402486430604 * erfc (sqrt (t));
  r = z / 8;
  f = c0 + c1 * r;
  ## Evaluate the recursion
  for n = 2:500
    c = 1 / (n - 1) * ((n - 3 / 2 - t) * c1 + t * c0);
    r = r * (z / 8) * (1 / n);
    fn = f + c * r;
    c0 = c1;
    c1 = c;
    if (f == fn)
      return;
    endif
    f = fn;
  endfor
endfunction

## An improved version of the Petitt method for the composite normal case.
function CVs = computeCriticalValues_norm (n)
  CVs = [1.5649,  1.4407,  1.3699,  1.3187,  1.1556,  1.0339,  0.8733, ...
         0.7519,  0.6308,  0.5598,  0.5092,  0.4694,  0.4366,  0.4084, ...
         0.3835,  0.3611,  0.3405,  0.3212,  0.3029,  0.2852,  0.2679, ...
         0.2506,  0.2330,  0.2144,  0.1935,  0.1673,  0.1296] + ...
       [-0.9362, -0.9029, -0.8906, -0.8865, -0.8375, -0.7835, -0.6746, ...
        -0.5835, -0.4775, -0.4094, -0.3679, -0.3327, -0.3099, -0.2969, ...
        -0.2795, -0.2623, -0.2464, -0.2325, -0.2164, -0.1994, -0.1784, ...
        -0.1569, -0.1377, -0.1201, -0.0989, -0.0800, -0.0598] ./ n + ...
       [-8.3249, -6.6022, -5.6461, -4.9685, -3.2208, -2.1647, -1.2460, ...
        -0.7803, -0.4627, -0.3672, -0.2833, -0.2349, -0.1442, -0.0229, ...
         0.0377,  0.0817,  0.1150,  0.1583,  0.1801,  0.1887,  0.1695, ...
         0.1513,  0.1533,  0.1724,  0.2027,  0.3158,  0.6431] ./ n ^ 2;
endfunction

## An improved version of the Petitt method for the composite exponential case.
function CVs = computeCriticalValues_exp (n)
  CVs = [3.2371,  2.9303,  2.7541,  2.6307,  2.2454,  1.9621,  1.5928, ...
         1.3223,  1.0621,  0.9153,  0.8134,  0.7355,  0.6725,  0.6194, ...
         0.5734,  0.5326,  0.4957,  0.4617,  0.4301,  0.4001,  0.3712, ...
         0.3428,  0.3144,  0.2849,  0.2527,  0.2131,  0.1581] + ...
        [1.6146,  0.8716,  0.4715,  0.2066, -0.4682, -0.7691, -0.7388, ...
        -0.5758, -0.4036, -0.3142, -0.2564, -0.2152, -0.1845, -0.1607, ...
        -0.1409, -0.1239, -0.1084, -0.0942, -0.0807, -0.0674, -0.0537, ...
        -0.0401, -0.0261, -0.0116,  0.0047,  0.0275,  0.0780] ./ n;
endfunction

## An improved version of the Petitt method for the composite extreme value case
function CVs = computeCriticalValues_ev(n)
  CVs = [1.6473,  1.5095,  1.4301,  1.3742,  1.1974,  1.0667,  0.8961, ...
         0.7683,  0.6416,  0.5680,  0.5156,  0.4744,  0.4405,  0.4115, ...
         0.3858,  0.3626,  0.3415,  0.3217,  0.3029,  0.2848,  0.2672, ...
         0.2496,  0.2315,  0.2124,  0.1909,  0.1633,  0.1223] + ...
       [-0.7097, -0.5934, -0.5328, -0.4930, -0.3708, -0.2973, -0.2075, ...
        -0.1449, -0.0892, -0.0619, -0.0442, -0.0302, -0.0196, -0.0112, ...
        -0.0039,  0.0024,  0.0074,  0.0122,  0.0167,  0.0207,  0.0245, ...
         0.0282,  0.0323,  0.0371,  0.0436,  0.0549,  0.0813] ./ n .^ (1 / 2);
endfunction

## Find rows of critical values at relevant significance levels
function [CVs, sampleSizes] = findCriticalValues
  CVs = [7.2943, 6.6014, 6.1962, 5.9088, 4.9940, 4.3033, 3.3946, 2.7142, ...
         2.0470, 1.6682, 1.4079, 1.2130, 1.0596, 0.9353, 0.8326, 0.7465, ...
         0.6740, 0.6126, 0.5606, 0.5170, 0.4806, 0.4508, 0.4271, 0.4091, ...
         NaN, NaN, NaN; ...  %n=1
         7.6624, 6.4955, 5.9916, 5.6682, 4.7338, 4.0740, 3.2247, 2.5920, ...
         1.9774, 1.6368, 1.4078, 1.2329, 1.0974, 0.9873, 0.8947, 0.8150, ...
         0.7448, 0.6820, 0.6251, 0.5727, 0.5240, 0.4779, 0.4337, 0.3903, ...
         0.3462, 0.3030, 0.2558; ...  %n=2
         7.2278, 6.3094, 5.8569, 5.5557, 4.6578, 4.0111, 3.1763, 2.5581, ...
         1.9620, 1.6314, 1.4079, 1.2390, 1.1065, 0.9979, 0.9060, 0.8264, ...
         0.7560, 0.6928, 0.6350, 0.5816, 0.5314, 0.4835, 0.4371, 0.3907, ...
         0.3424, 0.2885, 0.2255; ...  %n=3
         7.0518, 6.2208, 5.7904, 5.4993, 4.6187, 3.9788, 3.1518, 2.5414, ...
         1.9545, 1.6288, 1.4080, 1.2416, 1.1104, 1.0025, 0.9110, 0.8315, ...
         0.7611, 0.6977, 0.6397, 0.5859, 0.5352, 0.4868, 0.4395, 0.3920, ...
         0.3421, 0.2845, 0.2146; ...  %n=4
         6.9550, 6.1688, 5.7507, 5.4653, 4.5949, 3.9591, 3.1370, 2.5314, ...
         1.9501, 1.6272, 1.4080, 1.2430, 1.1126, 1.0051, 0.9138, 0.8344, ...
         0.7640, 0.7005, 0.6424, 0.5884, 0.5375, 0.4888, 0.4411, 0.3930, ...
         0.3424, 0.2833, 0.2097; ...  %n=5
         6.8935, 6.1345, 5.7242, 5.4426, 4.5789, 3.9459, 3.1271, 2.5248, ...
         1.9472, 1.6262, 1.4081, 1.2439, 1.1140, 1.0067, 0.9156, 0.8362, ...
         0.7658, 0.7023, 0.6441, 0.5901, 0.5391, 0.4901, 0.4422, 0.3938, ...
         0.3427, 0.2828, 0.2071; ...  %n=6
         6.8509, 6.1102, 5.7053, 5.4264, 4.5674, 3.9364, 3.1201, 2.5201, ...
         1.9451, 1.6255, 1.4081, 1.2445, 1.1149, 1.0079, 0.9168, 0.8375, ...
         0.7671, 0.7036, 0.6454, 0.5912, 0.5401, 0.4911, 0.4430, 0.3944, ...
         0.3430, 0.2826, 0.2056; ...  %n=7
         6.8196, 6.0920, 5.6912, 5.4142, 4.5588, 3.9293, 3.1148, 2.5166, ...
         1.9436, 1.6249, 1.4081, 1.2450, 1.1156, 1.0087, 0.9177, 0.8384, ...
         0.7681, 0.7045, 0.6463, 0.5921, 0.5409, 0.4918, 0.4436, 0.3949, ...
         0.3433, 0.2825, 0.2046; ...  %n=8
         6.7486, 6.0500, 5.6582, 5.3856, 4.5384, 3.9124, 3.1024, 2.5084, ...
         1.9400, 1.6237, 1.4081, 1.2460, 1.1171, 1.0106, 0.9197, 0.8406, ...
         0.7702, 0.7066, 0.6483, 0.5941, 0.5428, 0.4935, 0.4451, 0.3961, ...
         0.3441, 0.2826, 0.2029; ...  %n=12
         6.7140, 6.0292, 5.6417, 5.3713, 4.5281, 3.9040, 3.0962, 2.5044, ...
         1.9382, 1.6230, 1.4081, 1.2465, 1.1179, 1.0115, 0.9207, 0.8416, ...
         0.7712, 0.7077, 0.6493, 0.5950, 0.5437, 0.4944, 0.4459, 0.3968, ...
         0.3445, 0.2827, 0.2023; ...  %n=16
         6.6801, 6.0084, 5.6252, 5.3569, 4.5178, 3.8955, 3.0900, 2.5003, ...
         1.9365, 1.6224, 1.4081, 1.2470, 1.1186, 1.0123, 0.9217, 0.8426, ...
         0.7723, 0.7087, 0.6503, 0.5960, 0.5446, 0.4952, 0.4466, 0.3974, ...
         0.3450, 0.2829, 0.2019; ...  %n=24
         6.6468, 5.9877, 5.6087, 5.3425, 4.5075, 3.8869, 3.0837, 2.4963, ...
         1.9347, 1.6218, 1.4082, 1.2474, 1.1193, 1.0132, 0.9226, 0.8436, ...
         0.7732, 0.7097, 0.6513, 0.5969, 0.5455, 0.4960, 0.4474, 0.3980, ...
         0.3455, 0.2832, 0.2016; ...  %n=48
         6.6634, 5.9980, 5.6169, 5.3497, 4.5127, 3.8912, 3.0868, 2.4983, ...
         1.9356, 1.6221, 1.4081, 1.2472, 1.1190, 1.0128, 0.9222, 0.8431, ...
         0.7728, 0.7092, 0.6508, 0.5965, 0.5451, 0.4956, 0.4470, 0.3977, ...
         0.3453, 0.2830, 0.2017; ...  %n=32
         6.6385, 5.9825, 5.6046, 5.3389, 4.5049, 3.8848, 3.0822, 2.4953, ...
         1.9343, 1.6217, 1.4082, 1.2475, 1.1195, 1.0134, 0.9228, 0.8438, ...
         0.7735, 0.7099, 0.6516, 0.5972, 0.5458, 0.4962, 0.4476, 0.3982, ...
         0.3456, 0.2833, 0.2016; ...  %n=64
         6.6318, 5.9783, 5.6012, 5.3360, 4.5028, 3.8830, 3.0809, 2.4944, ...
         1.9339, 1.6215, 1.4082, 1.2476, 1.1197, 1.0136, 0.9230, 0.8440, ...
         0.7737, 0.7101, 0.6517, 0.5974, 0.5459, 0.4964, 0.4477, 0.3984, ...
         0.3457, 0.2833, 0.2015; ...  %n=88
         6.6297, 5.9770, 5.6001, 5.3350, 4.5021, 3.8825, 3.0805, 2.4942, ...
         1.9338, 1.6215, 1.4082, 1.2476, 1.1197, 1.0136, 0.9231, 0.8441, ...
         0.7738, 0.7102, 0.6518, 0.5974, 0.5460, 0.4965, 0.4478, 0.3984, ...
         0.3458, 0.2834, 0.2015; ...  %n=100
         6.6262, 5.9748, 5.5984, 5.3335, 4.5010, 3.8816, 3.0798, 2.4937, ...
         1.9336, 1.6214, 1.4082, 1.2477, 1.1198, 1.0137, 0.9232, 0.8442, ...
         0.7739, 0.7103, 0.6519, 0.5975, 0.5461, 0.4966, 0.4479, 0.3985, ...
         0.3458, 0.2834, 0.2015; ...  %n=128
         6.6201, 5.9709, 5.5953, 5.3308, 4.4990, 3.8800, 3.0787, 2.4930, ...
         1.9333, 1.6213, 1.4082, 1.2478, 1.1199, 1.0139, 0.9234, 0.8443, ...
         0.7740, 0.7105, 0.6521, 0.5977, 0.5463, 0.4967, 0.4480, 0.3986, ...
         0.3459, 0.2834, 0.2015; ...  %n=256
         6.6127, 5.9694, 5.5955, 5.3314, 4.4982, 3.8781, 3.0775, 2.4924, ...
         1.9330, 1.6212, 1.4082, 1.2479, 1.1201, 1.0140, 0.9235, 0.8445, ...
         0.7742, 0.7106, 0.6523, 0.5979, 0.5464, 0.4969, 0.4481, 0.3987, ...
         0.3460, 0.2835, 0.2015];     %n=Inf

  sampleSizes = [1 2 3 4 5 6 7 8 12 16 24 32 48 64 88 100 128 256 Inf];
endfunction
% ------------------------------------------
function critVals = findAsymptoticDistributionCriticalValues
  critVals = [6.6127034546551, 5.9694013422151, 5.5954643397078, ...
              5.3313658857909, 4.4981996466091, 3.8781250216054, ...
              3.0774641787107, 2.4923671600494, 1.9329578327416, ...
              1.6212385363175, 1.4081977005506, 1.2478596347253, ...
              1.1200136586965, 1.0140004020016, 0.9235137094902, ...
              0.8445069178452, 0.7742142410993, 0.7106405935247, ...
              0.6522701010084, 0.5978828157471, 0.5464229310982, ...
              0.4968804113119, 0.4481425895777, 0.3987228486242, ...
              0.3460480234939, 0.2835161264344, 0.2014922164166]; %n=Inf
endfunction

## Make sure alpha is within the lookup table
function validateAlpha (alpha, alphas)
  if (alpha < alphas(1) || alpha > alphas(end))
    error (strcat (["adtest: out of range invalid alpha -"], ...
                   sprintf (" lower limit: %g", alphas(1)),...
                   sprintf (" upper limit: %g", alphas(end))));
  endif
endfunction

## Simulated critical values and p-values for Anderson-Darling test
function [crit, p] = adtestMC (ADStat, n, alpha, distribution, MCTol)
  ## Initial values vartol = mctol^2;
  crit = 0;
  p = 0;
  mcRepsTot = 0;
  mcRepsMin = 1000;
  ## Monte Carlo loop
  while true
    mcRepsOld = mcRepsTot;
    mcReps = ceil(mcRepsMin - mcRepsOld);
    ADstatMC = zeros(mcReps,1);
    ## Switch to selected distribution
    switch distribution
      case "norm"
        mu0 = 0;
        sigma0 = 1;
        for rep = 1:length (ADstatMC)
          x = normrnd (mu0, sigma0, n, 1);
          xCDF = sort (x);
          nullCDF = normcdf (xCDF, mean (x), std (x));
          w = 2 * (1:n) - 1 ;
          ADstatMC(rep) = - w * (log (nullCDF) + ...
                                 log (1 - nullCDF(end:-1:1))) / n - n;
        endfor
      case "exp"
        beta0 = 1;
        for rep = 1:length (ADstatMC)
          x = exprnd (beta0, n, 1);
          xCDF = sort (x);
          nullCDF = expcdf (xCDF, mean (x));
          w = 2 * (1:n) - 1 ;
          ADstatMC(rep) = - w * (log (nullCDF) + ...
                                 log (1 - nullCDF(end:-1:1))) / n - n;
        endfor
      case "ev"
        mu0 = 0;
        sigma0 = 1;
        for rep = 1:length (ADstatMC)
          x = evrnd (mu0, sigma0, n, 1);
          pHat = evfit (x);
          xCDF = sort (x);
          nullCDF = evcdf (xCDF, pHat(1), pHat(2));
          w = 2 * (1:n) - 1 ;
          ADstatMC(rep) = - w * (log (nullCDF) + ...
                                 log (1 - nullCDF(end:-1:1))) / n - n;
        endfor
      case "unif"
        for rep = 1:length(ADstatMC)
          z = sort (rand (n, 1));
          w = 2 * (1:n) - 1 ;
          ADstatMC(rep) = - w * (log (z) + ...
                                 log (1 - z(end:-1:1))) / n - n;
        endfor
    endswitch

    critMC = prctile (ADstatMC, 100 * (1 - alpha));
    pMC = sum (ADstatMC > ADStat) ./ mcReps;

    mcRepsTot = mcRepsOld + mcReps;
    crit = (mcRepsOld * crit + mcReps * critMC) / mcRepsTot;
    p = (mcRepsOld * p + mcReps * pMC) / mcRepsTot;

    ## Compute a std err for p, with lower bound (1/N)*(1-1/N)/N when p==0.
    sepsq = max (p * (1 - p) / mcRepsTot, 1 / mcRepsTot ^ 2);
    if (sepsq < MCTol ^ 2)
      break
    endif
    ## Based on the current estimate, find the number of trials needed to
    ## make the MC std err less than the specified tolerance.
    mcRepsMin = 1.2 * (mcRepsTot * sepsq) / (MCTol ^ 2);
  endwhile
endfunction

## Test input validation
%!error<Invalid call to adtest.  Correct usage is:> adtest ();
%!error<adtest: X must be a vector of real numbers.> adtest (ones (20,2));
%!error<adtest: X must be a vector of real numbers.> adtest ([1+i,0-3i]);
%!error<adtest: invalid distribution family in char string.> ...
%! adtest (ones (20,1), "Distribution", "normal");
%!error<adtest: invalid distribution family in cell array.> ...
%! adtest (rand (20,1), "Distribution", {"normal", 5, 3});
%!error<adtest: invalid distribution parameters in cell array.> ...
%! adtest (rand (20,1), "Distribution", {"norm", 5});
%!error<adtest: invalid distribution parameters in cell array.> ...
%! adtest (rand (20,1), "Distribution", {"exp", 5, 4});
%!error<adtest: invalid distribution parameters in cell array.> ...
%! adtest (rand (20,1), "Distribution", {"ev", 5});
%!error<adtest: invalid distribution parameters in cell array.> ...
%! adtest (rand (20,1), "Distribution", {"logn", 5, 3, 2});
%!error<adtest: invalid distribution parameters in cell array.> ...
%! adtest (rand (20,1), "Distribution", {"Weibull", 5});
%!error<adtest: invalid distribution option.> ...
%! adtest (rand (20,1), "Distribution", 35);
%!error<adtest: invalid Name argument.> ...
%! adtest (rand (20,1), "Name", "norm");
%!error<adtest: invalid Name argument.> ...
%! adtest (rand (20,1), "Name", {"norm", 75, 10});
%!error<adtest: asymptotic option is not valid for the composite> ...
%! adtest (rand (20,1), "Distribution", "norm", "Asymptotic", true);
%!error<adtest: asymptotic option is not valid for the composite> ...
%! adtest (rand (20,1), "MCTol", 0.001, "Asymptotic", true);
%!error<adtest: asymptotic option is not valid for the Monte Carlo> ...
%! adtest (rand (20,1), "Distribution", {"norm", 5, 3}, "MCTol", 0.001, ...
%!         "Asymptotic", true);
%!error<adtest: out of range invalid alpha - lower limit: 0.0005 upper> ...
%! [h, pval, ADstat, CV] = adtest (ones (20,1), "Distribution", {"norm",5,3},...
%!                                 "Alpha", 0.000000001);
%!error<adtest: out of range invalid alpha - lower limit: 0.0005 upper> ...
%! [h, pval, ADstat, CV] = adtest (ones (20,1), "Distribution", {"norm",5,3},...
%!                                 "Alpha", 0.999999999);
%!error<adtest: not enough data for composite testing.> ...
%! adtest (10);

## Test warnings
%!warning<adtest: out of range min p-value:> ...
%! randn ("seed", 34);
%! adtest (ones (20,1), "Alpha", 0.000001);
%!warning<adtest: alpha not within the lookup table.> ...
%! randn ("seed", 34);
%! adtest (normrnd(0,1,100,1), "Alpha", 0.99999);
%!warning<adtest: alpha not within the lookup table.> ...
%! randn ("seed", 34);
%! adtest (normrnd(0,1,100,1), "Alpha", 0.00001);

## Test results
%!test
%! load examgrades
%! x = grades(:,1);
%! [h, pval, adstat, cv] = adtest (x);
%! assert (h, false);
%! assert (pval, 0.1854, 1e-4);
%! assert (adstat, 0.5194, 1e-4);
%! assert (cv, 0.7470, 1e-4);
%!test
%! load examgrades
%! x = grades(:,1);
%! [h, pval, adstat, cv] = adtest (x, "Distribution", "ev");
%! assert (h, false);
%! assert (pval, 0.071363, 1e-6);
%!test
%! load examgrades
%! x = grades(:,1);
%! [h, pval, adstat, cv] = adtest (x, "Distribution", {"norm", 75, 10});
%! assert (h, false);
%! assert (pval, 0.4687, 1e-4);
