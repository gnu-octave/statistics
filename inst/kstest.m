## Copyright (C) 2022-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{h} =} kstest (@var{x})
## @deftypefnx {statistics} {@var{h} =} kstest (@var{x}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{h}, @var{p}] =} kstest (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{ksstat}, @var{cv}] =} kstest (@dots{})
##
## Single sample Kolmogorov-Smirnov (K-S) goodness-of-fit hypothesis test.
##
## @code{@var{h} = kstest (@var{x})} performs a Kolmogorov-Smirnov (K-S) test to
## determine if a random sample @var{x} could have come from a standard normal
## distribution.  @var{h} indicates the results of the null hypothesis test.
##
## @itemize
## @item @var{h} = 0 => Do not reject the null hypothesis at the 5% significance
## @item @var{h} = 1 => Reject the null hypothesis at the 5% significance
## @end itemize
##
## @var{x} is a vector representing a random sample from some unknown
## distribution with a cumulative distribution function F(X).  Missing values
## declared as NaNs in @var{x} are ignored.
##
## @code{@var{h} = kstest (@var{x}, @var{name}, @var{value})} returns
## a test decision for a single-sample K-S test with additional options
## specified by one or more @var{Name}-@var{Value} pair arguments as shown
## below.
##
## @multitable @columnfractions 0.15 0.05 0.8
## @headitem Name @tab @tab Value
## @item @qcode{"alpha"} @tab @tab A numeric scalar between 0 and 1 specifying th
## the significance level.  Default is 0.05 for 5% significance.
##
## @item @qcode{"CDF"} @tab @tab The hypothesized CDF under the null hypothesis.
## It can be specified as a function handle of an existing cdf function, a
## character vector defining a probability distribution with default parameters,
## a probability distribution object, or a two-column matrix.  If not provided,
## the default is the standard normal, @math{N(0,1)}.  The one-sample
## Kolmogorov-Smirnov test is only valid for continuous cumulative distribution
## functions, and requires the CDF to be predetermined.  The result is not
## accurate if CDF is estimated from the data.
##
## @item @qcode{"tail"} @tab @tab A string indicating the type of test:
## @end multitable
## @multitable @columnfractions 0.2 0.15 0.05 0.5
## @item @tab @qcode{"unequal"} @tab @tab "F(X) not equal to CDF(X)" (two-sided)
## (Default)
##
## @item @tab @qcode{"larger"} @tab @tab "F(X) > CDF(X)" (one-sided)
##
## @item @tab @qcode{"smaller"} @tab @tab "F(X) < CDF(X)" (one-sided)
## @end multitable
##
## Let S(X) be the empirical c.d.f. estimated from the sample vector @var{x},
## F(X) be the corresponding true (but unknown) population c.d.f., and CDF be
## the known input c.d.f. specified under the null hypothesis.
## For @code{tail} = "unequal", "larger", and "smaller", the test statistics are
## max|S(X) - CDF(X)|, max[S(X) - CDF(X)], and max[CDF(X) - S(X)], respectively.
##
## @code{[@var{h}, @var{p}] = kstest (@dots{})} also returns the asymptotic
## p-value @var{p}.
##
## @code{[@var{h}, @var{p}, @var{ksstat}] = kstest (@dots{})} returns the K-S
## test statistic @var{ksstat} defined above for the test type indicated by the
## "tail" option
##
## In the matrix version of CDF, column 1 contains the x-axis data and column 2
## the corresponding y-axis c.d.f data.  Since the K-S test statistic will
## occur at one of the observations in @var{x}, the calculation is most
## efficient when CDF is only specified at the observations in @var{x}.  When
## column 1 of CDF represents x-axis points independent of @var{x}, CDF is
## linearly interpolated at the observations found in the vector @var{x}.  In
## this case, the interval along the x-axis (the column 1 spread of CDF) must
## span the observations in @var{x} for successful interpolation.
##
## The decision to reject the null hypothesis is based on comparing the p-value
## @var{p} with the "alpha" value, not by comparing the statistic @var{ksstat}
## with the critical value @var{cv}.  @var{cv} is computed separately using an
## approximate formula or by interpolation using Miller's approximation table.
## The formula and table cover the range 0.01 <= "alpha" <= 0.2 for two-sided
## tests and 0.005 <= "alpha" <= 0.1 for one-sided tests.  CV is returned as NaN
## if "alpha" is outside this range.  Since CV is approximate, a comparison of
## @var{ksstat} with @var{cv} may occasionally lead to a different conclusion
## than a comparison of @var{p} with "alpha".
##
## @seealso{kstest2, cdfplot}
## @end deftypefn

function [H, pValue, ksstat, cV] = kstest (x, varargin)

  ## Check input parameters
  if (nargin < 1)
    error ("kstest: too few inputs.");
  endif
  if (! isvector (x) || ! isreal (x))
    error ("kstest: X must be a vector of real numbers.");
  endif

  ## Add defaults
  alpha = 0.05;
  tail = "unequal";
  CDF = [];

  ## Parse extra parameters
  if (length (varargin) > 0 && mod (numel (varargin), 2) == 0)
    [~, prop] = parseparams (varargin);
    while (!isempty (prop))
      switch (lower (prop{1}))
        case "alpha"
          alpha = prop{2};
        case "tail"
          tail = prop{2};
        case "cdf"
          CDF = prop{2};
        otherwise
          error ("kstest: unknown option '%s'.", prop{1});
      endswitch
      prop = prop(3:end);
    endwhile
  elseif (mod (numel (varargin), 2) != 0)
    error ("kstest: optional parameters must be in name/value pairs.");
  endif

  ## Check for valid alpha and tail parameters
  if (! isnumeric (alpha) || isnan (alpha) || ! isscalar (alpha) ...
                          || alpha <= 0 || alpha >= 1)
    error ("kstest: alpha must be a numeric scalar in the range (0,1).");
  endif
  if (! isa (tail, 'char'))
    error ("kstest: tail argument must be a string.");
  elseif (sum (strcmpi (tail, {"unequal", "larger", "smaller"})) < 1)
    error ("kstest: tail value must be either 'both', right' or 'left'.");
  endif

  ## Remove NaNs, get sample size and compute empirical cdf
  x(isnan (x)) = [];
  n = length (x);
  [sampleCDF, x] = ecdf (x);

  ## Remove 1st element
  x = x(2:end);

  ## Check the hypothesized CDF specified under the null hypothesis.
  ## No CDF was provided (use default)
  if (isempty (CDF))
    xCDF = x;
    yCDF = normcdf (x, 0, 1);

  ## If CDF is a function handle
  elseif (isa (CDF, "function_handle"))
    xCDF = x;
    yCDF = feval (CDF, x);
    if (! isequal (size (xCDF), size (yCDF)))
      error ("kstest: invalid function handle.");
    endif

  ## If CDF is character vector
  elseif (isa (CDF, "char") && isvector (CDF))
    ## Check for supported distibutions
    PDO = makedist ();
    if (! any (strcmpi (PDO, CDF)))
      error ("kstest: '%s' is not a supported distribution.", CDF);
    endif
    pd = makedist (CDF);
    xCDF = x;
    yCDF = pd.cdf (x);

  ## If CDF is a probability distribution object
  elseif (isobject (CDF))
    PDO = makedist ();
    PDO = cellfun (@(x) sprintf ("%sDistribution", x), PDO, ...
                   "UniformOutput", false);

    if (! any (isa (CDF, PDO)))
      error ("kstest: 'CDF' must be a probability distribution object.");
    endif
    xCDF = x;
    yCDF = CDF.cdf (x);

  ## If CDF is numerical
  elseif (! isempty (CDF) && isnumeric (CDF))
    if (size (CDF, 2) != 2)
      error ("kstest: numerical CDF should have only 2 columns.");
    endif
    CDF(isnan (sum (CDF, 2)),:) = [];
    if (size (CDF, 1) == 0)
      error ("kstest: numerical CDF should have at least one row.");
    endif
    ## Sort numerical CDF
    [xCDF, i] =  sort (CDF(:,1));
    yCDF = CDF(i,2);
    ## Check that numerical CDF is incrementally sorted
    ydiff = diff (yCDF);
    if (any (ydiff < 0))
      error("kstest: non-incrementing numerical CDF.");
    endif
    ## Remove duplicates. Check for consistency
    rd = find (diff (xCDF) == 0);
    if (! isempty (rd))
      if (! all (ydiff(rd) == 0))
        error ("kstest: wrong duplicates in numerical CDF.");
      endif
      xCDF(rd) = [];
      yCDF(rd) = [];
    endif

  ## Invalid value parsed as CDF optinonal argument
  else
    error ("kstest: invalid value parsed as CDF optinonal argument.");
  endif

  ## Check if CDF is specified at the observations in X and assign 2nd column
  ## of numerical CDF to null CDF
  if (isequal (x, xCDF))
    nCDF = yCDF;
  ## Otherwise interpolate the numerical CDF to assign values to the null CDF
  else
    ## Check that 1st column range bounds the observations in X
    if (x(1) < xCDF(1) || x(end) > xCDF(end))
      error ("kstest: wrong span in CDF.");
    endif
    nCDF  =  interp1 (xCDF, yCDF, x);
  endif

  ## Calculate the suitable KS statistic according to tail
  switch (tail)
    case "unequal"    # 2-sided test: T = max|S(x) - CDF(x)|.
      delta1    =  sampleCDF(1:end - 1) - nCDF;
      delta2    =  sampleCDF(2:end) - nCDF;
      deltaCDF  =  abs ([delta1; delta2]);
    case "smaller"    # 1-sided test: T = max[CDF(x) - S(x)].
      delta1    =  nCDF - sampleCDF(1:end - 1);
      delta2    =  nCDF - sampleCDF(2:end);
      deltaCDF  =  [delta1; delta2];
    case "larger"     # 1-sided test: T = max[S(x) - CDF(x)].
      delta1    =  sampleCDF(1:end - 1) - nCDF;
      delta2    =  sampleCDF(2:end) - nCDF;
      deltaCDF  =  [delta1; delta2];
  endswitch
  ksstat = max (deltaCDF);

  ## Compute the asymptotic P-value approximation
  if (strcmpi (tail, "unequal"))    # 2-sided test
    s = n * ksstat ^ 2;
    ## For d values that are in the far tail of the distribution (i.e.
    ## p-values > .999), the following lines will speed up the computation
    ## significantly, and provide accuracy up to 7 digits.
    if ((s > 7.24) || ((s > 3.76) && (n > 99)))
      pValue = 2 * exp (-(2.000071 + 0.331 / sqrt (n) + 1.409 / n) * s);
    else
      ## Express d as d = (k-h)/n, where k is a +ve integer and 0 < h < 1.
      k = ceil (ksstat * n);
      h = k - ksstat * n;
      m = 2 * k - 1;
      ## Create the H matrix, according to Marsaglia et al.
      if (m > 1)
        c = 1 ./ gamma ((1:m)' + 1);
        r = zeros (1,m);
        r(1) = 1;
        r(2) = 1;
        T = toeplitz (c, r);
        T(:,1) = T(:,1) - (h .^ (1:m)') ./ gamma ((1:m)' + 1);
        T(m,:) = fliplr (T(:,1)');
        T(m,1) = (1 - 2 * h ^ m + max (0, 2 * h - 1) ^m) / gamma (m+1);
      else
        T = (1 - 2 * h ^ m + max (0, 2 * h - 1) ^ m) / gamma (m+1);
      endif
      ## Scaling before raising the matrix to a power
      if (! isscalar (T))
        lmax = max (eig (T));
        T = (T ./ lmax) ^ n;
      else
        lmax = 1;
      endif
      pValue = 1 - exp (gammaln (n+1) + n * log (lmax) - n * log (n)) * T(k,k);
    endif
  else                              # 1-sided test
    t = n * ksstat;
    k = ceil (t):n;
    pValue = sum (exp (log (t) - n * log (n) + gammaln (n + 1) ...
             - gammaln (k + 1) - gammaln (n - k + 1) + k .* log (k - t) ...
             + (n - k - 1) .* log (t + n - k)));
  endif

  ## Return hypothesis test
  H = (pValue < alpha);

  ## Calculate critical Value (cV) if requested
  if (nargout > 3)
    ## The critical value table used below is expressed in reference to a
    ## 1-sided significance level. Hence alpha is halved for a two-sided test.
    if (strcmpi (tail, "unequal"))  # 2-sided test
        alpha1 = alpha / 2;
    else                            # 1-sided test
        alpha1 = alpha;
    endif

    if ((alpha1 >= 0.005) && (alpha1 <= 0.10))
      ## If the sample size 'n' is greater than 20, use Miller's approximation
      ## Otherwise interpolate into his 'exact' table.

      if (n <= 20)                  # Small sample exact values.
        % Exact K-S test critical values based on Miller's approximation.
        a1    = [0.00500, 0.01000, 0.02500, 0.05000, 0.10000]';

        exact = [0.99500, 0.99000, 0.97500, 0.95000, 0.90000; ...
                 0.92929, 0.90000, 0.84189, 0.77639, 0.68377; ...
                 0.82900, 0.78456, 0.70760, 0.63604, 0.56481; ...
                 0.73424, 0.68887, 0.62394, 0.56522, 0.49265; ...
                 0.66853, 0.62718, 0.56328, 0.50945, 0.44698; ...
                 0.61661, 0.57741, 0.51926, 0.46799, 0.41037; ...
                 0.57581, 0.53844, 0.48342, 0.43607, 0.38148; ...
                 0.54179, 0.50654, 0.45427, 0.40962, 0.35831; ...
                 0.51332, 0.47960, 0.43001, 0.38746, 0.33910; ...
                 0.48893, 0.45662, 0.40925, 0.36866, 0.32260; ...
                 0.46770, 0.43670, 0.39122, 0.35242, 0.30829; ...
                 0.44905, 0.41918, 0.37543, 0.33815, 0.29577; ...
                 0.43247, 0.40362, 0.36143, 0.32549, 0.28470; ...
                 0.41762, 0.38970, 0.34890, 0.31417, 0.27481; ...
                 0.40420, 0.37713, 0.33760, 0.30397, 0.26588; ...
                 0.39201, 0.36571, 0.32733, 0.29472, 0.25778; ...
                 0.38086, 0.35528, 0.31796, 0.28627, 0.25039; ...
                 0.37062, 0.34569, 0.30936, 0.27851, 0.24360; ...
                 0.36117, 0.33685, 0.30143, 0.27136, 0.23735; ...
                 0.35241, 0.32866, 0.29408, 0.26473, 0.23156];

        cV  =  spline (a1 , exact(n,:)' , alpha1);

      else                          # Large sample approximate values.

        A = 0.09037 * (-log10 (alpha1)) .^ 1.5 + 0.01515 * ...
            log10 (alpha1) .^ 2 - 0.08467 * alpha1 - 0.11143;
        asymptoticStat = sqrt (-0.5 * log (alpha1) ./ n);
        cV = asymptoticStat - 0.16693 ./ n - A ./ n .^ 1.5;
        cV = min (cV, 1 - alpha1);
      endif

    else
      cV = NaN;
    endif
  endif
endfunction

%!demo
%! ## Use the stock return data set to test the null hypothesis that the data
%! ## come from a standard normal distribution against the alternative
%! ## hypothesis that the population CDF of the data is larger that the
%! ## standard normal CDF.
%!
%! load stockreturns;
%! x = stocks(:,2);
%! [h, p, k, c] = kstest (x, "Tail", "larger")
%!
%! ## Compute the empirical CDF and plot against the standard normal CDF
%! [f, x_values] = ecdf (x);
%! h1 = plot (x_values, f);
%! hold on;
%! h2 = plot (x_values, normcdf (x_values), 'r--');
%! set (h1, "LineWidth", 2);
%! set (h2, "LineWidth", 2);
%! legend ([h1, h2], "Empirical CDF", "Standard Normal CDF", ...
%!         "Location", "southeast");
%! title ("Empirical CDF of stock return data against standard normal CDF")

## Test input
%!error<kstest: too few inputs.> kstest ()
%!error<kstest: X must be a vector of real numbers.> kstest (ones (2, 4))
%!error<kstest: X must be a vector of real numbers.> kstest ([2, 3, 5, 3+3i])
%!error<kstest: unknown option 'opt'.> kstest ([2, 3, 4, 5, 6], "opt", 0.51)
%!error<kstest: optional parameters must be in name/value pairs.> ...
%! kstest ([2, 3, 4, 5, 6], "tail")
%!error<kstest: alpha must be a numeric scalar in the range> ...
%! kstest ([2,3,4,5,6],"alpha", [0.05, 0.05])
%!error<kstest: alpha must be a numeric scalar in the range> ...
%! kstest ([2, 3, 4, 5, 6], "alpha", NaN)
%!error<kstest: tail argument must be a string.> ...
%! kstest ([2, 3, 4, 5, 6], "tail", 0)
%!error<kstest: tail value must be either 'both', right' or 'left'.> ...
%! kstest ([2,3,4,5,6], "tail", "whatever")
%!error<kstest: invalid function handle.> ...
%! kstest ([1, 2, 3, 4, 5], "CDF", @(x) repmat (x, 2, 3))
%!error<kstest: 'somedist' is not a supported distribution.> ...
%! kstest ([1, 2, 3, 4, 5], "CDF", "somedist")
%!error<kstest: 'CDF' must be a probability distribution object.> ...
%! kstest ([1, 2, 3, 4, 5], "CDF", cvpartition (5))
%!error<kstest: numerical CDF should have only 2 columns.> ...
%! kstest ([2, 3, 4, 5, 6], "alpha", 0.05, "CDF", [2, 3, 4; 1, 3, 4; 1, 2, 1])
%!error<kstest: numerical CDF should have at least one row.> ...
%! kstest ([2, 3, 4, 5, 6], "alpha", 0.05, "CDF", nan (5, 2))
%!error<kstest: non-incrementing numerical CDF.> ...
%! kstest ([2, 3, 4, 5, 6], "CDF", [2, 3; 1, 4; 3, 2])
%!error<kstest: wrong duplicates in numerical CDF.> ...
%! kstest ([2, 3, 4, 5, 6], "CDF", [2, 3; 2, 4; 3, 5])
%!error<kstest: invalid value parsed as CDF optinonal argument.> ...
%! kstest ([2, 3, 4, 5, 6], "CDF", {1, 2, 3, 4, 5})

## Test results
%!test
%! load examgrades
%! [h, p] = kstest (grades(:,1));
%! assert (h, true);
%! assert (p, 7.58603305206105e-107, 1e-14);
%!test
%! load examgrades
%! [h, p] = kstest (grades(:,1), "CDF", @(x) normcdf(x, 75, 10));
%! assert (h, false);
%! assert (p, 0.5612, 1e-4);
%!test
%! load examgrades
%! x = grades(:,1);
%! test_cdf = makedist ("tlocationscale", "mu", 75, "sigma", 10, "nu", 1);
%! [h, p] = kstest (x, "alpha", 0.01, "CDF", test_cdf);
%! assert (h, true);
%! assert (p, 0.0021, 1e-4);
%!test
%! load stockreturns
%! x = stocks(:,3);
%! [h,p,k,c] = kstest (x, "Tail", "larger");
%! assert (h, true);
%! assert (p, 5.085438806199252e-05, 1e-14);
%! assert (k, 0.2197, 1e-4);
%! assert (c, 0.1207, 1e-4);
