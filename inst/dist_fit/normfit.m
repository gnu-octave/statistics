## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{muhat} =} normfit (@var{x})
## @deftypefnx {statistics} {[@var{muhat}, @var{sigmahat}] =} normfit (@var{x})
## @deftypefnx {statistics} {[@var{muhat}, @var{sigmahat}, @var{muci}] =} normfit (@var{x})
## @deftypefnx {statistics} {[@var{muhat}, @var{sigmahat}, @var{muci}, @var{sigmaci}] =} normfit (@var{x})
## @deftypefnx {statistics} {[@dots{}] =} normfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} normfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} normfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} normfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for the normal distribution.
##
## @code{[@var{muhat}, @var{sigmahat}] = normfit (@var{x})} estimates the
## parameters of the normal distribution given the data in @var{x}.  @var{muhat}
## is an estimate of the mean, and @var{sigmahat} is an estimate of the standard
## deviation.
##
## @code{[@var{muhat}, @var{sigmahat}, @var{muci}, @var{sigmaci}] = normfit
## (@var{x})} returns the 95% confidence intervals for the mean and standard
## deviation estimates in the arrays @var{muci} and @var{sigmaci}, respectively.
##
## @itemize
## @item
## @var{x} can be a vector or a matrix.  When @var{x} is a matrix, the parameter
## estimates and their confidence intervals are computed for each column.  In
## this case, @code{normfit} supports only 2 input arguments, @var{x} and
## @var{alpha}.  Optional arguments @var{censor}, @var{freq}, and @var{options}
## can be used only when @var{x} is a vector.
##
## @item
## @var{alpha} is a scalar value in the range @math{(0,1)} specifying the
## confidence level for the confidence intervals calculated as
## @math{100x(1 – alpha)%}.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @item
## @var{censor} is a logical vector of the same length as @var{x} specifying
## whether each value in @var{x} is right-censored or not.  1 indicates
## observations that are right-censored and 0 indicates observations that are
## fully observed.  With censoring, @var{muhat} and @var{sigmahat} are the
## maximum likelihood estimates (MLEs).  If empty, the default is an array of
## 0s, meaning that all observations are fully observed.
##
## @item
## @var{freq} is a vector of the same length as @var{x} and it typically
## contains non-negative integer counts of the corresponding elements in
## @var{x}.  If empty, the default is an array of 1s, meaning one observation
## per element of @var{x}.  To obtain the weighted MLEs for a data set with
## censoring, specify weights of observations, normalized to the number of
## observations in @var{x}.  However, when there is no censored data (default),
## the returned estimate for standard deviation is not exactly the WMLE.  To
## compute the weighted MLE, multiply the value returned in @var{sigmahat} by
## @code{(SUM (@var{freq}) - 1) / SUM (@var{freq})}.  This correction is needed
## because @code{normfit} normally computes @var{sigmahat} using an unbiased
## variance estimator when there is no censored data.  When there is censoring
## in the data, the correction is not needed, since @code{normfit} does not use
## the unbiased variance estimator in that case.
##
## @item
## @var{options} is a structure with the control parameters for
## @code{fminsearch} which is used internally to compute MLEs for censored data.
## By default, it uses the following options:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.MaxFunEvals = 400}
## @item @qcode{@var{options}.MaxIter = 200}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
## @end itemize
##
## @seealso{normcdf, norminv, normpdf, normrnd, normlike, normstat}
## @end deftypefn

function [muhat, sigmahat, muci, sigmaci] = normfit (x, alpha, censor, freq, options)

  ## Check for valid number of input arguments
  narginchk (1, 5);

  ## Check X for being a vector or a matrix
  if (ndims (x) != 2)
    error ("normfit: X must not be a multi-dimensional array.");
  endif
  if (! isvector (x))
    if (nargin < 3)
      [n, ncols] = size (x);
    else
      error ("normfit: matrix data acceptable only under 2-arg syntax.");
    endif
  else
    n = numel (x);
    ncols = 1;
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("normfit: Wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = 0;
  elseif (! isequal (size (x), size (censor)))
    error ("normfit: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = 1;
  elseif (isequal (size (x), size (freq)))
    n = sum (freq);
    is_zero = find (freq == 0);
    if (numel (is_zero) > 0)
      x(is_zero) = [];
      if (numel (censor) == numel (freq))
        censor(is_zero) = [];
      endif
      freq(is_zero) = [];
    end
  else
    error ("normfit: X and FREQ vector mismatch.");
  endif

  ## Check options structure or add defaults
  if (nargin > 4 && ! isempty (options))
    if (! isstruct (options) || ! isfield (options, "Display") ||
        ! isfield (options, "MaxFunEvals") || ! isfield (options, "MaxIter")
                                           || ! isfield (options, "TolX"))
      error (strcat (["normfit: 'options' 5th argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  else
    options.Display = "off";
    options.MaxFunEvals = 400;
    options.MaxIter = 200;
    options.TolX = 1e-6;
  endif

  ## Get number of censored and uncensored elements
  n_censored = sum(freq.*censor); % a scalar in all cases
  n_uncensored = n - n_censored; % a scalar in all cases

  ## Compute total sum in X
  totalsum = sum(freq.*x);

  ## Check cases that cannot make a fit.
  ## 1. Handle Infs and NaNs
  if (! isfinite (totalsum))
    muhat = totalsum;
    sigmahat = NaN ("like", x);
    muci = NaN (2, 1, "like", x);
    sigmaci = NaN (2, 1, "like", x);
    return
  endif

  ## 2. All observations are censored or empty data
  if (n == 0 || n_uncensored == 0)
    muhat = NaN (1, ncols,'like',x);
    sigmahat = NaN (1, ncols,'like',x);
    muci = NaN (2, ncols,'like',x);
    sigmaci = NaN (2, ncols,'like',x);
    return
  endif

  ## 3. No censored values, compute parameter estimates explicitly.
  if (n_censored == 0)
    muhat = totalsum ./ n;
    if (n > 1)
      if numel(muhat) == 1  # X is a vector
        xc = x - muhat;
      else                  # X is a matrix
        xc = x - repmat (muhat, [n, 1]);
      endif
      sigmahat = sqrt (sum (conj (xc) .* xc .* freq) ./ (n - 1));
    else
      sigmahat = zeros (1, ncols, "like", x);
    endif

    if (nargout > 2)
      if (n > 1)
        paramhat = [muhat; sigmahat];
        ci = norm_ci (paramhat, [], alpha, x, [], freq);
        muci = ci(:,:,1);
        sigmaci = ci(:,:,2);
      else
        muci = [-Inf; Inf] * ones (1, ncols, "like", x);
        sigmaci = [0; Inf] * ones (1, ncols, "like", x);
      endif
    endif
    return
  endif

  ## 4. Αll uncensored observations equal and greater than all the
  ## censored observations
  x_uncensored = x(censor == 0);
  range_x_uncensored = range (x_uncensored);
  if (range_x_uncensored < realmin (class (x)))
    if (x_uncensored(1) == max (x))
      muhat = x_uncensored(1);
      sigmahat = zeros ('like',x);
      if (n_uncensored > 1)
        muci = [muhat; muhat];
        sigmaci = zeros (2, 1, "like", x);
      else
        muci = cast ([-Inf; Inf], "like", x);
        sigmaci = cast ([0; Inf], "like", x);
      endif
      return
    endif
  endif

  ## Get an initial estimate for parameters using the "least squares" method
  if (range_x_uncensored > 0)
    if (numel (freq) == numel (x))
      [p,q] = ecdf (x, "censoring", censor, "frequency", freq);
    else
      [p,q] = ecdf (x, "censoring", censor);
    endif
    pmid = (p(1:(end-1)) + p(2:end)) / 2;
    linefit = polyfit (-sqrt (2) * erfcinv (2 * pmid), q(2:end), 1);
    paramhat = linefit ([2 1]);
  else  # only one uncensored element in X
    paramhat = [x_uncensored(1) 1];
  endif

  ## Optimize the parameters as doubles, regardless of input data type
  paramhat = cast (paramhat, "double");

  ## Search for parameter that minimize the negative log likelihood function
  [paramhat, ~, err, output] = fminsearch ...
                    (@(ph) norm_nlogl (ph, x, censor, freq), paramhat, options);

  ## Handle errors
  if (err == 0)
    if (output.funcCount >= options.MaxFunEvals)
      warning ("normfit: maximum number of function evaluations are exceeded.");
    elseif (output.iterations >= options.MaxIter)
      warning ("normfit: maximum number of iterations are exceeded.");
    endif
  elseif (err < 0)
    error ("normfit: NoSolution.");
  endif

  ## Make sure the outputs match the input data type
  muhat = cast (paramhat(1), "like", x);
  sigmahat = cast (paramhat(2), "like", x);

  if (nargout > 2)
    paramhat = paramhat(:);
    if (numel (freq) == numel (x))
      [~, avar] = normlike (paramhat, x, censor, freq);
    else
      [~, avar] = normlike (paramhat, x, censor);
    endif
    ci = norm_ci (paramhat, avar, alpha, x, censor, freq);
    muci = ci(:,:,1);
    sigmaci = ci(:,:,2);
  endif

endfunction

## Negative log-likelihood function and gradient for normal distribution.
function [nlogL, avar] = norm_nlogl (params, x, censor, freq)
  ## Get mu and sigma values
  mu = params(1);
  sigma = params(2);
  ## Compute the individual log-likelihood terms.  Force a log(0)==-Inf for
  ## data from extreme right tail, instead of getting exp(Inf-Inf)==NaN.
  z = (x - mu) ./ sigma;
  L = -0.5 .* z .^ 2 - log (sqrt (2 .* pi) .* sigma);
  if (any (censor))
    censored = censor == 1;
    z_censor = z(censored);
    S_censor = 0.5 * erfc (z_censor / sqrt (2));
    L(censored) = log (S_censor);
  endif
  ## Neg-log-like is the sum of the individual contributions
  nlogL = -sum (freq .* L);
  ## Compute the negative hessian at the parameter values.
  ## Invert to get the observed information matrix.
  if (nargout == 2)
    dL11 = -ones (size (z), class (z));
    dL12 = -2 .* z;
    dL22 = 1 - 3 .* z .^ 2;
    if (any (censor))
      dlogScen = exp (-0.5 .* z_censor .^ 2) ./ (sqrt (2 * pi) .* S_censor);
      d2logScen = dlogScen .* (dlogScen - z_censor);
      dL11(censored) = -d2logScen;
      dL12(censored) = -dlogScen - z_censor .* d2logScen;
      dL22(censored) = -z_censor .* (2 .* dlogScen + z_censor .* d2logScen);
    endif
    nH11 = -sum (freq .* dL11);
    nH12 = -sum (freq .* dL12);
    nH22 = -sum (freq .* dL22);
    avar =  (sigma .^ 2) * [nH22, -nH12; -nH12, nH11] / ...
                           (nH11 * nH22 - nH12 * nH12);
  endif
endfunction

## Confidence intervals for normal distribution
function ci = norm_ci (paramhat, cv, alpha, x, censor, freq)
  ## Check for missing input arguments
  if (nargin < 6 || isempty (freq))
    freq = ones (size (x));
  endif
  if (nargin < 5 || isempty (censor))
    censor = false (size (x));
  endif
  if (isvector (paramhat))
    paramhat = paramhat(:);
  endif
  muhat = paramhat(1,:);
  sigmahat = paramhat(2,:);
  ## Get number of elements
  if (isempty (freq) || isequal (freq, 1))
    if (isvector (x))
      n = length (x);
    else
      n = size (x, 1);
    endif
  else
    n = sum (freq);
  endif
  ## Get number of censored and uncensored elements
  n_censored = sum (freq .* censor);
  n_uncensored = n - n_censored;
  ## Just in case
  if (any (censor) && (n == 0 || n_uncensored == 0 || ! isfinite (paramhat(1))))
    ## X is a vector
    muci = NaN(2,1);
    sigmaci = NaN(2,1);
    ci = cast (cat (3, muci, sigmaci), "like", x);
    return
  endif
  ## Get confidence intervals for each parameter
  if ((isempty (censor) || ! any (censor(:))) && ! isequal (cv,zeros(2,2)))
    ## Use exact formulas
    tcrit = tinv ([alpha/2, 1-alpha/2], n-1);
    muci = [muhat+tcrit(1)*sigmahat/sqrt(n); muhat+tcrit(2)*sigmahat/sqrt(n)];
    chi2crit = chi2inv ([alpha/2, 1-alpha/2], n-1);
    sigmaci = [sigmahat*sqrt((n-1)./chi2crit(2)); ...
               sigmahat*sqrt((n-1)./chi2crit(1))];
  else
    ## Use normal approximation
    probs = [alpha/2; 1-alpha/2];
    se = sqrt (diag (cv))';
    z = norminv (probs);

    ## Compute the CI for mu using a normal distribution for muhat.
    muci = muhat + se(1) .* z;
    ## Compute the CI for sigma using a normal approximation for
    ## log(sigmahat), and transform back to the original scale.
    logsigci = log (sigmahat) + (se(2) ./ sigmahat) .* z;
    sigmaci = exp (logsigci);
  endif
  ## Return as a single array
  ci = cat (3, muci, sigmaci);
endfunction

%!demo
%! ## Sample 3 populations from 3 different normal distibutions
%! randn ("seed", 1);    # for reproducibility
%! r1 = normrnd (2, 5, 5000, 1);
%! randn ("seed", 2);    # for reproducibility
%! r2 = normrnd (5, 2, 5000, 1);
%! randn ("seed", 3);    # for reproducibility
%! r3 = normrnd (9, 4, 5000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 15, 0.4);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! hold on
%!
%! ## Estimate their mu and sigma parameters
%! [muhat, sigmahat] = normfit (r);
%!
%! ## Plot their estimated PDFs
%! x = [min(r(:)):max(r(:))];
%! y = normpdf (x, muhat(1), sigmahat(1));
%! plot (x, y, "-pr");
%! y = normpdf (x, muhat(2), sigmahat(2));
%! plot (x, y, "-sg");
%! y = normpdf (x, muhat(3), sigmahat(3));
%! plot (x, y, "-^c");
%! ylim ([0, 0.5])
%! xlim ([-20, 20])
%! hold off
%! legend ({"Normalized HIST of sample 1 with mu=2, σ=5", ...
%!          "Normalized HIST of sample 2 with mu=5, σ=2", ...
%!          "Normalized HIST of sample 3 with mu=9, σ=4", ...
%!          sprintf("PDF for sample 1 with estimated mu=%0.2f and σ=%0.2f", ...
%!                  muhat(1), sigmahat(1)), ...
%!          sprintf("PDF for sample 2 with estimated mu=%0.2f and σ=%0.2f", ...
%!                  muhat(2), sigmahat(2)), ...
%!          sprintf("PDF for sample 3 with estimated mu=%0.2f and σ=%0.2f", ...
%!                  muhat(3), sigmahat(3))}, "location", "northwest")
%! title ("Three population samples from different normal distibutions")
%! hold off

## Test output
%!test
%! load lightbulb
%! idx = find (lightbulb(:,2) == 0);
%! censoring = lightbulb(idx,3) == 1;
%! [muHat, sigmaHat] = normfit (lightbulb(idx,1), [], censoring);
%! assert (muHat, 9496.595867378575, 1e-12);
%! assert (sigmaHat, 3064.021012796456, 1e-12);
%!test
%! x = normrnd (3, 5, [1000, 1]);
%! [muHat, sigmaHat, muCI, sigmaCI] = normfit (x, 0.01);
%! assert (muCI(1) < 3);
%! assert (muCI(2) > 3);
%! assert (sigmaCI(1) < 5);
%! assert (sigmaCI(2) > 5);

## Test input validation
%!error<normfit: X must not be a multi-dimensional array.> ...
%! normfit (ones (3,3,3))
%!error<normfit: matrix data acceptable only under 2-arg syntax.> ...
%! normfit (ones (20,3), [], zeros (20,1))
%!error<normfit: Wrong value for ALPHA.> ...
%! normfit (ones (20,1), 0)
%!error<normfit: Wrong value for ALPHA.> ...
%! normfit (ones (20,1), -0.3)
%!error<normfit: Wrong value for ALPHA.> ...
%! normfit (ones (20,1), 1.2)
%!error<normfit: Wrong value for ALPHA.> ...
%! normfit (ones (20,1), [0.05 0.1])
%!error<normfit: Wrong value for ALPHA.> ...
%! normfit (ones (20,1), 0.02+i)
%!error<normfit: X and CENSOR vector mismatch.> ...
%! normfit (ones (20,1), [], zeros(15,1))
%!error<normfit: X and FREQ vector mismatch.> ...
%! normfit (ones (20,1), [], zeros(20,1), ones(25,1))
%!error<normfit: > normfit (ones (20,1), [], zeros(20,1), ones(20,1), "options")
