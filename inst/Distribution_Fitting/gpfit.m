## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{paramhat} =} gpfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gpfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gpfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gpfit (@var{x}, @var{alpha}, @var{options})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gpfit (@var{x}, @var{alpha}, @var{options}, @var{freq})
##
## Estimate parameters and confidence intervals for the generalized Pareto
## distribution.
##
## @code{@var{paramhat} = gpfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the generalized Pareto distribution given the
## data in @var{x}.  @qcode{@var{paramhat}(1)} is the shape parameter, @var{k},
## and @qcode{@var{paramhat}(2)} is the scale parameter, @var{sigma}.
##
## @code{gpfit} does not estimate the location parameter @var{theta} and assumes
## it to be zero, so @var{x} must not contain negative values.  To fit data with
## a known nonzero @var{theta}, subtract it from @var{x} before calling
## @code{gpfit}; the estimates of @var{k} and @var{sigma} are unchanged by the
## shift.
##
## @code{[@var{paramhat}, @var{paramci}] = gpfit (@var{x})} returns the 95%
## confidence intervals for the estimated parameters @var{k} and @var{sigma} as
## a @math{2}-by-@math{2} matrix whose first row holds the lower bounds and
## whose second row holds the upper bounds.
##
## @code{[@dots{}] = gpfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = gpfit (@var{x}, @var{alpha}, @var{options})}
## specifies control parameters for the iterative algorithm used to compute ML
## estimates with the @code{fminsearch} function.  @var{options} is a structure
## with the following fields and their default values:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.MaxFunEvals = 400}
## @item @qcode{@var{options}.MaxIter = 200}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
##
## @code{[@dots{}] = gpfit (@var{x}, @var{alpha}, @var{options}, @var{freq})}
## accepts a vector of the same size as @var{x} giving the number of times each
## element of @var{x} was observed.  This fourth argument is an Octave
## extension; MATLAB's @code{gpfit} takes three inputs at most.
##
## When the shape parameter falls below @math{-1} the likelihood is unbounded:
## the density at the upper endpoint of the support diverges as that endpoint
## closes onto the largest observation, so no maximum likelihood estimate
## exists and whatever is returned is an arbitrary point on that ridge.
## @code{gpfit} warns in this case and returns @code{NaN} confidence intervals.
## The estimate it does return always keeps every observation strictly inside
## the fitted support, since the likelihood is infinite outside it.  This is a
## deliberate deviation: MATLAB has been measured returning parameters for such
## data under which the largest observation has zero density and its own
## @code{gplike} returns @code{Inf}.
##
## When @qcode{@var{k} = 0} and @qcode{@var{theta} = 0}, the Generalized Pareto
## is equivalent to the exponential distribution.  When @qcode{@var{k} > 0} and
## @code{@var{theta} = @var{k} / @var{k}} the Generalized Pareto is equivalent
## to the Pareto distribution.  The mean of the Generalized Pareto is not finite
## when @qcode{@var{k} >= 1} and the variance is not finite when
## @qcode{@var{k} >= 1/2}.  When @qcode{@var{k} >= 0}, the Generalized Pareto
## has positive density for @qcode{@var{x} > @var{theta}}, or, when
## @qcode{@var{theta} < 0}, for
## @qcode{0 <= (@var{x} - @var{theta}) / @var{sigma} <= -1 / @var{k}}.
##
## Further information about the generalized Pareto distribution can be found at
## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
##
## @seealso{gpcdf, gpinv, gppdf, gprnd, gplike, gpstat}
## @end deftypefn

function [paramhat, paramci] = gpfit (x, alpha, options, freq)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("gpfit: function called with too few input arguments.");
  endif

  ## Check X for being a vector
  if (isempty (x))
    paramhat = nan (1, 2, class (x));
    paramci = nan (2, 2, class (x));
    return
  elseif (! isvector (x) || ! isreal (x))
    error ("gpfit: X must be a vector of real values.");
  endif

  ## The location parameter is assumed to be zero, so no observation may fall
  ## below it.  Data with a known nonzero location is shifted by the caller.
  if (any (x < 0))
    error ("gpfit: X must not contain negative values.");
  endif

  ## Parse ALPHA argument or add default
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("gpfit: wrong value for ALPHA.");
  endif

  ## Parse FREQ argument or add default
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("gpfit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("gpfit: FREQ must not contain negative values.");
  endif

  ## Expand frequency vector (if necessary)
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  ## Get options structure or add defaults
  if (nargin < 3 || isempty (options))
    options.Display = 'off';
    options.MaxFunEvals = 400;
    options.MaxIter = 200;
    options.TolX = 1e-6;
  else
    if (! isstruct (options) || ! isfield (options, 'Display') ||
        ! isfield (options, 'MaxFunEvals') || ! isfield (options, 'MaxIter')
                                           || ! isfield (options, 'TolX'))
      error (strcat ("gpfit: 'options' 3rd argument must be a", ...
                     " structure with 'Display', 'MaxFunEvals',", ...
                     " 'MaxIter', and 'TolX' fields present."));
    endif
  endif

  ## Non-finite data is not removed: it propagates into the estimates, as it
  ## does in MATLAB and in every iterative fitter of this package.  Dropping
  ## missing observations belongs to the wrappers a user hands raw data to,
  ## FITDIST and MLE, not to the estimator.

  ## Get sample size, max and range of X
  x_max = max (x);
  x_size = length (x);
  x_range = range (x);

  ## Check for appropriate sample size or all observations being equal
  if (x_size == 0)
    paramhat = NaN (1, 2);
    paramci = NaN (2, 2);
    warning ("gpfit: X contains no data.");
    return
  elseif (x_range < realmin (class (x)))
    paramhat = cast ([NaN, 0], class (x));
    paramci = [paramhat; paramhat];
    warning ("gpfit: X contains constant data.");
    return
  endif

  ## Make an initial guess
  x_mean = mean (x);
  x_var = var (x);
  k0 = -0.5 .* (x_mean .^ 2 ./ x_var - 1);
  s0 = 0.5 .* x_mean .* (x_mean .^ 2 ./ x_var + 1);
  ## If initial guess fails, start with an exponential fit
  if (k0 < 0 && (x_max >= -s0 / k0))
    k0 = 0;
    s0 = x_mean;
  endif
  paramhat = [k0, log(s0)];

  ## Maximize the log-likelihood with respect to shape and log_scale.
  f = @(paramhat) negloglike (paramhat, x);
  [paramhat, ~, err, output] = fminsearch (f, paramhat, options);
  paramhat(2) = exp (paramhat(2));

  ## Check output of fminsearch and produce warnings or errors if applicable
  if (err == 0)
    if (output.funcCount >= options.MaxFunEvals)
      warning ("gpfit: reached evaluation limit.");
    else
      warning ("gpfit: reached iteration limit.");
    endif
  elseif (err < 0)
    error ("gpfit: no solution.");
  endif

  ## Check if converged to boundaries
  if ((paramhat(1) < 0) && (x_max > -paramhat(2)/paramhat(1) - options.TolX))
    warning (strcat ("gpfit: the fitted upper bound of the support has", ...
                     " closed onto the largest observation, a boundary of", ...
                     " the parameter space, so the estimates are unreliable", ...
                     " and no confidence intervals are computed."));
    reachedBnd = true;
  elseif (paramhat(1) <= -1 / 2)
    warning (strcat ("gpfit: the shape parameter has converged to", ...
                     " K <= -1/2, where the maximum likelihood estimator is", ...
                     " not regular, so standard errors and confidence", ...
                     " intervals cannot be computed reliably."));
    reachedBnd = true;
  else
    reachedBnd = false;
  endif

  ## If second output argument is requested
  if (nargout > 1)
    if (! reachedBnd)
      probs = [alpha/2; 1-alpha/2];
      [~, acov] = gplike (paramhat, x);
      se = sqrt (diag (acov))';
      ## Compute the CI for shape using a normal distribution for khat.
      kci = norminv (probs, paramhat(1), se(1));
      ## Compute the CI for scale using a normal approximation for
      ## log(sigmahat), and transform back to the original scale.
      lnsigci = norminv (probs, log (paramhat(2)), se(2) ./ paramhat(2));
      paramci = [kci, exp(lnsigci)];
    else
      paramci = [NaN, NaN; NaN, NaN];
    endif
  endif

endfunction

## Negative log-likelihood for the GP
function nll = negloglike (paramhat, data)

  shape = paramhat(1);
  log_scale = paramhat(2);
  scale = exp (log_scale);
  sample_size = numel (data);
  z = data ./ scale;
  if (abs (shape) > eps)
    if (shape > 0 || max (z) < -1 / shape)
      nll = sample_size * log_scale + (1 + 1/shape) * sum (log1p (shape .* z));
    else
      nll = Inf;
    endif
  else
    nll = sample_size * log_scale + sum (z);
  endif

endfunction

%!demo
%! ## Sample 2 populations from different generalized Pareto distributions
%! ## Assume location parameter θ is known
%! rng (42);
%! theta = 0;
%! r1 = gprnd (1, 2, theta, 20000, 1);
%! r2 = gprnd (3, 1, theta, 20000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0.1:0.2:100], 5);
%! h = findobj (gca, 'Type', 'patch');
%! set (h(1), 'facecolor', 'r');
%! set (h(2), 'facecolor', 'c');
%! ylim ([0, 1]);
%! xlim ([0, 5]);
%! hold on
%!
%! ## Estimate their α and β parameters
%! k_sigmaA = gpfit (r(:,1));
%! k_sigmaB = gpfit (r(:,2));
%!
%! ## Plot their estimated PDFs
%! x = [0.01, 0.1:0.2:18];
%! y = gppdf (x, k_sigmaA(1), k_sigmaA(2), theta);
%! plot (x, y, '-pc');
%! y = gppdf (x, k_sigmaB(1), k_sigmaB(2), theta);
%! plot (x, y, '-sr');
%! hold off
%! legend ({'Normalized HIST of sample 1 with k=1 and σ=2', ...
%!          'Normalized HIST of sample 2 with k=2 and σ=2', ...
%!          sprintf("PDF for sample 1 with estimated k=%0.2f and σ=%0.2f", ...
%!                  k_sigmaA(1), k_sigmaA(2)), ...
%!          sprintf("PDF for sample 3 with estimated k=%0.2f and σ=%0.2f", ...
%!                  k_sigmaB(1), k_sigmaB(2))})
%! title ('Two population samples from different generalized Pareto distributions')
%! text (2, 0.7, 'Known location parameter θ = 0')
%! hold off

## Test output
## Values below are R2024a's, measured 2026-08-17.  The estimates agree to the
## convergence tolerance of fminsearch, which is what the 1e-4 covers.
%!shared x
%! x = [2.2196, 11.9301, 4.3673, 1.0949, 6.5626, ...
%!      1.2109, 1.8576, 1.0039, 12.7917, 2.2590];
%!test
%! [hat, ci] = gpfit (x);
%! assert_equal (hat, [-0.163107819293798, 5.305483917184919], 1e-4);
%! assert_equal (ci, [-1.174106637867584,  1.627748133572634; ...
%!                     0.847890999279987, 17.292699659699391], 1e-4);
%!test
%! [hat, ci] = gpfit (x, 0.10);
%! assert_equal (ci, [-1.011564773958192,  1.968276868273005; ...
%!                     0.685349135370595, 14.300914698146826], 1e-4);
%!test
%! ## a known location is fitted by shifting the data, and only shifts the fit
%! [hat, ci] = gpfit (x - 1);
%! assert_equal (hat, [0.893710299404345, 1.322962458731574], 1e-6);
%! assert_equal (ci, [-0.774991092191746, 0.243695078371714; ...
%!                     2.562411691000436, 7.182047659343478], 1e-5);
%!assert_equal (size (gpfit (x)), [1, 2])
%!test
%! [~, ci] = gpfit (x);
%! assert_equal (size (ci), [2, 2]);
%!test
%! ## the default confidence level is 95%
%! [h1, c1] = gpfit (x);
%! [h2, c2] = gpfit (x, 0.05);
%! assert_equal (h1, h2);
%! assert_equal (c1, c2);
%!test
%! ## FREQ counts repeated observations
%! assert_equal (gpfit (x, [], [], [2, ones(1,9)]), gpfit ([x(1), x]), 1e-10);
%!test
%! ## non-finite data propagates into the estimates instead of being dropped
%! assert_equal (gpfit ([x, NaN]), [NaN, NaN]);
%!test
%! assert_equal (gpfit ([x, Inf]), [NaN, NaN]);
%!test
%! assert_equal (gpfit ([x, NaN, Inf]), [NaN, NaN]);
%!test
%! ## below a shape of -1 the likelihood is unbounded, but the estimate still
%! ## keeps every observation strictly inside the fitted support
%! xb = [1.2 2.3 0.5 3.1 2.2 1.8 0.9 2.7 1.1 3.3];
%! warning ('off', 'all');
%! p = gpfit (xb);
%! warning ('on', 'all');
%! assert_equal (p(1) < -1, true);
%! assert_equal (max (xb) < -p(2) / p(1), true);
%! assert_equal (isfinite (gplike (p, xb)), true);
%!test
%! ## the confidence intervals are withheld there
%! warning ('off', 'all');
%! [~, ci] = gpfit ([1.2 2.3 0.5 3.1 2.2 1.8 0.9 2.7 1.1 3.3]);
%! warning ('on', 'all');
%! assert_equal (ci, [NaN, NaN; NaN, NaN]);

%!warning<gpfit: the fitted upper bound of the support has closed onto the largest observation, a boundary of the parameter space, so the estimates are unreliable and no confidence intervals are computed.> ...
%! gpfit ([1.2 2.3 0.5 3.1 2.2 1.8 0.9 2.7 1.1 3.3]);

## Test input validation
%!error<gpfit: function called with too few input arguments.> gpfit ()
%!error<gpfit: X must be a vector of real values.> gpfit ([0.2, 0.5+i]);
%!error<gpfit: X must be a vector of real values.> gpfit (ones (2,2) * 0.5);
%!error<gpfit: X must not contain negative values.> gpfit ([-1, 2, 3]);
%!error<gpfit: wrong value for ALPHA.> gpfit ([0.01:0.1:0.99], 1.2);
%!error<gpfit: wrong value for ALPHA.> gpfit ([0.01:0.1:0.99], i);
%!error<gpfit: wrong value for ALPHA.> gpfit ([0.01:0.1:0.99], -1);
%!error<gpfit: wrong value for ALPHA.> gpfit ([0.01:0.1:0.99], [0.05, 0.01]);
%!error<gpfit: X and FREQ vectors mismatch.> ...
%! gpfit ([1 2 3], [], [], [1 5])
%!error<gpfit: FREQ must not contain negative values.> ...
%! gpfit ([1 2 3], [], [], [1 5 -1])
%!error<gpfit: 'options' 3rd argument must be a structure> ...
%! gpfit ([1:10], 0.05, 5)
