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
## @deftypefn  {statistics} {@var{paramhat} =} gpfit (@var{x}, @var{theta})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gpfit (@var{x}, @var{theta})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gpfit (@var{x}, @var{theta}, @var{alpha})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gpfit (@var{x}, @var{theta}, @var{alpha}, @var{options})
##
## Estimate parameters and confidence intervals for the generalized Pareto
## distribution.
##
## @code{@var{paramhat} = gpfit (@var{x}, @var{theta})} returns the maximum
## likelihood estimates of the parameters of the generalized Pareto distribution
## given the data in @var{x} and the location parameter @var{theta}.
## @qcode{@var{paramhat}(1)} is the shape parameter, @var{k},
## @qcode{@var{paramhat}(2)} is the scale parameter, @var{sigma}, and
## @qcode{@var{paramhat}(3)} is the location parameter, @var{theta}.  Although
## @var{theta} is returned in the estimated @var{paramhat}, @code{gpfit} does
## not estimate the location parameter @var{theta}, and it must be assumed to be
## known, given as a fixed parameter in input argument @var{theta}.
##
## @code{[@var{paramhat}, @var{paramci}] = gpfit (@var{x}, @var{theta})} returns
## the 95% confidence intervals for the estimated parameter @var{k} and
## @var{sigma}. The third column of @var{paramci} includes the location
## parameter @var{theta} without any confidence bounds.
##
## @code{[@dots{}] = gpfit (@var{x}, @var{theta}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = gpfit (@var{x}, @var{theta}, @var{alpha}, @var{options})}
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

function [paramhat, paramci] = gpfit (x, theta, alpha, freq, options)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("gpfit: function called with too few input arguments.");
  endif

  ## Check X for being a vector
  if (isempty (x))
    phat = nan (1, 3, class (x));
    pci = nan (3, 3, class (x));
    return
  elseif (! isvector (x) || ! isreal (x))
    error ("gpfit: X must be a vector of real values.");
  endif

  ## Check for THETA being a scalar real value
  if (! isscalar (theta) || ! isreal (theta))
    error ("gpfit: THETA must be a real scalar value.");
  endif

  ## Check X >= THETA
  if (any (x < theta))
    error ("gpfit: X cannot contain values less than THETA.");
  endif

  ## Parse ALPHA argument or add default
  if (nargin < 3 || isempty (alpha))
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
  if (nargin < 5)
    options.Display = "off";
    options.MaxFunEvals = 400;
    options.MaxIter = 200;
    options.TolX = 1e-6;
  else
    if (! isstruct (options) || ! isfield (options, "Display") ||
        ! isfield (options, "MaxFunEvals") || ! isfield (options, "MaxIter")
                                           || ! isfield (options, "TolX"))
      error (strcat (["gpfit: 'options' 5th argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  endif

  ## Remove NaN from data and make a warning
  if (any (isnan (x)))
    x(isnan(x)) = [];
    warning ("gpfit: X contains NaN values, which are ignored.");
  endif

  ## Remove Inf from data and make a warning
  if (any (isinf (x)))
    x(isnan(x)) = [];
    warning ("gpfit: X contains Inf values, which are ignored.");
  endif

  ## Get sample size, max and range of X
  x = x - theta;
  x_max = max (x);
  x_size = length (x);
  x_range = range (x);

  ## Check for appropriate sample size or all observations being equal
  if (x_size == 0)
    paramhat = [NaN(1, 2), theta];
    paramci = [NaN(2, 2), [theta; theta]];
    warning ("gpfit: X contains no data.");
    return
  elseif (x_range < realmin (class (x)))
    paramhat = cast ([NaN, 0, theta], class (x));
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
  paramhat = [paramhat, theta];

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
    warning ("gpfit: converged to boundary 1.");
    reachedBnd = true;
  elseif (paramhat(1) <= -1 / 2)
    warning ("gpfit: converged to boundary 2.");
    reachedBnd = true;
  else
    reachedBnd = false;
  endif

  ## If second output argument is requested
  if (nargout > 1)
    if (! reachedBnd)
      probs = [alpha/2; 1-alpha/2];
      [~, acov] = gplike (paramhat, x + theta);
      se = sqrt (diag (acov))';
      ## Compute the CI for shape using a normal distribution for khat.
      kci = norminv (probs, paramhat(1), se(1));
      ## Compute the CI for scale using a normal approximation for
      ## log(sigmahat), and transform back to the original scale.
      lnsigci = norminv (probs, log (paramhat(2)), se(2) ./ paramhat(2));
      paramci = [kci, exp(lnsigci), [theta; theta]];
    else
      paramci = [NaN, NaN, theta; NaN, NaN, theta];
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
%! theta = 0;
%! rand ("seed", 5);    # for reproducibility
%! r1 = gprnd (1, 2, theta, 20000, 1);
%! rand ("seed", 2);    # for reproducibility
%! r2 = gprnd (3, 1, theta, 20000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0.1:0.2:100], 5);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "r");
%! set (h(2), "facecolor", "c");
%! ylim ([0, 1]);
%! xlim ([0, 5]);
%! hold on
%!
%! ## Estimate their α and β parameters
%! k_sigmaA = gpfit (r(:,1), theta);
%! k_sigmaB = gpfit (r(:,2), theta);
%!
%! ## Plot their estimated PDFs
%! x = [0.01, 0.1:0.2:18];
%! y = gppdf (x, k_sigmaA(1), k_sigmaA(2), theta);
%! plot (x, y, "-pc");
%! y = gppdf (x, k_sigmaB(1), k_sigmaB(2), theta);
%! plot (x, y, "-sr");
%! hold off
%! legend ({"Normalized HIST of sample 1 with k=1 and σ=2", ...
%!          "Normalized HIST of sample 2 with k=2 and σ=2", ...
%!          sprintf("PDF for sample 1 with estimated k=%0.2f and σ=%0.2f", ...
%!                  k_sigmaA(1), k_sigmaA(2)), ...
%!          sprintf("PDF for sample 3 with estimated k=%0.2f and σ=%0.2f", ...
%!                  k_sigmaB(1), k_sigmaB(2))})
%! title ("Two population samples from different generalized Pareto distributions")
%! text (2, 0.7, "Known location parameter θ = 0")
%! hold off

## Test output
%!test
%! k = 0.8937; sigma = 1.3230; theta = 1;
%! x = [2.2196, 11.9301, 4.3673, 1.0949, 6.5626, ...
%!      1.2109, 1.8576, 1.0039, 12.7917, 2.2590];
%! [hat, ci] = gpfit (x, theta);
%! assert (hat, [k, sigma, theta], 1e-4);
%! assert (ci, [-0.7750, 0.2437, 1; 2.5624, 7.1820, 1], 1e-4);

## Test input validation
%!error<gpfit: function called with too few input arguments.> gpfit ()
%!error<gpfit: function called with too few input arguments.> gpfit (1)
%!error<gpfit: X must be a vector of real values.> gpfit ([0.2, 0.5+i], 0);
%!error<gpfit: X must be a vector of real values.> gpfit (ones (2,2) * 0.5, 0);
%!error<gpfit: THETA must be a real scalar value.> ...
%! gpfit ([0.5, 1.2], [0, 1]);
%!error<gpfit: THETA must be a real scalar value.> ...
%! gpfit ([0.5, 1.2], 5+i);
%!error<gpfit: X cannot contain values less than THETA.> ...
%! gpfit ([1:5], 2);
%!error<gpfit: wrong value for ALPHA.> gpfit ([0.01:0.1:0.99], 0, 1.2);
%!error<gpfit: wrong value for ALPHA.> gpfit ([0.01:0.1:0.99], 0, i);
%!error<gpfit: wrong value for ALPHA.> gpfit ([0.01:0.1:0.99], 0, -1);
%!error<gpfit: wrong value for ALPHA.> gpfit ([0.01:0.1:0.99], 0, [0.05, 0.01]);
%!error<gpfit: X and FREQ vectors mismatch.>
%! gpfit ([1 2 3], 0, [], [1 5])
%!error<gpfit: FREQ must not contain negative values.>
%! gpfit ([1 2 3], 0, [], [1 5 -1])
%!error<gpfit: 'options' 5th argument must be a structure> ...
%! gpfit ([1:10], 1, 0.05, [], 5)
