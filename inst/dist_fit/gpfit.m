## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefnx {statistics} {[@dots{}] =} gpfit (@var{x}, @var{alpha}, @var{options})
##
## Estimate parameters and confidence intervals for the generalized Pareto
## distribution.
##
## @code{@var{paramhat} = gpfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the generalized Pareto distribution given the
## data in @var{x}.  @qcode{@var{paramhat}(1)} is the shape parameter, @var{k},
## and @qcode{@var{paramhat}(2)} is the scale parameter, @var{sigma}.  Other
## functions for the generalized Pareto, such as @code{gpcdf}, allow a location
## parameter, @var{mu}.  However, @code{gpfit} does not estimate a location
## parameter, and it must be assumed known, and subtracted from @var{x} before
## calling @code{gpfit}.
##
## @code{[@var{paramhat}, @var{paramci}] = gpfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = gpfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = gpfit (@var{x}, @var{alpha}, @var{options})} specifies
## control parameters for the iterative algorithm used to compute ML estimates
## with the @code{fminsearch} function.  @var{options} is a structure with the
## following fields and their default values:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.MaxFunEvals = 400}
## @item @qcode{@var{options}.MaxIter = 200}
## @item @qcode{@var{options}.TolBnd = 1e-6}
## @item @qcode{@var{options}.TolFun = 1e-6}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
##
## When @qcode{@var{k} = 0} and @qcode{@var{mu} = 0}, the Generalized Pareto CDF
## is equivalent to the exponential distribution.  When @qcode{@var{k} > 0} and
## @code{@var{mu} = @var{k} / @var{k}} the Generalized Pareto is equivalent to
## the Pareto distribution.  The mean of the Generalized Pareto is not finite
## when @qcode{@var{k} >= 1} and the variance is not finite when
## @qcode{@var{k} >= 1/2}.  When @qcode{@var{k} >= 0}, the Generalized Pareto
## has positive density for @qcode{@var{x} > @var{mu}}, or, when
## @qcode{@var{mu} < 0}, for
## @qcode{0 <= (@var{x} - @var{mu}) / @var{sigma} <= -1 / @var{k}}.
##
## Further information about the generalized Pareto distribution can be found at
## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
##
## @seealso{gpcdf, gpinv, gppdf, gprnd, gplike, gpstat}
## @end deftypefn

function [paramhat, paramci] = gpfit (x, alpha, options)

  ## Check input arguments, X must be a vector of positive values
  if (! isvector (x))
    error ("gpfit: X must be a vector.");
  endif
  if (any (x <= 0))
    error ("gpfit: X must contain only positive values.");
  endif

  ## Add default value for alpha if not supplied
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  endif

  ## Check for valid value of alpha
  if (! isscalar (alpha) || ! isnumeric (alpha) || alpha <= 0 || alpha >= 1)
     error ("gpfit: wrong value for alpha.");
  endif

  ## Add default values to OPTIONS structure
  if (nargin < 3 || isempty (options))
    options.Display = "off";
    options.MaxFunEvals = 400;
    options.MaxIter = 200;
    options.TolBnd = 1e-6;
    options.TolFun = 1e-6;
    options.TolX = 1e-6;
  elseif (! isstruct (options))
    error ("gpfit: OPTIONS must be a structure for 'fminsearch' function.");
  endif
  if (! isfield (options, "Display"))
    options.Display = "off";
  endif
  if (! isfield (options, "MaxFunEvals"))
    options.MaxFunEvals = 400;
  endif
  if (! isfield (options, "MaxIter"))
    options.MaxIter = 200;
  endif
  if (! isfield (options, "TolBnd"))
    options.TolBnd = 1e-6;
  endif
  if (! isfield (options, "TolFun"))
    options.TolFun = 1e-6;
  endif
  if (! isfield (options, "TolX"))
    options.TolX = 1e-6;
  endif

  ## Get class of X
  is_class = class (x);
  if (strcmp (is_class, "single"))
    x = double (x);
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
  x_max = max (x);
  x_size = length (x);
  x_range = range (x);

  ## Check for appropriate sample size or all observations being equal
  if (x_size == 0)
    paramhat = NaN (1, 2, is_class);
    paramci = NaN (2, 2, is_class);
    warning ("gpfit: X contains no data.");
    return
  elseif (x_range < realmin (is_class))
    if (x_max <= sqrt (realmax (is_class)))
      paramhat = cast ([NaN 0], is_class);
    endif
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
  [paramhat, ~, err, output] = fminsearch (@negloglike, paramhat, options, x);
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
  if ((paramhat(1) < 0) && (x_max > -paramhat(2)/paramhat(1) - options.TolBnd))
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

  ## Preserve class tu output arguments
  paramhat = cast (paramhat, is_class);
  if (nargout > 1)
    paramci = cast (paramci, is_class);
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
%! ## Sample 2 populations from different generalized Pareto distibutions
%! ## Assume location parameter is known
%! mu = 0;
%! rand ("seed", 5);    # for reproducibility
%! r1 = gprnd (1, 2, mu, 20000, 1);
%! rand ("seed", 2);    # for reproducibility
%! r2 = gprnd (3, 1, mu, 20000, 1);
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
%! k_sigmaA = gpfit (r(:,1));
%! k_sigmaB = gpfit (r(:,2));
%!
%! ## Plot their estimated PDFs
%! x = [0.01, 0.1:0.2:18];
%! y = gppdf (x, k_sigmaA(1), k_sigmaA(2), mu);
%! plot (x, y, "-pc");
%! y = gppdf (x, k_sigmaB(1), k_sigmaB(2), mu);
%! plot (x, y, "-sr");
%! hold off
%! legend ({"Normalized HIST of sample 1 with k=1 and σ=2", ...
%!          "Normalized HIST of sample 2 with k=2 and σ=2", ...
%!          sprintf("PDF for sample 1 with estimated k=%0.2f and σ=%0.2f", ...
%!                  k_sigmaA(1), k_sigmaA(2)), ...
%!          sprintf("PDF for sample 3 with estimated k=%0.2f and σ=%0.2f", ...
%!                  k_sigmaB(1), k_sigmaB(2))})
%! title ("Three population samples from different generalized Pareto distibutions")
%! text (2, 0.7, "Known location parameter μ = 0")
%! hold off

## Test output
%!test
%! shape = 5; scale = 2;
%! x = gprnd (shape, scale, 0, 1, 100000);
%! [hat, ci] = gpfit (x);
%! assert (hat, [shape, scale], 1e-1);
%! assert (ci, [shape, scale; shape, scale], 2e-1);
%!test
%! shape = 1; scale = 1;
%! x = gprnd (shape, scale, 0, 1, 100000);
%! [hat, ci] = gpfit (x);
%! assert (hat, [shape, scale], 1e-1);
%! assert (ci, [shape, scale; shape, scale], 1e-1);
%!test
%! shape = 3; scale = 2;
%! x = gprnd (shape, scale, 0, 1, 100000);
%! [hat, ci] = gpfit (x);
%! assert (hat, [shape, scale], 1e-1);
%! assert (ci, [shape, scale; shape, scale], 1e-1);

## Test input validation
%!error<gpfit: X must be a vector.> gpfit (ones (2))
%!error<gpfit: X must contain only positive values.> gpfit ([-1, 2])
%!error<gpfit: X must contain only positive values.> gpfit ([0, 1, 2])
%!error<gpfit: wrong value for alpha.> gpfit ([1, 2], 0)
%!error<gpfit: wrong value for alpha.> gpfit ([1, 2], 1.2)
%!error<gpfit: OPTIONS must be a structure for 'fminsearch' function.> ...
%! gpfit ([1:10], 0.05, 5)
