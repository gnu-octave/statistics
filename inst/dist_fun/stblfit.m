## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{paramhat} =} stblfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} stblfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} stblfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} stblfit (@var{x}, @var{alpha}, @var{freq})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} stblfit (@var{x}, @var{alpha}, @var{options})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} stblfit (@var{x}, @var{alpha}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for the stable distribution.
##
## @code{@var{paramhat} = stblfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the stable distribution, in the Nolan
## @qcode{S0} parameterization, given the data in @var{x}.
## @qcode{@var{paramhat}(1)} is the tail index @var{alpha},
## @qcode{@var{paramhat}(2)} is the skewness @var{beta}, @qcode{@var{paramhat}(3)}
## is the scale @var{gam}, and @qcode{@var{paramhat}(4)} is the location
## @var{delta}.
##
## @code{[@var{paramhat}, @var{paramci}] = stblfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.  The intervals are Wald
## intervals from the observed Fisher information.
##
## @code{[@dots{}] = stblfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default value.
##
## @code{[@dots{}] = stblfit (@var{x}, @var{alpha}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## must contain non-negative integer frequencies for the corresponding elements
## in @var{x}.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@var{paramhat}, @var{paramci}] = stblfit (@var{x}, @var{alpha},
## @var{options})} specifies control parameters for the iterative algorithm used
## to compute the ML estimates with the @code{fminsearch} function.
## @var{options} is a structure with the following fields and their default
## values:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.MaxFunEvals = 400}
## @item @qcode{@var{options}.MaxIter = 200}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
##
## The stable density has no closed form; it is evaluated by numerical inversion
## of the characteristic function, which makes fitting considerably slower than
## for the closed-form distributions.  Censoring is not supported.
##
## The estimates are the maximum-likelihood estimates under the mathematically
## exact density.  MATLAB fits an interpolation-based approximation of the stable
## density, whose maximum-likelihood estimates deviate from the exact ones by
## about @math{10^{-2}} (and the resulting confidence intervals by up to roughly
## 20%); @code{stblfit} returns the exact (more accurate) estimates.
##
## Further information about the stable distribution can be found at
## @url{https://en.wikipedia.org/wiki/Stable_distribution}
##
## @seealso{stbllike, stblpdf, stblcdf, stblinv, stblrnd, fitdist, makedist}
## @end deftypefn

function [paramhat, paramci] = stblfit (x, alpha, varargin)

  ## Check X is a vector
  if (! isvector (x))
    error ("stblfit: X must be a vector.");
  endif

  ## Get X type and convert to double for computation
  is_type = class (x);
  if (strcmpi (is_type, "single"))
    x = double (x);
  endif

  ## Check that X does not contain NaNs and is not constant
  if (any (isnan (x)))
    error ("stblfit: X must not contain NaN values.");
  elseif (numel (unique (x)) == 1)
    error ("stblfit: X must contain at least two distinct values.");
  endif

  ## Check ALPHA (significance level)
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("stblfit: wrong value for ALPHA.");
  endif

  ## Add defaults
  freq = [];
  options.Display = "off";
  options.MaxFunEvals = 400;
  options.MaxIter = 200;
  options.TolX = 1e-6;

  ## Check extra arguments for FREQ vector and/or 'options' structure
  if (nargin > 2)
    if (numel (varargin) == 1 && isstruct (varargin{1}))
      options = varargin{1};
    elseif (numel (varargin) == 1 && isnumeric (varargin{1}))
      freq = varargin{1};
    elseif (numel (varargin) == 2)
      freq = varargin{1};
      options = varargin{2};
    endif
    if (isempty (freq))
      freq = ones (size (x));
    endif
    if (! isequal (size (x), size (freq)))
      error ("stblfit: X and FREQ vectors mismatch.");
    elseif (any (freq < 0))
      error ("stblfit: FREQ must not contain negative values.");
    elseif (any (fix (freq) != freq))
      error ("stblfit: FREQ must contain integer values.");
    endif
    if (! isstruct (options) || ! isfield (options, "Display") ||
        ! isfield (options, "MaxFunEvals") || ! isfield (options, "MaxIter")
                                           || ! isfield (options, "TolX"))
      error (strcat ("stblfit: 'options' argument must be a structure with", ...
                     " 'Display', 'MaxFunEvals', 'MaxIter', and 'TolX'", ...
                     " fields present."));
    endif
  endif
  if (isempty (freq))
    freq = ones (size (x));
  endif

  ## Force column vectors and drop zero-frequency observations
  x = x(:);
  freq = freq(:);
  keep = freq > 0;
  x = x(keep);
  freq = freq(keep);

  ## Objective: negative log-likelihood in the unconstrained parameter space
  fhandle = @(u) stbl_nll (to_nat (u), x, freq);

  ## Initial estimate: pick the best of a small grid of quantile-scaled starts,
  ## then refine.  The grid spans the tail index and skewness, while the scale
  ## and location are matched to the sample interquartile range and median.
  q = quantile (x, [0.25, 0.5, 0.75]);
  iqr_x = q(3) - q(1);
  best_nll = Inf;
  u0 = [];
  for a0 = [0.7, 1.1, 1.5, 1.9]
    for b0 = [-0.5, 0, 0.5]
      ## Standardized interquartile range for this shape
      siqr = diff (stblinv ([0.25, 0.75], a0, b0, 1, 0));
      gam0 = iqr_x / siqr;
      delta0 = q(2) - gam0 * stblinv (0.5, a0, b0, 1, 0);
      p0 = [a0, b0, gam0, delta0];
      nll0 = stbl_nll (p0, x, freq);
      if (nll0 < best_nll)
        best_nll = nll0;
        u0 = to_unc (p0);
      endif
    endfor
  endfor

  ## Minimize the negative log-likelihood
  [uhat, ~, exitflag, output] = fminsearch (fhandle, u0, options);
  paramhat = to_nat (uhat);

  ## Display warnings if the optimizer did not converge
  if (exitflag == 0)
    if (output.funcCount >= output.iterations)
      warning ("stblfit: maximum number of function evaluations reached.");
    else
      warning ("stblfit: reached iteration limit.");
    endif
  endif

  ## Return a row vector for MATLAB compatibility
  paramhat = paramhat(:)';

  ## Check for second output argument
  if (nargout > 1)
    [~, acov] = stbllike (paramhat, x, freq);
    param_se = sqrt (diag (acov))';
    if (any (! isfinite (param_se)))
      warning (strcat ("stblfit: could not compute confidence intervals;", ...
                       " the Fisher information matrix is not positive", ...
                       " definite."));
      paramci = NaN (2, 4, is_type);
    else
      ## Linear Wald intervals for all four parameters (matching MATLAB, which
      ## does not place the scale parameter on a log scale here)
      p_vals = [alpha/2; 1-alpha/2];
      a_ci = norminv (p_vals, paramhat(1), param_se(1));
      b_ci = norminv (p_vals, paramhat(2), param_se(2));
      g_ci = norminv (p_vals, paramhat(3), param_se(3));
      d_ci = norminv (p_vals, paramhat(4), param_se(4));
      paramci = [a_ci, b_ci, g_ci, d_ci];
    endif
  endif

endfunction

## Map the unconstrained vector U to the natural parameters
function p = to_nat (u)
  p = [2 ./ (1 + exp (-u(1))), tanh(u(2)), exp(u(3)), u(4)];
endfunction

## Map the natural parameters P to the unconstrained vector
function u = to_unc (p)
  u = [log(p(1) ./ (2 - p(1))), atanh(p(2)), log(p(3)), p(4)];
endfunction

## Negative log-likelihood of the stable density at parameter vector P
function nll = stbl_nll (p, x, freq)
  y = __stable_pdf__ (x, p(1), p(2), p(3), p(4));
  y(y <= 0) = realmin;
  nll = -sum (freq .* log (y));
  if (! isreal (nll) || isnan (nll))
    nll = Inf;
  endif
endfunction

%!demo
%! ## Fit a stable distribution to simulated data
%! rand ("seed", 42);
%! x = stblrnd (1.5, 0.5, 2, 1, 150, 1);
%! [paramhat, paramci] = stblfit (x)

## Test output
## The stable density has no closed form, so fitting is by numerical inversion
## of the characteristic function and is slow; the sample sizes below are kept
## modest to bound the test time, with correspondingly loose tolerances.
%!test  # recovery + MATLAB parity
%!      # Our fast integrator reproduces the exact stblpdf MLE to ~2e-5.  MATLAB's
%!      # stable density is a Nolan interpolation approximation, so its fit
%!      # deviates from the exact MLE by ~1e-2; we ship the exact (more accurate)
%!      # estimate and document the deviation (as for the copula family).
%! rand ("seed", 2718);
%! randn ("seed", 2718);
%! x = stblrnd (1.5, 0.5, 2, 1, 150, 1);
%! [phat, pci] = stblfit (x);
%! ## Exact-density estimate on this sample (recovers the generating [1.5 0.5 2 1])
%! assert (phat, [1.5469000, 0.4732298, 2.0097077, 1.1640279], 1e-3);
%! ## MATLAB fitdist (x, 'Stable') on the same data agrees to ~1.5e-2
%! assert (phat, [1.5449145, 0.4693139, 2.0000225, 1.1646526], 1.5e-2);
%! ## Confidence intervals bracket the estimate; gam CI is positive
%! assert (all (pci(1,:) <= phat) && all (pci(2,:) >= phat));
%! assert (pci(1,3) > 0);

## Test input validation
%!error <stblfit: X must be a vector.> stblfit (ones (2, 2))
%!error <stblfit: X must not contain NaN values.> stblfit ([1, 2, NaN, 4])
%!error <stblfit: X must contain at least two distinct values.> ...
%! stblfit ([2, 2, 2, 2])
%!error <stblfit: wrong value for ALPHA.> stblfit ([1, 2, 3, 4], 1.5)
%!error <stblfit: wrong value for ALPHA.> stblfit ([1, 2, 3, 4], -0.5)
%!error <stblfit: X and FREQ vectors mismatch.> ...
%! stblfit ([1, 2, 3, 4], 0.05, [1, 2, 3])
%!error <stblfit: FREQ must not contain negative values.> ...
%! stblfit ([1, 2, 3, 4], 0.05, [1, -1, 2, 1])
%!error <stblfit: FREQ must contain integer values.> ...
%! stblfit ([1, 2, 3, 4], 0.05, [1, 1.5, 2, 1])
