## Copyright (C) 2012-2021 Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{paramhat} =} gevfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gevfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gevfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gevfit (@var{x}, @var{alpha}, @var{freq})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gevfit (@var{x}, @var{alpha}, @var{options})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gevfit (@var{x}, @var{alpha}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for the generalized extreme
## value (GEV) distribution.
##
## @code{@var{paramhat} = gevfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the GEV distribution given the data in
## @var{x}.  @qcode{@var{paramhat}(1)} is the shape parameter, @var{k}, and
## @qcode{@var{paramhat}(2)} is the scale parameter, @var{sigma}, and
## @qcode{@var{paramhat}(3)} is the location parameter, @var{mu}.
##
## @code{[@var{paramhat}, @var{paramci}] = gevfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = gevfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = gevfit (@var{params}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## must contain non-negative integer frequencies for the corresponding elements
## in @var{x}.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@var{paramhat}, @var{paramci}] = gevfit (@var{x}, @var{alpha},
## @var{options})} specifies control parameters for the iterative algorithm used
## to compute ML estimates with the @code{fminsearch} function.  @var{options}
## is a structure with the following fields and their default values:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.MaxFunEvals = 400}
## @item @qcode{@var{options}.MaxIter = 200}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
##
## When @qcode{@var{k} < 0}, the GEV is the type III extreme value distribution.
## When @qcode{@var{k} > 0}, the GEV distribution is the type II, or Frechet,
## extreme value distribution.  If @var{W} has a Weibull distribution as
## computed by the @code{wblcdf} function, then @qcode{-@var{W}} has a type III
## extreme value distribution and @qcode{1/@var{W}} has a type II extreme value
## distribution.  In the limit as @var{k} approaches @qcode{0}, the GEV is the
## mirror image of the type I extreme value distribution as computed by the
## @code{evcdf} function.
##
## The mean of the GEV distribution is not finite when @qcode{@var{k} >= 1}, and
## the variance is not finite when @qcode{@var{k} >= 1/2}.  The GEV distribution
## has positive density only for values of @var{x} such that
## @qcode{@var{k} * (@var{x} - @var{mu}) / @var{sigma} > -1}.
##
## Further information about the generalized extreme value distribution can be
## found at
## @url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
##
## @seealso{gevcdf, gevinv, gevpdf, gevrnd, gevlike, gevstat}
## @end deftypefn

function [paramhat, paramci] = gevfit (x, alpha, varargin)

  ## Check X is vector
  if (! isvector (x))
    error ("gevfit: X must be a vector.");
  endif

  ## Get X type and convert to double for computation
  is_type = class (x);
  if (strcmpi (is_type, "single"))
    x = double (x);
  endif

  ## Check that X is not constant and does not contain NaNs
  sample_size = length (x);
  if (sample_size == 0 || any (isnan (x)))
    paramhat = NaN (1,3, is_type);
    paramci = NaN (2,3, is_type);
    warning ("gevfit: X contains NaNs.");
    return
  elseif (numel (unique (x)) == 1)
    paramhat = cast ([0, 0, unique(x)], is_type);
    if (length (x) == 1)
      paramci = cast ([-Inf, 0, -Inf; Inf, Inf, Inf], is_type);
    else
      paramci = [paramhat; paramhat];
    endif
    warning ("gevfit: X is a constant vector.");
    return
  endif

  ## Check ALPHA
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("gevfit: wrong value for ALPHA.");
    endif
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
    ## Check for valid freq vector
    if (! isequal (size (x), size (freq)))
      error ("gevfit: X and FREQ vectors mismatch.");
    elseif (any (freq < 0))
      error ("gevfit: FREQ must not contain negative values.");
    elseif (any (fix (freq) != freq))
      error ("gevfit: FREQ must contain integer values.");
    endif
    ## Check for valid options structure
    if (! isstruct (options) || ! isfield (options, "Display") ||
        ! isfield (options, "MaxFunEvals") || ! isfield (options, "MaxIter")
                                           || ! isfield (options, "TolX"))
      error (strcat (["gevfit: 'options' argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  endif

  ## Expand frequency
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  ## Force to column vector
  x = x(:);

  ## Compute initial parameters
  F = (0.5:1:(sample_size - 0.5))' ./ sample_size;
  k_0 = fminsearch (@(k) 1 - corr (x, gevinv (F, k, 1, 0)), 0);
  paramguess = [k_0, polyfit(gevinv(F,k_0,1,0),x',1)];

  ## Check if x support initial parameters or fall back to unbounded evfit
  if (k_0 < 0 && (max (x) > - paramguess(2) / k_0 + paramguess(3)) || ...
      k_0 > 0 && (min (x) < - paramguess(2) / k_0 + paramguess(3)))
    paramguess = [evfit(x), 0];
    paramguess = flip (paramguess);
  endif

  ## Minimize the negative log-likelihood according to initial parameters
  paramguess(2) = log (paramguess(2));
  fhandle = @(paramguess) nll (paramguess, x);
  [paramhat, ~, exitflag, output] = fminsearch (fhandle, paramguess, options);
  paramhat(2) = exp (paramhat(2));

  ## Display errors and warnings if any
  if (exitflag == 0)
    if (output.funcCount >= output.iterations)
      warning ("gevfit: maximum number of evaluations reached");
    else
      warning ("gevfit: reached iteration limit");
    endif
  elseif (exitflag == -1)
    error ("gevfit: No solution");
  endif

  ## Return a row vector for Matlab compatibility
  paramhat = paramhat(:)';

  ## Check for second output argument
  if (nargout > 1)
  	[~, acov] = gevlike (paramhat, x);
  	param_se = sqrt (diag (acov))';
    if (any (iscomplex (param_se)))
      warning (["gevfit: Fisher information matrix not positive definite;", ...
                " parameter optimization likely did not converge"]);
      paramci = NaN (2, 3, is_type);
    else
      p_vals = [alpha/2; 1-alpha/2];
      k_ci = norminv (p_vals, paramhat(1), param_se(1));
      s_ci = exp (norminv (p_vals, log (paramhat(2)), param_se(2) ./ paramhat(2)));
      m_ci = norminv (p_vals, paramhat(3), param_se(3));
      paramci = [k_ci, s_ci, m_ci];
    endif
  endif
endfunction

## Negative log-likelihood for the GEV (log(sigma) parameterization)
function out = nll (parms, x)
  k_0 = parms(1);
  log_sigma = parms(2);
  sigma = exp (log_sigma);
  mu = parms(3);
  n = numel (x);
  z = (x - mu) ./ sigma;
  if abs(k_0) > eps
      u = 1 + k_0.*z;
      if min(u) > 0
          lnu = log1p (k_0 .* z);
          out = n * log_sigma + sum (exp ((-1 / k_0) * lnu)) + ...
                                (1 + 1 / k_0) * sum (lnu);
      else
          out = Inf;
      endif
  else
      out = n * log_sigma + sum (exp (-z) + z);
  endif
endfunction

%!demo
%! ## Sample 2 populations from 2 different exponential distibutions
%! rand ("seed", 1);   # for reproducibility
%! r1 = gevrnd (-0.5, 1, 2, 5000, 1);
%! rand ("seed", 2);   # for reproducibility
%! r2 = gevrnd (0, 1, -4, 5000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 50, 5);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! hold on
%!
%! ## Estimate their k, sigma, and mu parameters
%! k_sigma_muA = gevfit (r(:,1));
%! k_sigma_muB = gevfit (r(:,2));
%!
%! ## Plot their estimated PDFs
%! x = [-10:0.5:20];
%! y = gevpdf (x, k_sigma_muA(1), k_sigma_muA(2), k_sigma_muA(3));
%! plot (x, y, "-pr");
%! y = gevpdf (x, k_sigma_muB(1), k_sigma_muB(2), k_sigma_muB(3));
%! plot (x, y, "-sg");
%! ylim ([0, 0.7])
%! xlim ([-7, 5])
%! legend ({"Normalized HIST of sample 1 with ξ=-0.5, σ=1, μ=2", ...
%!          "Normalized HIST of sample 2 with ξ=0, σ=1, μ=-4",
%!     sprintf("PDF for sample 1 with estimated ξ=%0.2f, σ=%0.2f, μ=%0.2f", ...
%!                 k_sigma_muA(1), k_sigma_muA(2), k_sigma_muA(3)), ...
%!     sprintf("PDF for sample 3 with estimated ξ=%0.2f, σ=%0.2f, μ=%0.2f", ...
%!                 k_sigma_muB(1), k_sigma_muB(2), k_sigma_muB(3))})
%! title ("Two population samples from different exponential distibutions")
%! hold off

## Test output
%!test
%! x = 1:50;
%! [pfit, pci] = gevfit (x);
%! pfit_out = [-0.4407, 15.1923, 21.5309];
%! pci_out = [-0.7532, 11.5878, 16.5686; -0.1282, 19.9183, 26.4926];
%! assert (pfit, pfit_out, 1e-3);
%! assert (pci, pci_out, 1e-3);
%!test
%! x = 1:2:50;
%! [pfit, pci] = gevfit (x);
%! pfit_out = [-0.4434, 15.2024, 21.0532];
%! pci_out = [-0.8904, 10.3439, 14.0168; 0.0035, 22.3429, 28.0896];
%! assert (pfit, pfit_out, 1e-3);
%! assert (pci, pci_out, 1e-3);

## Test input validation
%!error<gevfit: X must be a vector.> gevfit (ones (2,5));
%!error<gevfit: wrong value for ALPHA.> gevfit ([1, 2, 3, 4, 5], 1.2);
%!error<gevfit: wrong value for ALPHA.> gevfit ([1, 2, 3, 4, 5], 0);
%!error<gevfit: wrong value for ALPHA.> gevfit ([1, 2, 3, 4, 5], "alpha");
%!error<gevfit: X and FREQ vectors mismatch.> ...
%! gevfit ([1, 2, 3, 4, 5], 0.05, [1, 2, 3, 2]);
%!error<gevfit: FREQ must not contain negative values.> ...
%! gevfit ([1, 2, 3, 4, 5], 0.05, [1, 2, 3, 2, -1]);
%!error<gevfit: FREQ must contain integer values.> ...
%! gevfit ([1, 2, 3, 4, 5], 0.05, [1, 2, 3, 2, 1.5]);
%!error<gevfit: 'options' argument must be a structure> ...
%! gevfit ([1, 2, 3, 4, 5], 0.05, struct ("option", 234));
%!error<gevfit: 'options' argument must be a structure> ...
%! gevfit ([1, 2, 3, 4, 5], 0.05, ones (1,5), struct ("option", 234));
