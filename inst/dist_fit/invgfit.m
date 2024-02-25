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
## @deftypefn  {statistics} {@var{paramhat} =} invgfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} invgfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} invgfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} invgfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} invgfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} invgfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate mean and confidence intervals for the inverse Gaussian distribution.
##
## @code{@var{mu0} = invgfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the inverse Gaussian distribution given the
## data in @var{x}.  @qcode{@var{paramhat}(1)} is the scale parameter, @var{mu},
## and @qcode{@var{paramhat}(2)} is the shape parameter, @var{lambda}.
##
## @code{[@var{paramhat}, @var{paramci}] = invgfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = invgfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = invgfit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = invgfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = invgfit (@dots{}, @var{options})} specifies control
## parameters for the iterative algorithm used to compute ML estimates with the
## @code{fminsearch} function.  @var{options} is a structure with the following
## fields and their default values:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.MaxFunEvals = 400}
## @item @qcode{@var{options}.MaxIter = 200}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
##
## Further information about the inverse Gaussian distribution can be found at
## @url{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}
##
## @seealso{invgcdf, invginv, invgpdf, invgrnd, invglike, invgstat}
## @end deftypefn

function [paramhat, paramci] = invgfit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("invgfit: X must be a vector.");
  elseif (any (x <= 0))
    error ("invgfit: X must contain only positive values.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("invgfit: wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("invgfit: X and CENSOR vectors mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("invgfit: X and FREQ vectors mismatch.");
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
      error (strcat (["invgfit: 'options' 5th argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  endif

  n_censored = sum (freq .* censor);
  ## Compute parameters for uncensored data
  if (n_censored == 0)
    ## Expand frequency vector (MATLAB does not do this, in R2018 at least)
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    xbar = mean (xf);
    paramhat = [xbar, (1 ./ mean (1 ./ x - 1 ./ xbar))];
  else
    ## Use MLEs of the uncensored data as initial searching values
    x_uncensored = x(censor == 0);
    mu0 = mean (x_uncensored);
    lambda0 = 1 ./ mean (1 ./ x_uncensored - 1 ./ mu0);
    x0 = [mu0, lambda0];

    ## Minimize negative log-likelihood to estimate parameters
    f = @(params) invglike (params, x, censor, freq);
    [paramhat, ~, err, output] = fminsearch (f, x0, options);
    ## Force positive parameter values
    paramhat = abs (paramhat);

    ## Handle errors
    if (err == 0)
      if (output.funcCount >= options.MaxFunEvals)
        msg = "invgfit: maximum number of function evaluations are exceeded.";
        warning (msg);
      elseif (output.iterations >= options.MaxIter)
        warning ("invgfit: maximum number of iterations are exceeded.");
      endif
    elseif (err < 0)
      error ("invgfit: no solution.");
    endif
  endif

  ## Compute CIs using a log normal approximation for parameters.
  if (nargout > 1)
    ## Compute asymptotic covariance
    [~, acov] = invglike (paramhat, x, censor, freq);
    ## Get standard errors
    stderr = sqrt (diag (acov))';
    ## Get normal quantiles
    probs = [alpha/2; 1-alpha/2];
    ## Compute CI
    paramci = norminv([probs, probs], [paramhat; paramhat], [stderr; stderr]);
 endif

endfunction

%!demo
%! ## Sample 3 populations from different inverse Gaussian distibutions
%! rand ("seed", 5); randn ("seed", 5);   # for reproducibility
%! r1 = invgrnd (1, 0.2, 2000, 1);
%! rand ("seed", 2); randn ("seed", 2);   # for reproducibility
%! r2 = invgrnd (1, 3, 2000, 1);
%! rand ("seed", 7); randn ("seed", 7);   # for reproducibility
%! r3 = invgrnd (3, 1, 2000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0.1:0.1:3.2], 9);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 3]);
%! xlim ([0, 3]);
%! hold on
%!
%! ## Estimate their MU and LAMBDA parameters
%! mu_lambdaA = invgfit (r(:,1));
%! mu_lambdaB = invgfit (r(:,2));
%! mu_lambdaC = invgfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0:0.1:3];
%! y = invgpdf (x, mu_lambdaA(1), mu_lambdaA(2));
%! plot (x, y, "-pr");
%! y = invgpdf (x, mu_lambdaB(1), mu_lambdaB(2));
%! plot (x, y, "-sg");
%! y = invgpdf (x, mu_lambdaC(1), mu_lambdaC(2));
%! plot (x, y, "-^c");
%! hold off
%! legend ({"Normalized HIST of sample 1 with μ=1 and λ=0.5", ...
%!          "Normalized HIST of sample 2 with μ=2 and λ=0.3", ...
%!          "Normalized HIST of sample 3 with μ=4 and λ=0.5", ...
%!          sprintf("PDF for sample 1 with estimated μ=%0.2f and λ=%0.2f", ...
%!                  mu_lambdaA(1), mu_lambdaA(2)), ...
%!          sprintf("PDF for sample 2 with estimated μ=%0.2f and λ=%0.2f", ...
%!                  mu_lambdaB(1), mu_lambdaB(2)), ...
%!          sprintf("PDF for sample 3 with estimated μ=%0.2f and λ=%0.2f", ...
%!                  mu_lambdaC(1), mu_lambdaC(2))})
%! title ("Three population samples from different inverse Gaussian distibutions")
%! hold off

## Test output
%!test
%! paramhat = invgfit ([1:50]);
%! paramhat_out = [25.5, 19.6973];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = invgfit ([1:5]);
%! paramhat_out = [3, 8.1081];
%! assert (paramhat, paramhat_out, 1e-4);

## Test input validation
%!error<invgfit: X must be a vector.> invgfit (ones (2,5));
%!error<invgfit: X must contain only positive values.> invgfit ([-1 2 3 4]);
%!error<invgfit: wrong value for ALPHA.> invgfit ([1, 2, 3, 4, 5], 1.2);
%!error<invgfit: wrong value for ALPHA.> invgfit ([1, 2, 3, 4, 5], 0);
%!error<invgfit: wrong value for ALPHA.> invgfit ([1, 2, 3, 4, 5], "alpha");
%!error<invgfit: X and CENSOR vectors mismatch.> ...
%! invgfit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<invgfit: X and CENSOR vectors mismatch.> ...
%! invgfit ([1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<invgfit: X and FREQ vectors mismatch.> ...
%! invgfit ([1, 2, 3, 4, 5], 0.05, zeros (1,5), [1 1 0]);
%!error<invgfit: X and FREQ vectors mismatch.> ...
%! invgfit ([1, 2, 3, 4, 5], [], [], [1 1 0 1 1]');
%!error<invgfit: 'options' 5th argument must be a structure> ...
%! invgfit ([1, 2, 3, 4, 5], 0.05, [], [], 2);
