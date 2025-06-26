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
## @deftypefn  {statistics} {@var{paramhat} =} nakafit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} nakafit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} nakafit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} nakafit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} nakafit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} nakafit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate mean and confidence intervals for the Nakagami distribution.
##
## @code{@var{mu0} = nakafit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the Nakagami distribution given the data in
## @var{x}.  @qcode{@var{paramhat}(1)} is the shape parameter, @var{mu}, and
## @qcode{@var{paramhat}(2)} is the spread parameter, @var{omega}.
##
## @code{[@var{paramhat}, @var{paramci}] = nakafit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = nakafit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = nakafit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = nakafit (@var{params}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} must contain non-negative integer frequencies for the
## corresponding elements in @var{x}.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = nakafit (@dots{}, @var{options})} specifies control
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
## Further information about the Nakagami distribution can be found at
## @url{https://en.wikipedia.org/wiki/Nakagami_distribution}
##
## @seealso{nakacdf, nakainv, nakapdf, nakarnd, nakalike, nakastat}
## @end deftypefn

function [paramhat, paramci] = nakafit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("nakafit: X must be a vector.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("nakafit: wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("nakafit: X and CENSOR vectors mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("nakafit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("nakafit: FREQ must not contain negative values.");
  elseif (any (fix (freq) != freq))
    error ("nakafit: FREQ must contain integer values.");
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
      error (strcat (["nakafit: 'options' 5th argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  endif

  ## Expand frequency and censor vectors (if necessary)
  if (! all (freq == 1))
    xf = [];
    cf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
      cf = [cf, repmat(censor(i), 1, freq(i))];
    endfor
    x = xf;
    freq = ones (size (x));
    censor = cf;
  endif

  ## Get parameter estimates from the Gamma distribution
  paramhat = gamfit (x .^ 2, alpha, censor, freq, options);
  ## Transform back to Nakagami parameters
  paramhat(2) = paramhat(1) .* paramhat(2);

  ## Compute CIs using a log normal approximation for parameters.
  if (nargout > 1)
    ## Compute asymptotic covariance
    [~, acov] = nakalike (paramhat, x, censor, freq);
    ## Get standard errors
    se = sqrt (diag (acov))';
    ## Get normal quantiles
    probs = [alpha/2; 1-alpha/2];
    ## Compute muci using a normal approximation
    paramci(:,1) = norminv (probs, paramhat(1), se(1));
    ## Compute omegaci using a normal approximation for log (omega)
    paramci(:,2) = exp (norminv (probs, log (paramhat(2)), log (se(2))));
 endif

endfunction

%!demo
%! ## Sample 3 populations from different Nakagami distributions
%! randg ("seed", 5)  # for reproducibility
%! r1 = nakarnd (0.5, 1, 2000, 1);
%! randg ("seed", 2)   # for reproducibility
%! r2 = nakarnd (5, 1, 2000, 1);
%! randg ("seed", 7)   # for reproducibility
%! r3 = nakarnd (2, 2, 2000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0.05:0.1:3.5], 10);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 2.5]);
%! xlim ([0, 3.0]);
%! hold on
%!
%! ## Estimate their MU and LAMBDA parameters
%! mu_omegaA = nakafit (r(:,1));
%! mu_omegaB = nakafit (r(:,2));
%! mu_omegaC = nakafit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0.01:0.1:3.01];
%! y = nakapdf (x, mu_omegaA(1), mu_omegaA(2));
%! plot (x, y, "-pr");
%! y = nakapdf (x, mu_omegaB(1), mu_omegaB(2));
%! plot (x, y, "-sg");
%! y = nakapdf (x, mu_omegaC(1), mu_omegaC(2));
%! plot (x, y, "-^c");
%! legend ({"Normalized HIST of sample 1 with μ=0.5 and ω=1", ...
%!          "Normalized HIST of sample 2 with μ=5 and ω=1", ...
%!          "Normalized HIST of sample 3 with μ=2 and ω=2", ...
%!          sprintf("PDF for sample 1 with estimated μ=%0.2f and ω=%0.2f", ...
%!                  mu_omegaA(1), mu_omegaA(2)), ...
%!          sprintf("PDF for sample 2 with estimated μ=%0.2f and ω=%0.2f", ...
%!                  mu_omegaB(1), mu_omegaB(2)), ...
%!          sprintf("PDF for sample 3 with estimated μ=%0.2f and ω=%0.2f", ...
%!                  mu_omegaC(1), mu_omegaC(2))})
%! title ("Three population samples from different Nakagami distributions")
%! hold off

## Test output
%!test
%! paramhat = nakafit ([1:50]);
%! paramhat_out = [0.7355, 858.5];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = nakafit ([1:5]);
%! paramhat_out = [1.1740, 11];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = nakafit ([1:6], [], [], [1 1 1 1 1 0]);
%! paramhat_out = [1.1740, 11];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = nakafit ([1:5], [], [], [1 1 1 1 2]);
%! paramhat_out = nakafit ([1:5, 5]);
%! assert (paramhat, paramhat_out, 1e-4);

## Test input validation
%!error<nakafit: X must be a vector.> nakafit (ones (2,5));
%!error<nakafit: wrong value for ALPHA.> nakafit ([1, 2, 3, 4, 5], 1.2);
%!error<nakafit: wrong value for ALPHA.> nakafit ([1, 2, 3, 4, 5], 0);
%!error<nakafit: wrong value for ALPHA.> nakafit ([1, 2, 3, 4, 5], "alpha");
%!error<nakafit: X and CENSOR vectors mismatch.> ...
%! nakafit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<nakafit: X and CENSOR vectors mismatch.> ...
%! nakafit ([1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<nakafit: X and FREQ vectors mismatch.> ...
%! nakafit ([1, 2, 3, 4, 5], 0.05, zeros (1,5), [1 1 0]);
%!error<nakafit: X and FREQ vectors mismatch.> ...
%! nakafit ([1, 2, 3, 4, 5], [], [], [1 1 0 1 1]');
%!error<nakafit: FREQ must not contain negative values.> ...
%! nakafit ([1, 2, 3, 4, 5], [], [], [1 1 -1 1 1]);
%!error<nakafit: FREQ must contain integer values.> ...
%! nakafit ([1, 2, 3, 4, 5], [], [], [1 1 1.5 1 1]);
%!error<nakafit: 'options' 5th argument must be a structure> ...
%! nakafit ([1, 2, 3, 4, 5], 0.05, [], [], 2);
