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
## @deftypefn  {statistics} {@var{paramhat} =} logifit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} logifit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} logifit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} logifit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} logifit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} logifit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate mean and confidence intervals for the logistic distribution.
##
## @code{@var{mu0} = logifit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the logistic distribution given the data in
## @var{x}.  @qcode{@var{paramhat}(1)} is the scale parameter, @var{mu}, and
## @qcode{@var{paramhat}(2)} is the shape parameter, @var{s}.
##
## @code{[@var{paramhat}, @var{paramci}] = logifit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = logifit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = logifit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = logifit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = logifit (@dots{}, @var{options})} specifies control
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
## Further information about the logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Logistic_distribution}
##
## @seealso{logicdf, logiinv, logipdf, logirnd, logilike, logistat}
## @end deftypefn

function [paramhat, paramci] = logifit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("logifit: X must be a vector.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("logifit: wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("logifit: X and CENSOR vectors mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("logifit: X and FREQ vectors mismatch.");
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
      error (strcat ("logifit: 'options' 5th argument must be a", ...
                     " structure with 'Display', 'MaxFunEvals',", ...
                     " 'MaxIter', and 'TolX' fields present."));
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

  ## Use MLEs of the uncensored data as initial searching values
  x_uncensored = x(censor == 0);
  mu0 = mean (x_uncensored);
  s0 = std (x_uncensored) .* sqrt (3) ./ pi;
  x0 = [mu0, s0];

  ## Minimize negative log-likelihood to estimate parameters
  f = @(params) logilike (params, x, censor, freq);
  [paramhat, ~, err, output] = fminsearch (f, x0, options);

  ## Handle errors
  if (err == 0)
    if (output.funcCount >= options.MaxFunEvals)
      msg = "logifit: maximum number of function evaluations are exceeded.";
      warning (msg);
    elseif (output.iterations >= options.MaxIter)
      warning ("logifit: maximum number of iterations are exceeded.");
    endif
  elseif (err < 0)
    error ("logifit: no solution.");
  endif

  ## Compute CIs using a log normal approximation for parameters.
  if (nargout > 1)
    ## Compute asymptotic covariance
    [~, acov] = logilike (paramhat, x, censor, freq);
    ## Get standard errors
    se = sqrt (diag (acov))';
    ## Get normal quantiles
    probs = [alpha/2; 1-alpha/2];
    ## Compute muci using a normal approximation
    paramci(:,1) = norminv (probs, paramhat(1), se(1));
    ## Compute sci using a normal approximation for log (s) and transform back
    paramci(:,2) = exp (norminv (probs, log (paramhat(2)), se(2) / paramhat(2)));
 endif

endfunction

%!demo
%! ## Sample 3 populations from different logistic distributions
%! rand ("seed", 5)  # for reproducibility
%! r1 = logirnd (2, 1, 2000, 1);
%! rand ("seed", 2)   # for reproducibility
%! r2 = logirnd (5, 2, 2000, 1);
%! rand ("seed", 7)   # for reproducibility
%! r3 = logirnd (9, 4, 2000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [-6:20], 1);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 0.3]);
%! xlim ([-5, 20]);
%! hold on
%!
%! ## Estimate their MU and LAMBDA parameters
%! mu_sA = logifit (r(:,1));
%! mu_sB = logifit (r(:,2));
%! mu_sC = logifit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [-5:0.5:20];
%! y = logipdf (x, mu_sA(1), mu_sA(2));
%! plot (x, y, "-pr");
%! y = logipdf (x, mu_sB(1), mu_sB(2));
%! plot (x, y, "-sg");
%! y = logipdf (x, mu_sC(1), mu_sC(2));
%! plot (x, y, "-^c");
%! hold off
%! legend ({"Normalized HIST of sample 1 with μ=1 and s=0.5", ...
%!          "Normalized HIST of sample 2 with μ=2 and s=0.3", ...
%!          "Normalized HIST of sample 3 with μ=4 and s=0.5", ...
%!          sprintf("PDF for sample 1 with estimated μ=%0.2f and s=%0.2f", ...
%!                  mu_sA(1), mu_sA(2)), ...
%!          sprintf("PDF for sample 2 with estimated μ=%0.2f and s=%0.2f", ...
%!                  mu_sB(1), mu_sB(2)), ...
%!          sprintf("PDF for sample 3 with estimated μ=%0.2f and s=%0.2f", ...
%!                  mu_sC(1), mu_sC(2))})
%! title ("Three population samples from different logistic distributions")
%! hold off

## Test output
%!test
%! paramhat = logifit ([1:50]);
%! paramhat_out = [25.5, 8.7724];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = logifit ([1:5]);
%! paramhat_out = [3, 0.8645];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = logifit ([1:6], [], [], [1 1 1 1 1 0]);
%! paramhat_out = [3, 0.8645];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = logifit ([1:5], [], [], [1 1 1 1 2]);
%! paramhat_out = logifit ([1:5, 5]);
%! assert (paramhat, paramhat_out, 1e-4);

## Test input validation
%!error<logifit: X must be a vector.> logifit (ones (2,5));
%!error<logifit: wrong value for ALPHA.> logifit ([1, 2, 3, 4, 5], 1.2);
%!error<logifit: wrong value for ALPHA.> logifit ([1, 2, 3, 4, 5], 0);
%!error<logifit: wrong value for ALPHA.> logifit ([1, 2, 3, 4, 5], "alpha");
%!error<logifit: X and CENSOR vectors mismatch.> ...
%! logifit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<logifit: X and CENSOR vectors mismatch.> ...
%! logifit ([1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<logifit: X and FREQ vectors mismatch.> ...
%! logifit ([1, 2, 3, 4, 5], 0.05, zeros (1,5), [1 1 0]);
%!error<logifit: X and FREQ vectors mismatch.> ...
%! logifit ([1, 2, 3, 4, 5], [], [], [1 1 0 1 1]');
%!error<logifit: 'options' 5th argument must be a structure> ...
%! logifit ([1, 2, 3, 4, 5], 0.05, [], [], 2);
