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
## @deftypefn  {statistics} {@var{paramhat} =} loglfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} loglfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} loglfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} loglfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} loglfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} loglfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate mean and confidence intervals for the log-logistic distribution.
##
## @code{@var{mu0} = loglfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the log-logistic distribution given the data
## in @var{x}.  @qcode{@var{paramhat}(1)} is the scale parameter, @var{a}, and
## @qcode{@var{paramhat}(2)} is the shape parameter, @var{b}.
##
## @code{[@var{paramhat}, @var{paramci}] = loglfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = loglfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = loglfit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = loglfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = loglfit (@dots{}, @var{options})} specifies control
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
## Further information about the log-logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-logistic_distribution}
##
## MATLAB compatibility: MATLAB uses an alternative parameterization given by
## the pair @math{μ, s}, i.e. @var{mu} and @var{s}, in analogy with the logistic
## distribution.  Their relation to the @var{a} and @var{b} parameters is given
## below:
##
## @itemize
## @item @qcode{@var{a} = exp (@var{mu})}
## @item @qcode{@var{b} = 1 / @var{s}}
## @end itemize
##
## @seealso{loglcdf, loglinv, loglpdf, loglrnd, logllike}
## @end deftypefn

function [paramhat, paramci] = loglfit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("loglfit: X must be a vector.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("loglfit: wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("loglfit: X and CENSOR vectors mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("loglfit: X and FREQ vectors mismatch.");
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
      error (strcat (["loglfit: 'options' 5th argument must be a"], ...
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

  ## Use MLEs of the uncensored data as initial searching values
  x_uncensored = x(censor == 0);
  a0 = mean (x_uncensored);
  b0 = 1 ./ (std (x_uncensored) .* sqrt (3) ./ pi);
  x0 = [a0, b0];

  ## Minimize negative log-likelihood to estimate parameters
  f = @(params) logllike (params, x, censor, freq);
  [paramhat, ~, err, output] = fminsearch (f, x0, options);
    ## Force positive parameter values
    paramhat = abs (paramhat);

  ## Handle errors
  if (err == 0)
    if (output.funcCount >= options.MaxFunEvals)
      msg = "loglfit: maximum number of function evaluations are exceeded.";
      warning (msg);
    elseif (output.iterations >= options.MaxIter)
      warning ("loglfit: maximum number of iterations are exceeded.");
    endif
  elseif (err < 0)
    error ("loglfit: no solution.");
  endif

  ## Compute CIs using a log normal approximation for parameters.
  if (nargout > 1)
    ## Compute asymptotic covariance
    [~, acov] = logllike (paramhat, x, censor, freq);
    ## Get standard errors
    se = sqrt (diag (acov))';
    ## Get normal quantiles
    probs = [alpha/2; 1-alpha/2];
    ## Compute muhat using a normal approximation
    paramci(:,1) = norminv (probs, paramhat(1), se(1));
    ## Compute shat using a normal approximation for log (s) and transform back
    paramci(:,2) = exp (norminv (probs, log (paramhat(2)), log (se(2))));
 endif

endfunction

%!demo
%! ## Sample 3 populations from different log-logistic distibutions
%! rand ("seed", 5)  # for reproducibility
%! r1 = loglrnd (1, 1, 2000, 1);
%! rand ("seed", 2)   # for reproducibility
%! r2 = loglrnd (1, 2, 2000, 1);
%! rand ("seed", 7)   # for reproducibility
%! r3 = loglrnd (1, 8, 2000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0.05:0.1:2.5], 10);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 3.5]);
%! xlim ([0, 2.0]);
%! hold on
%!
%! ## Estimate their MU and LAMBDA parameters
%! a_bA = loglfit (r(:,1));
%! a_bB = loglfit (r(:,2));
%! a_bC = loglfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0.01:0.1:2.01];
%! y = loglpdf (x, a_bA(1), a_bA(2));
%! plot (x, y, "-pr");
%! y = loglpdf (x, a_bB(1), a_bB(2));
%! plot (x, y, "-sg");
%! y = loglpdf (x, a_bC(1), a_bC(2));
%! plot (x, y, "-^c");
%! legend ({"Normalized HIST of sample 1 with α=1 and β=1", ...
%!          "Normalized HIST of sample 2 with α=1 and β=2", ...
%!          "Normalized HIST of sample 3 with α=1 and β=8", ...
%!          sprintf("PDF for sample 1 with estimated α=%0.2f and β=%0.2f", ...
%!                  a_bA(1), a_bA(2)), ...
%!          sprintf("PDF for sample 2 with estimated α=%0.2f and β=%0.2f", ...
%!                  a_bB(1), a_bB(2)), ...
%!          sprintf("PDF for sample 3 with estimated α=%0.2f and β=%0.2f", ...
%!                  a_bC(1), a_bC(2))})
%! title ("Three population samples from different log-logistic distibutions")
%! hold off

## Test output
%!test
%! paramhat = loglfit ([1:50]);
%! paramhat_out = [exp(3.097175), 1/0.468525];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = loglfit ([1:5]);
%! paramhat_out = [exp(1.01124), 1/0.336449];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = loglfit ([1:6], [], [], [1 1 1 1 1 0]);
%! paramhat_out = [exp(1.01124), 1/0.336449];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = loglfit ([1:5], [], [], [1 1 1 1 2]);
%! paramhat_out = loglfit ([1:5, 5]);
%! assert (paramhat, paramhat_out, 1e-4);

## Test input validation
%!error<loglfit: X must be a vector.> loglfit (ones (2,5));
%!error<loglfit: wrong value for ALPHA.> loglfit ([1, 2, 3, 4, 5], 1.2);
%!error<loglfit: wrong value for ALPHA.> loglfit ([1, 2, 3, 4, 5], 0);
%!error<loglfit: wrong value for ALPHA.> loglfit ([1, 2, 3, 4, 5], "alpha");
%!error<loglfit: X and CENSOR vectors mismatch.> ...
%! loglfit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<loglfit: X and CENSOR vectors mismatch.> ...
%! loglfit ([1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<loglfit: X and FREQ vectors mismatch.> ...
%! loglfit ([1, 2, 3, 4, 5], 0.05, zeros (1,5), [1 1 0]);
%!error<loglfit: X and FREQ vectors mismatch.> ...
%! loglfit ([1, 2, 3, 4, 5], [], [], [1 1 0 1 1]');
%!error<loglfit: 'options' 5th argument must be a structure> ...
%! loglfit ([1, 2, 3, 4, 5], 0.05, [], [], 2);
