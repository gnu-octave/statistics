## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
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
## @deftypefn  {statistics} {@var{paramhat} =} evfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} evfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} evfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} evfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} evfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} evfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for extreme value distribution.
##
## @code{@var{paramhat} = evfit (@var{x})} returns maximum likelihood
## estimates of the parameters of the type 1 extreme value distribution (also
## known as the Gumbel distribution) given in @var{x}.  @var{paramhat(1)} is
## the location parameter, mu, and @var{paramhat(2)} is the scale parameter,
## sigma.
##
## @code{[@var{paramhat}, @var{paramci}] = evfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = evfit (@var{x}, @var{alpha})} returns 100(1-@var{alpha})
## percent confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = evfit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector of the same size as @var{x} with 1 for observations that
## are right-censored and 0 for observations that are observed exactly.
##
## @code{[@dots{}] = evfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but may contain any non-integer non-negative values.
##
## @code{[@dots{}] = evfit (@dots{}, @var{options})}
##
## @seealso{evcdf, evinv, evpdf, evrnd, evlike, evstat}
## @end deftypefn

function [paramhat, paramci] = evfit (x, alpha, censor, freq, options)

  ## Check for valid number of input arguments
  narginchk (1, 5);
  ## Check X for being a double precision vector
  if (! isvector (x) || length (x) < 2 || ! isa (x, "double"))
    error ("evfit: X must be a double-precision vector.");
  endif
  ## If X is a column vector, make it a row vector
  if (size (x, 1) > 1)
    x = x.';
  endif
  ## Check that X does not contain missing values (NaNs)
  if (any (isnan (x)))
    error ("evfit: X must NOT contain missing values (NaNs).");
  endif
  ## Parse extra input arguments or add defaults
  if (nargin > 1)
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("evfit: Wrong value for ALPHA.");
    endif
  else
    alpha = 0.05;
  endif
  if (nargin > 2)
    if (! isempty (censor) && ! all (size (censor) == size (x)))
      error ("evfit: Censoring vector must match X in size.");
    endif
  else
    censor = zeros (size (x));
  endif
  if (nargin > 3)
    if (! isempty (freq) && ! all (size (freq) == size (x)))
      error ("evfit: Frequency vector must match X in size.");
    endif
    ## Remove elements with zero frequency (if applicable)
    rm = find (freq == 0);
    if (length (rm) > 0)
      x(rm) = [];
      censor(rm) = [];
      freq(rm) = [];
    endif
  else
    freq = ones (size (x));
  endif
  ## Get options structure or add defaults
  if (nargin > 4)
    if (! isstruct (options) || ! isfield (options, "Display") || ...
                                ! isfield (options, "TolX"))
      error (strcat (["evfit: 'options' 5th argument must be a structure"], ...
                     [" with 'Display' and 'TolX' fields present."]));
    endif
  else
    options.Display = "off";
    options.TolX = 1e-6;
  endif
  ## Censor x and get number of samples
  sample_size = sum (freq);
  censored_sample_size = sum (freq .* censor);
  uncensored_sample_size = sample_size - censored_sample_size;
  x_range = range (x);
  x_max = max (x);
  ## Check cases that cannot make a fit.
  ## 1. All observations are censored
  if (sample_size == 0 || uncensored_sample_size == 0 || ! isfinite (x_range))
    paramhat = NaN (1, 2);
    paramci = NaN (2, 2);
    return
  endif
  ## 2. Constant x in X
  if (censored_sample_size == 0 && x_range == 0)
    paramhat = [x(1), 0];
    if (sample_size == 1)
      paramci = [-Inf, 0; Inf, Inf];
    else
      paramci = [paramhat, paramhat];
    endif
    return
  elseif (censored_sample_size == 0 && x_range != 0)
    ## Data can fit, so preprocess them to make likelihood eqn more stable.
    ## Shift x to max(x) == 0, min(x) = -1.
    x_0 = (x - x_max) ./ x_range;
    ## Get a rough initial estimate for scale parameter
    initial_sigma_parm = (sqrt (6) * std (x_0)) / pi;
    uncensored_weights = sum (freq .* x_0) ./ sample_size;
  endif
  ## 3. All uncensored observations are equal and greater than all censored ones
  uncensored_x_range = range (x(censor == 0));
  uncensored_x = x(censor == 0);
  if (censored_sample_size > 0 && uncensored_x_range == 0 ...
                               && uncensored_x(1) >= x_max)
    paramhat = [uncensored_x(1), 0];
    if uncensored_sample_size == 1
      paramci = [-Inf, 0; Inf, Inf];
    else
      paramci = [paramhat; paramhat];
    end
    return
  else
    ## Data can fit, so preprocess them to make likelihood eqn more stable.
    ## Shift x to max(x) == 0, min(x) = -1.
    x_0 = (x - x_max) ./ x_range;
    ## Get a rough initial estimate for scale parameter
    if (uncensored_x_range > 0)
      [F_y, y] = ecdf (x_0, "censoring", censor', "frequency", freq');
      pmid = (F_y(1:(end-1)) + F_y(2:end)) / 2;
      linefit = polyfit (log (- log (1 - pmid)), y(2:end), 1);
      initial_sigma_parm = linefit(1);
    else
      initial_sigma_parm = 1;
    endif
    uncensored_weights = sum (freq .* x_0 .* (1 - censor)) ./ ...
                              uncensored_sample_size;
  endif
  ## Find lower and upper boundaries for bracketing the likelihood equation for
  ## the extreme value scale parameter assessed in fzero function later on
  if (evscale_lkeq (initial_sigma_parm, x_0, freq, uncensored_weights) > 0)
    upper = initial_sigma_parm;
    lower = 0.5 * upper;
    while (evscale_lkeq (lower, x_0, freq, uncensored_weights) > 0)
      upper = lower;
      lower = 0.5 * upper;
      if (lower <= realmin ("double"))
        error ("evfit: no solution for maximum likelihood estimates.");
      endif
    endwhile
    boundaries = [lower, upper];
  else
    lower = initial_sigma_parm;
    upper = 2 * lower;
    while (evscale_lkeq (upper, x_0, freq, uncensored_weights) < 0)
      lower = upper;
      upper = 2 * lower;
      if (upper > realmax ("double"))
        error ("evfit: no solution for maximum likelihood estimates.");
      endif
    endwhile
    boundaries = [lower, upper];
  endif
  ## Compute maximum likelihood for scale parameter as the root of the equation
  ## Custom code for finding the value within the boundaries [lower, upper] that
  ## evscale_lkeq function returns zero
  ## First check that there is a root within the boundaries
  new_lo = boundaries(1);
  new_up = boundaries(2);
  v_lower = evscale_lkeq (new_lo, x_0, freq, uncensored_weights);
  v_upper = evscale_lkeq (new_up, x_0, freq, uncensored_weights);
  if (! (sign (v_lower) * sign (v_upper) <= 0))
    error ("evfit: no solution for maximum likelihood estimates.");
  endif
  ## Get a value at mid boundary range
  old_sigma = new_lo;
  new_sigma = (new_lo + new_up) / 2;
  new_fzero = evscale_lkeq (new_sigma, x_0, freq, uncensored_weights);
  ## Start searching
  cur_iter = 0;
  max_iter = 1e+3;
  while (cur_iter < max_iter && abs (old_sigma - new_sigma) > options.TolX)
    cur_iter++;
    if (new_fzero < 0)
      old_sigma = new_sigma;
      new_lo = new_sigma;
      new_sigma = (new_lo + new_up) / 2;
      new_fzero = evscale_lkeq (new_sigma, x_0, freq, uncensored_weights);
    else
      old_sigma = new_sigma;
      new_up = new_sigma;
      new_sigma = (new_lo + new_up) / 2;
      new_fzero = evscale_lkeq (new_sigma, x_0, freq, uncensored_weights);
    endif
  endwhile
  ## Check for maximum number of iterations
  if (cur_iter == max_iter)
    warning (strcat (["evfit: maximum number of function "], ...
                     [" evaluations (1e+4) has been reached."]));
  endif
  ## Compute mu
  muhat = new_sigma .* log (sum (freq .* exp (x_0 ./ new_sigma)) ./ ...
                                             uncensored_sample_size);
  ## Transform mu and sigma back to original location and scale
  paramhat = [(x_range*muhat)+x_max, x_range*new_sigma];
  ## Compute the CI for mu and sigma
  if (nargout == 2)
    probs = [alpha/2; 1-alpha/2];
    [~, acov] = evlike (paramhat, x, censor, freq);
    transfhat = [paramhat(1), log(paramhat(2))];
    se = sqrt (diag (acov))';
    se(2) = se(2) ./ paramhat(2);
    paramci = norminv ([probs, probs], [transfhat; transfhat], [se; se]);
    paramci(:,2) = exp (paramci(:,2));
  endif
endfunction

## Likelihood equation for the extreme value scale parameter.
function v = evscale_lkeq (sigma, x, freq, x_weighted_uncensored)
  freq = freq .* exp (x ./ sigma);
  v = sigma + x_weighted_uncensored - sum (x .* freq) / sum (freq);
endfunction

%!demo
%! ## Sample 3 populations from 3 different Extreme Value distibutions
%! r = [evrnd(2, 5, 200, 1), evrnd(-5, 3, 200, 1), evrnd(14, 8, 200, 1)];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 12, 1);
%! h = findobj(gca,'Type','patch');
%! set(h(1),'facecolor',"c");
%! set(h(2),'facecolor',"g");
%! set(h(3),'facecolor',"r");
%! hold on
%!
%! ## Estimate their lambda parameter
%! mu_sigmaA = evfit (r(:,1));
%! mu_sigmaB = evfit (r(:,2));
%! mu_sigmaC = evfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [min(r(:)):max(r(:))];
%! y = evpdf (x, mu_sigmaA(1), mu_sigmaA(2));
%! plot (x, y, "-pr");
%! y = evpdf (x, mu_sigmaB(1), mu_sigmaB(2));
%! plot (x, y, "-sg");
%! y = evpdf (x, mu_sigmaC(1), mu_sigmaC(2));
%! plot (x, y, "-^c");
%! hold off
%! legend ({"Normalized HIST of sample 1 with μ=2 and σ=5", ...
%!          "Normalized HIST of sample 2 with μ=-5 and σ=3", ...
%!          "Normalized HIST of sample 3 with μ=14 and σ=8", ...
%!          sprintf("PDF for sample 1 with estimated μ=%0.2f and σ=%0.2f", ...
%!                  mu_sigmaA(1), mu_sigmaA(2)), ...
%!          sprintf("PDF for sample 2 with estimated μ=%0.2f and σ=%0.2f", ...
%!                  mu_sigmaB(1), mu_sigmaB(2)), ...
%!          sprintf("PDF for sample 3 with estimated μ=%0.2f and σ=%0.2f", ...
%!                  mu_sigmaC(1), mu_sigmaC(2))})
%! title ("Three population samples from different Extreme Value distibutions")
%! hold off

## test results
%!test
%! x = 1:50;
%! [paramhat, paramci] = evfit (x);
%! paramhat_out = [32.6811, 13.0509];
%! paramci_out = [28.8504, 10.5294; 36.5118, 16.1763];
%! assert (paramhat, paramhat_out, 1e-4);
%! assert (paramci, paramci_out, 1e-4);
%!test
%! x = 1:50;
%! [paramhat, paramci] = evfit (x, 0.01);
%! paramci_out = [27.6468, 9.8426; 37.7155, 17.3051];
%! assert (paramci, paramci_out, 1e-4);

## test input validation
%!error<evfit: X must be a double-precision vector.> evfit (ones (2,5));
%!error<evfit: X must be a double-precision vector.> evfit (single (ones (1,5)));
%!error<evfit: X must NOT contain missing values> evfit ([1, 2, 3, 4, NaN]);
%!error<evfit: Wrong value for ALPHA.> evfit ([1, 2, 3, 4, 5], 1.2);
%!error<evfit: Censoring vector must> evfit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<evfit: Frequency vector must> evfit ([1, 2, 3, 4, 5], 0.05, [], [1 1 0]);
%!error<evfit: 'options' 5th argument> evfit ([1, 2, 3, 4, 5], 0.05, [], [], 2);
