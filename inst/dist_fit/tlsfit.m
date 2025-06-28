## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{paramhat} =} tlsfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} tlsfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} tlsfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} tlsfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} tlsfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} tlsfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for the Location-scale Student's
## T distribution.
##
## @code{@var{muhat} = tlsfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the location-scale T distribution
## given the data in @var{x}.  @qcode{@var{paramhat}(1)} is the location
## parameter, @math{mu}, @qcode{@var{paramhat}(2)} is the scale parameter,
## @math{sigma}, and @qcode{@var{paramhat}(3)} is the degrees of freedom,
## @math{nu}.
##
## @code{[@var{paramhat}, @var{paramci}] = tlsfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = tlsfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = tlsfit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = tlsfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = tlsfit (@dots{}, @var{options})} specifies control
## parameters for the iterative algorithm used to compute ML estimates with the
## @code{fminsearch} function.  @var{options} is a structure with the following
## fields and their default values:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
##
## Further information about the location-scale Student's T distribution can be
## found at @url{https://en.wikipedia.org/wiki/Student%27s_t-distribution#Location-scale_t_distribution}
##
## @seealso{tlscdf, tlsinv, tlspdf, tlsrnd, tlslike, tlsstat}
## @end deftypefn

function [paramhat, paramci] = tlsfit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("tlsfit: X must be a vector.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("tlsfit: wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("tlsfit: X and CENSOR vectors mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("tlsfit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("tlsfit: FREQ cannot have negative values.");
  endif

  ## Get options structure or add defaults
  if (nargin < 5)
    options.Display = "off";
    options.TolX = 1e-6;
  else
    if (! isstruct (options) || ! isfield (options, "Display") || ...
                                ! isfield (options, "TolX"))
      error (strcat (["tlsfit: 'options' 5th argument must be a structure"], ...
                     [" with 'Display' and 'TolX' fields present."]));
    endif
  endif

  ## Starting points as robust estimators for MU and SIGMA, and method
  ## of moments for DF
  x_uncensored = x(censor==0);
  mu = median (x_uncensored);
  sigma = 1.253 * mad (x_uncensored);
  mom = kurtosis (x_uncensored);
  mom(mom < 4) = 4;
  nu = 2 * (2 * mom - 3) / (mom - 3);
  x0 = [mu, sigma, nu];

  ## Minimize negative log-likelihood to estimate parameters
  f = @(params) tlslike (params, x, censor, freq);
  [paramhat, ~, err, output] = fminsearch (f, x0, options);

  ## Handle errors
  if (err == 0)
    if (output.funcCount >= options.MaxFunEvals)
      warning (strcat (["tlsfit: maximum number of function"], ...
                       [" evaluations are exceeded."]));
    elseif (output.iterations >= options.MaxIter)
      warning ("tlsfit: maximum number of iterations are exceeded.");
    endif
  elseif (err < 0)
    error ("tlsfit: no solution.");
  endif

  ## Compute CIs using a log normal approximation for parameters.
  if (nargout > 1)
    ## Compute asymptotic covariance
    [~, acov] = tlslike (paramhat, x, censor, freq);
    ## Get standard errors
    stderr = sqrt (diag (acov))';
    mu_se = stderr(1);
    sigma_se = stderr(2) ./ paramhat(2);
    df_se = stderr(3) ./ paramhat(3);
    ## Apply log transform to SIGMA and DF
    sigma_log = log (paramhat(2));
    df_log = log (paramhat(3));
    ## Compute normal quantiles
    z = norminv (alpha / 2);
    ## Compute CI
    muci = [paramhat(1); paramhat(1)] + [mu_se; mu_se] .* [z; -z];
    sigmaci = [sigma_log; sigma_log] + [sigma_se; sigma_se] .* [z; -z];
    dfci = [df_log; df_log] + [df_se; df_se] .* [z; -z];
    ## Inverse log transform
    paramci = [muci, exp([sigmaci, dfci])];
 endif

endfunction

%!demo
%! ## Sample 3 populations from 3 different location-scale T distributions
%! randn ("seed", 1);    # for reproducibility
%! randg ("seed", 2);    # for reproducibility
%! r1 = tlsrnd (-4, 3, 1, 2000, 1);
%! randn ("seed", 3);    # for reproducibility
%! randg ("seed", 4);    # for reproducibility
%! r2 = tlsrnd (0, 3, 1, 2000, 1);
%! randn ("seed", 5);    # for reproducibility
%! randg ("seed", 6);    # for reproducibility
%! r3 = tlsrnd (5, 5, 4, 2000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [-21:21], [1, 1, 1]);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 0.25]);
%! xlim ([-20, 20]);
%! hold on
%!
%! ## Estimate their lambda parameter
%! mu_sigma_nuA = tlsfit (r(:,1));
%! mu_sigma_nuB = tlsfit (r(:,2));
%! mu_sigma_nuC = tlsfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [-20:0.1:20];
%! y = tlspdf (x, mu_sigma_nuA(1), mu_sigma_nuA(2), mu_sigma_nuA(3));
%! plot (x, y, "-pr");
%! y = tlspdf (x, mu_sigma_nuB(1), mu_sigma_nuB(2), mu_sigma_nuB(3));
%! plot (x, y, "-sg");
%! y = tlspdf (x, mu_sigma_nuC(1), mu_sigma_nuC(2), mu_sigma_nuC(3));
%! plot (x, y, "-^c");
%! hold off
%! legend ({"Normalized HIST of sample 1 with μ=0, σ=2 and nu=1", ...
%!          "Normalized HIST of sample 2 with μ=5, σ=2 and nu=1", ...
%!          "Normalized HIST of sample 3 with μ=3, σ=4 and nu=3", ...
%!          sprintf("PDF for sample 1 with estimated μ=%0.2f, σ=%0.2f, and ν=%0.2f", ...
%!                  mu_sigma_nuA(1), mu_sigma_nuA(2), mu_sigma_nuA(3)), ...
%!          sprintf("PDF for sample 2 with estimated μ=%0.2f, σ=%0.2f, and ν=%0.2f", ...
%!                  mu_sigma_nuB(1), mu_sigma_nuB(2), mu_sigma_nuB(3)), ...
%!          sprintf("PDF for sample 3 with estimated μ=%0.2f, σ=%0.2f, and ν=%0.2f", ...
%!                  mu_sigma_nuC(1), mu_sigma_nuC(2), mu_sigma_nuC(3))})
%! title ("Three population samples from different location-scale T distributions")
%! hold off

## Test output
%!test
%! x = [-1.2352, -0.2741, 0.1726, 7.4356, 1.0392, 16.4165];
%! [paramhat, paramci] = tlsfit (x);
%! paramhat_out = [0.035893, 0.862711, 0.649261];
%! paramci_out = [-0.949034, 0.154655, 0.181080; 1.02082, 4.812444, 2.327914];
%! assert (paramhat, paramhat_out, 1e-6);
%! assert (paramci, paramci_out, 1e-5);
%!test
%! x = [-1.2352, -0.2741, 0.1726, 7.4356, 1.0392, 16.4165];
%! [paramhat, paramci] = tlsfit (x, 0.01);
%! paramci_out = [-1.2585, 0.0901, 0.1212; 1.3303, 8.2591, 3.4771];
%! assert (paramci, paramci_out, 1e-4);

## Test input validation
%!error<tlsfit: X must be a vector.> tlsfit (ones (2,5));
%!error<tlsfit: wrong value for ALPHA.> tlsfit ([1, 2, 3, 4, 5], 1.2);
%!error<tlsfit: wrong value for ALPHA.> tlsfit ([1, 2, 3, 4, 5], 0);
%!error<tlsfit: wrong value for ALPHA.> tlsfit ([1, 2, 3, 4, 5], "alpha");
%!error<tlsfit: X and CENSOR vectors mismatch.> ...
%! tlsfit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<tlsfit: X and CENSOR vectors mismatch.> ...
%! tlsfit ([1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<tlsfit: X and FREQ vectors mismatch.> ...
%! tlsfit ([1, 2, 3, 4, 5], 0.05, zeros (1,5), [1 1 0]);
%!error<tlsfit: X and FREQ vectors mismatch.> ...
%! tlsfit ([1, 2, 3, 4, 5], [], [], [1 1 0 1 1]');
%!error<tlsfit: FREQ cannot have negative values.> ...
%! tlsfit ([1, 2, 3, 4, 5], [], [], [1 1 0 1 -1]);
%!error<tlsfit: 'options' 5th argument must be a structure> ...
%! tlsfit ([1, 2, 3, 4, 5], 0.05, [], [], 2);
