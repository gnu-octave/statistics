## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{paramhat} =} ricefit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} ricefit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} ricefit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} ricefit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} ricefit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} ricefit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for the Rician distribution.
##
## @code{@var{paramhat} = ricefit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the Rician distribution given the data in
## @var{x}.  @qcode{@var{paramhat}(1)} is the non-centrality (distance)
## parameter, @math{s}, and @qcode{@var{paramhat}(2)} is the scale parameter,
## @math{sigma}.
##
## @code{[@var{paramhat}, @var{paramci}] = ricefit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = ricefit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = ricefit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = ricefit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = ricefit (@dots{}, @var{options})} specifies control
## parameters for the iterative algorithm used to compute the maximum likelihood
## estimates.  @var{options} is a structure with the following field and its
## default value:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.MaxFunEvals = 1000}
## @item @qcode{@var{options}.MaxIter = 500}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
##
## Further information about the Rician distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rice_distribution}
##
## @seealso{ricecdf, ricepdf, riceinv, ricernd, ricelike, ricestat}
## @end deftypefn

function [paramhat, paramci] = ricefit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("ricefit: X must be a vector.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("ricefit: wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("ricefit: X and CENSOR vectors mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("ricefit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("ricefit: FREQ cannot have negative values.");
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
      error (strcat (["ricefit: 'options' 5th argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  endif

  ## Get sample size and data type
  cls = class (x);
  szx = sum (freq);
  ncen = sum (freq .* censor);
  nunc = szx - ncen;

  ## Check for illegal value in X
  if (any (x <= 0))
    error ("ricefit: X must contain positive values.");
  endif

  ## Handle ill-conditioned cases: no data or all censored
  if (szx == 0 || nunc == 0 || any (! isfinite (x)))
    paramhat = nan (1, 2, cls);
    paramci = nan (2, cls);
    return
  endif

  ## Check for identical data in X
  if (! isscalar (x) && max (abs (diff (x)) ./ x(2:end)) <= sqrt (eps))
    paramhat = cast ([Inf, 0], cls);
    paramci = cast ([Inf, 0; Inf, 0], cls);
    return
  endif

  ## Use 2nd and 4th Moment Estimators of uncensored data as starting point
  xsq_uncensored = x(censor == 0) .^ 2;
  meanxsq = mean (xsq_uncensored);
  meanx_4 = mean (xsq_uncensored .^ 2);
  if (meanxsq ^ 2 < meanx_4 && meanx_4 < 2 * meanxsq ^ 2)
    nu_4 = 2 * meanxsq ^ 2 - meanx_4;
    nusq = sqrt (nu_4);
    sigmasq = 0.5 * (meanxsq - nusq);
    params = cast ([sqrt(nusq), sqrt(sigmasq)], cls);
  else
    params = cast ([1, 1], cls);
  endif

  ## Minimize negative log-likelihood to estimate parameters
  f = @(params) ricelike (params, x, censor, freq);
  [paramhat, ~, err, output] = fminsearch (f, params, options);
  ## Force positive parameter values
  paramhat = abs (paramhat);

  ## Handle errors
  if (err == 0)
    if (output.funcCount >= options.MaxFunEvals)
      warning (strcat (["ricefit: maximum number of function"], ...
                       [" evaluations are exceeded."]));
    elseif (output.iterations >= options.MaxIter)
      warning ("ricefit: maximum number of iterations are exceeded.");
    endif
  elseif (err < 0)
    error ("ricefit: no solution.");
  endif

  ## Compute CIs using a log normal approximation for parameters.
  if (nargout > 1)
    ## Compute asymptotic covariance
    [~, acov] = ricelike (paramhat, x, censor, freq);
    ## Get standard errors
    stderr = sqrt (diag (acov))';
    stderr = stderr ./ paramhat;
    ## Apply log transform
    phatlog = log (paramhat);
    ## Compute normal quantiles
    z = probit (alpha / 2);
    ## Compute CI
    paramci = [phatlog; phatlog] + [stderr; stderr] .* [z, z; -z, -z];
    ## Inverse log transform
    paramci = exp (paramci);
 endif

endfunction

%!demo
<<<<<<< Updated upstream
%! ## Sample 3 populations from different Gamma distibutions
=======
%! ## Sample 3 populations from different Rician distributions
>>>>>>> Stashed changes
%! randg ("seed", 5);    # for reproducibility
%! randp ("seed", 6);
%! r1 = ricernd (1, 2, 3000, 1);
%! randg ("seed", 2);    # for reproducibility
%! randp ("seed", 8);
%! r2 = ricernd (2, 4, 3000, 1);
%! randg ("seed", 7);    # for reproducibility
%! randp ("seed", 9);
%! r3 = ricernd (7.5, 1, 3000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 75, 4);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 0.7]);
%! xlim ([0, 12]);
%! hold on
%!
%! ## Estimate their α and β parameters
%! s_sigmaA = ricefit (r(:,1));
%! s_sigmaB = ricefit (r(:,2));
%! s_sigmaC = ricefit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0.01,0.1:0.2:18];
%! y = ricepdf (x, s_sigmaA(1), s_sigmaA(2));
%! plot (x, y, "-pr");
%! y = ricepdf (x, s_sigmaB(1), s_sigmaB(2));
%! plot (x, y, "-sg");
%! y = ricepdf (x, s_sigmaC(1), s_sigmaC(2));
%! plot (x, y, "-^c");
%! hold off
%! legend ({"Normalized HIST of sample 1 with s=1 and σ=2", ...
%!          "Normalized HIST of sample 2 with s=2 and σ=4", ...
%!          "Normalized HIST of sample 3 with s=7.5 and σ=1", ...
%!          sprintf("PDF for sample 1 with estimated s=%0.2f and σ=%0.2f", ...
%!                  s_sigmaA(1), s_sigmaA(2)), ...
%!          sprintf("PDF for sample 2 with estimated s=%0.2f and σ=%0.2f", ...
%!                  s_sigmaB(1), s_sigmaB(2)), ...
%!          sprintf("PDF for sample 3 with estimated s=%0.2f and σ=%0.2f", ...
%!                  s_sigmaC(1), s_sigmaC(2))})
%! title ("Three population samples from different Rician distibutions")
%! hold off

## Test output
%!test
%! [paramhat, paramci] = ricefit ([1:50]);
%! assert (paramhat, [15.3057, 17.6668], 1e-4);
%! assert (paramci, [9.5468, 11.7802; 24.5383, 26.4952], 1e-4);
%!test
%! [paramhat, paramci] = ricefit ([1:50], 0.01);
%! assert (paramhat, [15.3057, 17.6668], 1e-4);
%! assert (paramci, [8.2309, 10.3717; 28.4615, 30.0934], 1e-4);
%!test
%! [paramhat, paramci] = ricefit ([1:5]);
%! assert (paramhat, [2.3123, 1.6812], 1e-4);
%! assert (paramci, [1.0819, 0.6376; 4.9424, 4.4331], 1e-4);
%!test
%! [paramhat, paramci] = ricefit ([1:5], 0.01);
%! assert (paramhat, [2.3123, 1.6812], 1e-4);
%! assert (paramci, [0.8521, 0.4702; 6.2747, 6.0120], 1e-4);
%!test
%! freq = [1 1 1 1 5];
%! [paramhat, paramci] = ricefit ([1:5], [], [], freq);
%! assert (paramhat, [3.5181, 1.5565], 1e-4);
%! assert (paramci, [2.5893, 0.9049; 4.7801, 2.6772], 1e-4);
%!test
%! censor = [1 0 0 0 0];
%! [paramhat, paramci] = ricefit ([1:5], [], censor);
%! assert (paramhat, [3.2978, 1.1527], 1e-4);
%! assert (paramci, [2.3192, 0.5476; 4.6895, 2.4261], 1e-4);

## Test class of input preserved
%!assert (class (ricefit (single ([1:50]))), "single")

## Test input validation
%!error<ricefit: X must be a vector.> ricefit (ones (2))
%!error<ricefit: wrong value for ALPHA.> ricefit ([1:50], 1)
%!error<ricefit: wrong value for ALPHA.> ricefit ([1:50], -1)
%!error<ricefit: wrong value for ALPHA.> ricefit ([1:50], {0.05})
%!error<ricefit: wrong value for ALPHA.> ricefit ([1:50], "k")
%!error<ricefit: wrong value for ALPHA.> ricefit ([1:50], i)
%!error<ricefit: wrong value for ALPHA.> ricefit ([1:50], [0.01 0.02])
%!error<ricefit: X and CENSOR vectors mismatch.> ricefit ([1:50], [], [1 1])
%!error<ricefit: X and FREQ vectors mismatch.> ricefit ([1:50], [], [], [1 1])
%!error<ricefit: FREQ cannot have negative values.> ...
%! ricefit ([1:5], [], [], [1, 1, 2, 1, -1])
%!error<ricefit: X must contain positive values.> ricefit ([1 2 3 -4])
%!error<ricefit: X must contain positive values.> ricefit ([1 2 0], [], [1 0 0])
