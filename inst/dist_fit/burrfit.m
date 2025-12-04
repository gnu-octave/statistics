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
## @deftypefn  {statistics} {@var{paramhat} =} burrfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} burrfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} burrfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} burrfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} burrfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} burrfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate mean and confidence intervals for the Burr type XII distribution.
##
## @code{@var{muhat} = burrfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the Burr type XII distribution given the data
## in @var{x}.  @qcode{@var{paramhat}(1)} is the scale parameter, @var{lambda},
## @qcode{@var{paramhat}(2)} is the first shape parameter, @var{c}, and
## @qcode{@var{paramhat}(3)} is the second shape parameter, @var{k}
##
## @code{[@var{paramhat}, @var{paramci}] = burrfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = burrfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = burrfit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = burrfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = burrfit (@dots{}, @var{options})} specifies control
## parameters for the iterative algorithm used to compute the maximum likelihood
## estimates.  @var{options} is a structure with the following field and its
## default value:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.MaxFunEvals = 400}
## @item @qcode{@var{options}.MaxIter = 200}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
##
## Further information about the Burr type XII distribution can be found at
## @url{https://en.wikipedia.org/wiki/Burr_distribution}
##
## @seealso{burrcdf, burrinv, burrpdf, burrrnd, burrlike, burrstat}
## @end deftypefn

function [paramhat, paramci] = burrfit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("burrfit: X must be a vector.");
  elseif (any (x <= 0))
    error ("burrfit: X must contain only positive values.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("burrfit: wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("burrfit: X and CENSOR vectors mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("burrfit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("burrfit: FREQ must not contain negative values.");
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
      error (strcat ("burrfit: 'options' 5th argument must be a", ...
                     " structure with 'Display', 'MaxFunEvals',", ...
                     " 'MaxIter', and 'TolX' fields present."));
    endif
  endif

  ## Force censoring vector into logical
  notc = ! censor;
  cens = ! notc;

  ## Check for identical data in X
  if (! isscalar (x) && max (abs (diff (x)) ./ x(2:end)) <= sqrt (eps))
    warning ("burrfit: X must not contain identical data.");

    ## Return some sensical values for estimated parameters
    lambda = x(1);
    c = Inf;
    k = sum (notc .* freq) / sum (freq) / log (2);
    paramhat = [lambda, c, k];
    if (nargout > 1)
      paramci = [paramhat; paramhat];
    endif
    return
  endif

  ## Fit a Pareto distribution
  [paramhat_prt, nlogL_prt] = prtfit (x, cens, freq);
  ## Fit a Weibull distribution
  paramhat_wbl = wblfit (x, alpha, cens, freq);
  nlogL_wbl = wbllike (paramhat_wbl, x, cens, freq);

  ## Calculate the discriminator
  x_lambda = x ./ paramhat_wbl(1);
  x_lambdk = x_lambda .^ paramhat_wbl(2);
  discrimi = sum (freq .* (0.5 * x_lambdk .^ 2 - x_lambdk .* notc));

  ## Compute Burr distribution
  if (discrimi > 0)

    ## Expand data (if necessary)
    if (any (freq != 1))
      ## Preserve class
      x_expand = zeros (1, sum (freq), class (x));
      id0 = 1;
      for idx = 1:numel (x)
        x_expand(id0:id0 + freq(idx) - 1) = x_lambda(idx);
        id0 += freq(idx);
      endfor
    else
      x_expand = x_lambda;
    endif

    ## Calculate median and 3rd quartile to estimate LAMBDA and C parameters
    Q = prctile (x_expand);
    xl_median = Q(3);
    xl_upperq = Q(4);

    ## Avoid median and upper quartile being too close together
    IRDdist = sqrt (eps (xl_median)) * xl_median;
    if ((xl_upperq - xl_median) < IRDdist)
      if (any (x_lambda > xl_upperq))
        xl_upperq = min (x_lambda(x_lambda > xl_upperq));
      elseif (any (x_lambda < xl_upperq))
        xl_median = max (x_lambda(x_lambda < xl_median));
      endif
    endif

    ## Compute starting LAMBDA and C, either directly or by minimization
    if (xl_median >= xl_upperq / xl_median)
      l0 = xl_median;
      c0 = log(3)/log(xl_upperq/xl_median);
    else
      l0 = 1;
      opts = optimset ("fzero");
      opts = optimset (opts, "Display", "off");
      cmax  = log(realmax)/(2*log(xl_upperq/xl_median));
      c0 = fzero (@(c)(xl_upperq/xl_median).^c-xl_median.^c-2, [0, cmax], opts);
    endif
    ## Calculate starting K from other starting parameters and scaled data
    k0 = exp (compute_logk (x_lambda, l0, c0, censor, freq));

    ## Estimate parameters by minimizing the negative log-likelihood function
    f = @(params) burrlike (params, x_lambda, censor, freq);
    [paramhat, ~, err, output] = fminsearch (f, [l0, c0, k0], options);
    ## Force positive parameter values
    paramhat = abs (paramhat);

    ## Handle errors
    if (err == 0)
      if (output.funcCount >= options.MaxFunEvals)
        warning (strcat ("burrfit: maximum number of function", ...
                         " evaluations are exceeded."));
      elseif (output.iterations >= options.MaxIter)
        warning ("burrfit: maximum number of iterations are exceeded.");
      endif
    endif

    ## Scale back LAMBDA parameter
    paramhat(1) = paramhat(1) * paramhat_wbl(1);

    ## Compute negative log-likelihood with estimated parameters
    nlogL_burr = burrlike (paramhat, x, censor, freq);

    ## Check if fitting a Burr distribution is better than fitting a Pareto
    ## according to step 5 of the algorithmic implementation in Shao, 2004
    if (paramhat(3) > 1e-6 && nlogL_burr < nlogL_prt)
      ## Compute CIs using a log normal approximation for phat.
      if (nargout > 1)
        ## Compute asymptotic covariance
        [~, acov] = burrlike (paramhat, x, censor, freq);
        ## Get standard errors
        stderr = sqrt (diag (acov))';
        stderr = stderr ./ paramhat;
        ## Apply log transform
        phatlog = log (paramhat);
        ## Compute normal quantiles
        z = norminv (alpha / 2);
        ## Compute CI
        paramci = [phatlog; phatlog] + ...
                  [stderr; stderr] .* [z, z, z; -z, -z, -z];
        ## Inverse log transform
        paramci = exp (paramci);
      endif
    else
      if (nlogL_prt < nlogL_wbl)
        error ("burrfit: Pareto distribution fits better in X.");
      else
        error ("burrfit: Weibull distribution fits better in X.");
      endif
    endif
  else
    if (nlogL_prt < nlogL_wbl)
      error ("burrfit: Pareto distribution fits better in X.");
    else
      error ("burrfit: Weibull distribution fits better in X.");
    endif
  endif

endfunction

## Helper function for fitting a Pareto distribution
function [paramhat, nlogL] = prtfit (x, censor, freq)
  ## Force censoring vector into logical
  notc = ! censor;
  cens = ! notc;
  ## Compute MLE for x_m
  xm = x(notc);
  xm = min (xm);
  ## Handle case with all data censored
  if (all (cens))
    paramhat = [max(x), NaN];
    nlogL_prt = 0;
    return
  endif
  ## Compute some values
  logx = log (x);
  suml = sum (freq .* (logx - log (xm)) .* (x > xm));
  sumf = sum (freq .* notc);
  ## Compute MLE for alpha
  a = sumf ./ suml;
  ## Add MLEs to returning vector
  paramhat = [xm, a];
  ## Compute negative log-likelihood
  nlogL = a .* suml + sum (freq .* notc .* logx) - log (a) .* sumf;
endfunction

## Helper function for computing K from X, LAMBDA, and C
function logk = compute_logk (x, lambda, c, censor, freq)
  ## Force censoring vector into logical
  notc = ! censor;
  cens = ! notc;
  ## Precalculate some values
  xl = x ./ lambda;
  l1_xlc = log1p (xl .^ c);
  ## Avoid realmax overflow by approximation
  is_inf = isinf (l1_xlc);
  l1_xlc(is_inf) = c .* log (xl(is_inf));
  if (sum (freq .* l1_xlc) < eps)
    lsxc = log (sum (freq .* (x .^ c)));
    if (isinf (lsxc))
      [maxx, idx] = max (x);
      lsxc = c * log (freq(idx) * maxx);
    endif
    logk = log (sum (freq .* notc)) + c * log (lambda) - lsxc;
  else
    logk = log (sum (freq .* notc)) - log (sum (freq .* l1_xlc));
  endif
endfunction


%!demo
%! ## Sample 3 populations from different Burr type XII distributions
%! rand ("seed", 4);    # for reproducibility
%! r1 = burrrnd (3.5, 2, 2.5, 10000, 1);
%! rand ("seed", 2);    # for reproducibility
%! r2 = burrrnd (1, 3, 1, 10000, 1);
%! rand ("seed", 9);    # for reproducibility
%! r3 = burrrnd (0.5, 2, 3, 10000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0.1:0.2:20], [18, 5, 3]);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 3]);
%! xlim ([0, 5]);
%! hold on
%!
%! ## Estimate their α and β parameters
%! lambda_c_kA = burrfit (r(:,1));
%! lambda_c_kB = burrfit (r(:,2));
%! lambda_c_kC = burrfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0.01:0.15:15];
%! y = burrpdf (x, lambda_c_kA(1), lambda_c_kA(2), lambda_c_kA(3));
%! plot (x, y, "-pr");
%! y = burrpdf (x, lambda_c_kB(1), lambda_c_kB(2), lambda_c_kB(3));
%! plot (x, y, "-sg");
%! y = burrpdf (x, lambda_c_kC(1), lambda_c_kC(2), lambda_c_kC(3));
%! plot (x, y, "-^c");
%! hold off
%! legend ({"Normalized HIST of sample 1 with λ=3.5, c=2, and k=2.5", ...
%!          "Normalized HIST of sample 2 with λ=1, c=3, and k=1", ...
%!          "Normalized HIST of sample 3 with λ=0.5, c=2, and k=3", ...
%!  sprintf("PDF for sample 1 with estimated λ=%0.2f, c=%0.2f, and k=%0.2f", ...
%!          lambda_c_kA(1), lambda_c_kA(2), lambda_c_kA(3)), ...
%!  sprintf("PDF for sample 2 with estimated λ=%0.2f, c=%0.2f, and k=%0.2f", ...
%!          lambda_c_kB(1), lambda_c_kB(2), lambda_c_kB(3)), ...
%!  sprintf("PDF for sample 3 with estimated λ=%0.2f, c=%0.2f, and k=%0.2f", ...
%!          lambda_c_kC(1), lambda_c_kC(2), lambda_c_kC(3))})
%! title ("Three population samples from different Burr type XII distributions")
%! hold off

## Test output
%!test
%! l = 1; c = 2; k = 3;
%! r = burrrnd (l, c, k, 100000, 1);
%! lambda_c_kA = burrfit (r);
%! assert (lambda_c_kA(1), l, 0.2);
%! assert (lambda_c_kA(2), c, 0.2);
%! assert (lambda_c_kA(3), k, 0.3);
%!test
%! l = 0.5; c = 1; k = 3;
%! r = burrrnd (l, c, k, 100000, 1);
%! lambda_c_kA = burrfit (r);
%! assert (lambda_c_kA(1), l, 0.2);
%! assert (lambda_c_kA(2), c, 0.2);
%! assert (lambda_c_kA(3), k, 0.3);
%!test
%! l = 1; c = 3; k = 1;
%! r = burrrnd (l, c, k, 100000, 1);
%! lambda_c_kA = burrfit (r);
%! assert (lambda_c_kA(1), l, 0.2);
%! assert (lambda_c_kA(2), c, 0.2);
%! assert (lambda_c_kA(3), k, 0.3);
%!test
%! l = 3; c = 2; k = 1;
%! r = burrrnd (l, c, k, 100000, 1);
%! lambda_c_kA = burrfit (r);
%! assert (lambda_c_kA(1), l, 0.2);
%! assert (lambda_c_kA(2), c, 0.2);
%! assert (lambda_c_kA(3), k, 0.3);
%!test
%! l = 4; c = 2; k = 4;
%! r = burrrnd (l, c, k, 100000, 1);
%! lambda_c_kA = burrfit (r);
%! assert (lambda_c_kA(1), l, 0.2);
%! assert (lambda_c_kA(2), c, 0.2);
%! assert (lambda_c_kA(3), k, 0.3);

## Test input validation
%!error<burrfit: X must be a vector.> burrfit (ones (2,5));
%!error<burrfit: X must contain only positive values.> burrfit ([-1 2 3 4]);
%!error<burrfit: wrong value for ALPHA.> burrfit ([1, 2, 3, 4, 5], 1.2);
%!error<burrfit: wrong value for ALPHA.> burrfit ([1, 2, 3, 4, 5], 0);
%!error<burrfit: wrong value for ALPHA.> burrfit ([1, 2, 3, 4, 5], "alpha");
%!error<burrfit: X and CENSOR vectors mismatch.> ...
%! burrfit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<burrfit: X and CENSOR vectors mismatch.> ...
%! burrfit ([1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<burrfit: X and FREQ vectors mismatch.>
%! burrfit ([1, 2, 3, 4, 5], 0.05, [], [1, 1, 5])
%!error<burrfit: FREQ must not contain negative values.>
%! burrfit ([1, 2, 3, 4, 5], 0.05, [], [1, 5, 1, 1, -1])
%!error<burrfit: 'options' 5th argument must be a structure> ...
%! burrfit ([1:10], 0.05, [], [], 5)
