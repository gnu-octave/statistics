## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{paramhat} =} nbinfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} nbinfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} nbinfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} nbinfit (@var{x}, @var{alpha}, @var{options})
##
## Estimate parameter and confidence intervals for the negative binomial
## distribution.
##
## @code{@var{paramhat} = nbinfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the negative binomial distribution given the
## data in vector @var{x}.  @qcode{@var{paramhat}(1)} is the number of successes
## until the experiment is stopped, @var{r}, and @qcode{@var{paramhat}(2)} is
## the probability of success in each experiment, @var{ps}.
##
## @code{[@var{paramhat}, @var{paramci}] = nbinfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@var{paramhat}, @var{paramci}] = nbinfit (@var{x}, @var{alpha})} also
## returns the @qcode{100 * (1 - @var{alpha})} percent confidence intervals of
## the estimated parameter.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.
##
## @code{[@var{paramhat}, @var{paramci}] = nbinfit (@var{x}, @var{alpha},
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
## When @var{r} is an integer, the negative binomial distribution is also known
## as the Pascal distribution and it models the number of failures in @var{x}
## before a specified number of successes is reached in a series of independent,
## identical trials.  Its parameters are the probability of success in a single
## trial, @var{ps}, and the number of successes, @var{r}.  A special case of the
## negative binomial distribution, when @qcode{@var{r} = 1}, is the geometric
## distribution, which models the number of failures before the first success.
##
## @var{r} can also have non-integer positive values, in which form the negative
## binomial distribution, also known as the Polya distribution, has no
## interpretation in terms of repeated trials, but, like the Poisson
## distribution, it is useful in modeling count data.  The negative binomial
## distribution is more general than the Poisson distribution because it has a
## variance that is greater than its mean, making it suitable for count data
## that do not meet the assumptions of the Poisson distribution.  In the limit,
## as @var{r} increases to infinity, the negative binomial distribution
## approaches the Poisson distribution.
##
## Further information about the negative binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Negative_binomial_distribution}
##
## @seealso{nbincdf, nbininv, nbinpdf, nbinrnd, nbinlike, nbinstat}
## @end deftypefn

function [paramhat, paramci] = nbinfit (x, alpha, options)

  ## Check data in X
  if (any (x < 0))
    error ("nbinfit: X cannot have negative values.");
  endif
  if (! isvector (x))
    error ("nbinfit: X must be a vector.");
  endif
  if (any (x < 0) || any (x != round (x)) || any (isinf (x)))
    error ("nbinfit: X must be a non-negative integer.");
  endif

  ## Check ALPHA
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("nbinfit: wrong value for ALPHA.");
  endif

  ## Get options structure or add defaults
  if (nargin < 3)
    options.Display = "off";
    options.MaxFunEvals = 400;
    options.MaxIter = 200;
    options.TolX = 1e-6;
  else
    if (! isstruct (options) || ! isfield (options, "Display") ||
        ! isfield (options, "MaxFunEvals") || ! isfield (options, "MaxIter")
                                           || ! isfield (options, "TolX"))
      error (strcat (["nbinfit: 'options' 3rd argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  endif

  ## Ensure that a negative binomial fit is valid.
  xbar = mean (x);
  varx = var (x);
  if (varx <= xbar)
    paramhat = cast ([Inf, 1.0], class (x));
    paramci = cast ([Inf, 1; Inf, 1], class (x));
    fprintf ("warning: nbinfit: mean exceeds variance.\n");
    return
  endif

  ## Use Method of Moments estimates as starting point for MLEs.
  rhat = (xbar .* xbar) ./ (varx - xbar);

  ## Minimize negative log-likelihood to estimate parameters by parameterizing
  ## with mu=r(1-p)/p, so it becomes 1-parameter search for rhat.
  f = @(rhat) nbinfit_search (rhat, x, numel (x), sum (x), options.TolX);
  [rhat, ~, err, output] = fminsearch (f, rhat, options);

  ## Handle errors
  if (err == 0)
    if (output.funcCount >= options.MaxFunEvals)
      warning (strcat (["nbinfit: maximum number of function"], ...
                       [" evaluations are exceeded."]));
    elseif (output.iterations >= options.MaxIter)
      warning ("nbinfit: maximum number of iterations are exceeded.");
    endif
  elseif (err < 0)
    error ("nbinfit: NoSolution.");
  endif

  ## Compute parameter estimates
  pshat = rhat ./ (xbar + rhat);
  paramhat = [rhat, pshat];

  ## Compute confidence interval
  if (nargout > 1)
    [~, avar] = nbinlike (paramhat, x);
    ## Get standard errors
    sigma = sqrt (diag (avar));
    ## Get normal quantiles
    probs = [alpha/2; 1-alpha/2];
    ## Compute paramci using a normal approximation
    paramci = norminv ([probs, probs], [paramhat; paramhat], [sigma'; sigma']);
    ## Restrict CI to valid values: r >= 0, 0 <= ps <= 1
    paramci(paramci < 0) = 0;
    if (paramci(2,2) > 1)
      paramci(2,2) = 1;
    endif
  endif

endfunction

## Helper function for minimizing the negative log-likelihood
function nll = nbinfit_search (r, x, nx, sx, tol)
  if (r < tol)
    nll = Inf;
  else
    xbar = sx / nx;
    nll = -sum (gammaln (r +x )) + nx * gammaln (r) ...
          -nx * r * log (r / (xbar + r)) - sx * log (xbar / (xbar + r));
  endif
endfunction

%!demo
%! ## Sample 2 populations from different negative binomial distibutions
%! randp ("seed", 5); randg ("seed", 5);    # for reproducibility
%! r1 = nbinrnd (2, 0.15, 5000, 1);
%! randp ("seed", 8); randg ("seed", 8);    # for reproducibility
%! r2 = nbinrnd (5, 0.2, 5000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0:51], 1);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! hold on
%!
%! ## Estimate their probability of success
%! r_psA = nbinfit (r(:,1));
%! r_psB = nbinfit (r(:,2));
%!
%! ## Plot their estimated PDFs
%! x = [0:40];
%! y = nbinpdf (x, r_psA(1), r_psA(2));
%! plot (x, y, "-pg");
%! x = [min(r(:,2)):max(r(:,2))];
%! y = nbinpdf (x, r_psB(1), r_psB(2));
%! plot (x, y, "-sc");
%! ylim ([0, 0.1])
%! xlim ([0, 50])
%! legend ({"Normalized HIST of sample 1 with r=2 and ps=0.15", ...
%!          "Normalized HIST of sample 2 with r=5 and ps=0.2", ...
%!          sprintf("PDF for sample 1 with estimated r=%0.2f and ps=%0.2f", ...
%!                  r_psA(1), r_psA(2)), ...
%!          sprintf("PDF for sample 2 with estimated r=%0.2f and ps=%0.2f", ...
%!                  r_psB(1), r_psB(2))})
%! title ("Two population samples from negative different binomial distibutions")
%! hold off

## Test output
%!test
%! [paramhat, paramci] = nbinfit ([1:50]);
%! assert (paramhat, [2.420857, 0.086704], 1e-6);
%! assert (paramci(:,1), [1.382702; 3.459012], 1e-6);
%! assert (paramci(:,2), [0.049676; 0.123732], 1e-6);
%!test
%! [paramhat, paramci] = nbinfit ([1:20]);
%! assert (paramhat, [3.588233, 0.254697], 1e-6);
%! assert (paramci(:,1), [0.451693; 6.724774], 1e-6);
%! assert (paramci(:,2), [0.081143; 0.428251], 1e-6);
%!test
%! [paramhat, paramci] = nbinfit ([1:10]);
%! assert (paramhat, [8.8067, 0.6156], 1e-4);
%! assert (paramci(:,1), [0; 30.7068], 1e-4);
%! assert (paramci(:,2), [0.0217; 1], 1e-4);

## Test input validation
%!error<nbinfit: X cannot have negative values.> nbinfit ([-1 2 3 3])
%!error<nbinfit: X must be a vector.> nbinfit (ones (2))
%!error<nbinfit: X must be a non-negative integer.> nbinfit ([1 2 1.2 3])
%!error<nbinfit: wrong value for ALPHA.> nbinfit ([1 2 3], 0)
%!error<nbinfit: wrong value for ALPHA.> nbinfit ([1 2 3], 1.2)
%!error<nbinfit: wrong value for ALPHA.> nbinfit ([1 2 3], [0.02 0.05])
%!error<nbinfit: 'options' 3rd argument must be a structure> ...
%! nbinfit ([1, 2, 3, 4, 5], 0.05, 2);
