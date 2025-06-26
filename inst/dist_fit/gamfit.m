## Copyright (C) 2019 Nir Krakauer <mail@nirkrakauer.net>
## Copyright (C) 2023-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Based on previous work by Martijn van Oosterhout <kleptog@svana.org>
## originally granted to the public domain.
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
## @deftypefn  {statistics} {@var{paramhat} =} gamfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gamfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gamfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} gamfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} gamfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} gamfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for the Gamma distribution.
##
## @code{@var{paramhat} = gamfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the Gamma distribution given the data in
## @var{x}.  @qcode{@var{paramhat}(1)} is the shape parameter, @var{a}, and
## @qcode{@var{paramhat}(2)} is the scale parameter, @var{b}.
##
## @code{[@var{paramhat}, @var{paramci}] = gamfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = gamfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = gamfit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = gamfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = gamfit (@dots{}, @var{options})} specifies control
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
## OCTAVE/MATLAB use the alternative parameterization given by the pair
## @math{α, β}, i.e. shape @var{a} and scale @var{b}.  In Wikipedia, the two
## common parameterizations use the pairs @math{k, θ}, as shape and scale, and
## @math{α, β}, as shape and rate, respectively.  The parameter names @var{a}
## and @var{b} used here (for MATLAB compatibility) correspond to the parameter
## notation @math{k, θ} instead of the @math{α, β} as reported in Wikipedia.
##
## Further information about the Gamma distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gamma_distribution}
##
## @seealso{gamcdf, gampdf, gaminv, gamrnd, gamlike}
## @end deftypefn

function [paramhat, paramci] = gamfit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("gamfit: X must be a vector.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("gamfit: wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("gamfit: X and CENSOR vectors mismatch.");
  endif

  ## Parse FREQ argument or add default
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("gamfit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("gamfit: FREQ must not contain negative values.");
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
      error (strcat (["gamfit: 'options' 5th argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  endif

  ## Remove zeros and NaNs from frequency vector (if necessary)
  if (! all (freq == 1))
    remove = freq == 0 | isnan (freq);
    x(remove) = [];
    censor(remove) = [];
    freq(remove) = [];
  endif

  ## Get sample size and data type
  cls = class (x);
  szx = sum (freq);
  ncen = sum (freq .* censor);
  nunc = szx - ncen;

  ## Check for illegal value in X
  if (ncen == 0 && any (x < 0))
    error ("gamfit: X cannot contain negative values.");
  endif
  if (ncen > 0 && any (x <= 0))
    error ("gamfit: X must contain positive values when censored.");
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

  ## When CENSOR and FREQ are default
  if (all (censor == 0) && all (freq == 1))
    ## Optimize with respect to log(a), since both A and B must be positive
    meanx = mean (x);
    x0 = 0;

    ## Minimize negative log-likelihood to estimate parameters
    f = @(loga) gamfit_search (loga, meanx, x);
    [loga, ~, err, output] = fminsearch (f, x0, options);

    ## Inverse log(a)
    a = exp (loga);
    b = meanx / a;
    paramhat = [a, b];

    ## Handle errors
    if (err == 0)
      if (output.funcCount >= options.MaxFunEvals)
        warning (strcat (["gamfit: maximum number of function"], ...
                         [" evaluations are exceeded."]));
      elseif (output.iterations >= options.MaxIter)
        warning ("gamfit: maximum number of iterations are exceeded.");
      endif
    elseif (err < 0)
      error ("gamfit: NoSolution.");
    endif
  endif

  ## No censoring
  if (all (censor == 0))

    ## Scale data to allow parameter estimation
    ## for extremely large or small values
    scale = sum (freq .* x) / szx;
    ## Check for all data being ~zero
    if (scale < realmin (cls))
      paramhat = cast ([NaN, 0], cls);
      paramci = cast ([NaN, 0; NaN, 0], cls);
      return
    endif
    scaledx = x / scale;

    ## Use Method of Moments for initial estimates
    meansqx = sum (freq .* (scaledx - 1) .^ 2) / szx;
    b = meansqx * szx / (szx - 1);
    a = 1 / b;

    ## Ensure that MLEs is possible, otherwise return initial estimates
    if (any (scaledx == 0))
      paramhat = [a, b*scale];
      paramci = nan (2, cls);
      warning ("gamfit: X contains zeros.");
      return

    ## Compute MLEs
    else

      ## Bracket the root of the scale parameter likelihood equation
      sumlogx = sum (freq .* log (scaledx));
      bracket = sumlogx / szx;
      if (lkeqn (a, bracket) > 0)
        upper = a;
        lower = 0.5 * upper;
        while (lkeqn (lower, bracket) > 0)
          upper = lower;
          lower = 0.5 * upper;
          if (lower < realmin (cls))
            error ("gamfit: no solution");
          endif
        endwhile
      else
        lower = a;
        upper = 2 * lower;
        while (lkeqn (upper, bracket) < 0)
          lower = upper;
          upper = 2 * lower;
          if (upper > realmax (cls))
            error ("gamfit: no solution");
          endif
        endwhile
      endif
      bounds = [lower upper];

      ## Find the root of the likelihood equation.
      opts = optimset ("fzero");
      opts = optimset (opts, "Display", "off");
      f = @(a) lkeqn (a, bracket);
      [a, lkeqnval, err] = fzero (f, bounds, opts);

      ## Rescale B
      paramhat = [a, (1/a)*scale];
    endif

  ## With censoring
  else

    ## Get uncensored data
    notc = ! censor;
    xunc = x(notc);
    freq_notc = freq(notc);

    ## Ensure that MLEs is possible and get initial estimates
    xuncbar = sum (freq_notc .* xunc) / nunc;
    s2unc = sum (freq_notc .* (xunc - xuncbar) .^ 2) / nunc;
    if s2unc <= 100.*eps(xuncbar.^2)

      ## When all uncensored observations are equal and greater than all
      ## the censored observations, the likelihood surface becomes infinite
      if (max (xunc) == max (x))
        paramhat = cast ([Inf, 0], cls);
        if (nunc > 1)
          paramci = cast ([Inf, 0; Inf, 0], cls);
        else
          paramci = cast ([0, 0; Inf, Inf], cls);
        endif
        return
      endif

      ## Set some default parameter estimates.
      x0 = [2, xuncbar./2];

    else
      ## Fit a Weibull distribution and equate the parameter estimates
      ## into a Gamma distribution
      wblphat = wblfit (x, alpha, censor, freq);
      [m, v] = wblstat (wblphat(1), wblphat(2));
      x0 = [m.*m./v, v./m];
    endif

    ## Minimize negative log-likelihood to estimate parameters
    f = @(params) gamlike (params, x, censor, freq);
    [paramhat, ~, err, output] = fminsearch (f, x0, options);
    ## Force positive parameter values
    paramhat = abs (paramhat);

    ## Handle errors
    if (err == 0)
      if (output.funcCount >= options.MaxFunEvals)
        warning (strcat (["gamfit: maximum number of function"], ...
                         [" evaluations are exceeded."]));
      elseif (output.iterations >= options.MaxIter)
        warning ("gamfit: maximum number of iterations are exceeded.");
      endif
    elseif (err < 0)
      error ("gamfit: no solution.");
    endif
  endif


  ## Compute CIs using a log normal approximation for parameters.
  if (nargout > 1)
    ## Compute asymptotic covariance
    [~, acov] = gamlike (paramhat, x, censor, freq);
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

## Helper function so we only have to minimize for one variable.
function nlogL = gamfit_search (loga, meanx, x)
  a = exp (loga);
  b = meanx / a;
  nlogL = gamlike ([a, b], x);
endfunction

## Helper function for MLE with no censoring
function v = lkeqn (a, bracket)
v = -bracket - log (a) + psi (a);
endfunction

%!demo
%! ## Sample 3 populations from different Gamma distributions
%! randg ("seed", 5);    # for reproducibility
%! r1 = gamrnd (1, 2, 2000, 1);
%! randg ("seed", 2);    # for reproducibility
%! r2 = gamrnd (2, 2, 2000, 1);
%! randg ("seed", 7);    # for reproducibility
%! r3 = gamrnd (7.5, 1, 2000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 75, 4);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 0.62]);
%! xlim ([0, 12]);
%! hold on
%!
%! ## Estimate their α and β parameters
%! a_bA = gamfit (r(:,1));
%! a_bB = gamfit (r(:,2));
%! a_bC = gamfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0.01,0.1:0.2:18];
%! y = gampdf (x, a_bA(1), a_bA(2));
%! plot (x, y, "-pr");
%! y = gampdf (x, a_bB(1), a_bB(2));
%! plot (x, y, "-sg");
%! y = gampdf (x, a_bC(1), a_bC(2));
%! plot (x, y, "-^c");
%! hold off
%! legend ({"Normalized HIST of sample 1 with α=1 and β=2", ...
%!          "Normalized HIST of sample 2 with α=2 and β=2", ...
%!          "Normalized HIST of sample 3 with α=7.5 and β=1", ...
%!          sprintf("PDF for sample 1 with estimated α=%0.2f and β=%0.2f", ...
%!                  a_bA(1), a_bA(2)), ...
%!          sprintf("PDF for sample 2 with estimated α=%0.2f and β=%0.2f", ...
%!                  a_bB(1), a_bB(2)), ...
%!          sprintf("PDF for sample 3 with estimated α=%0.2f and β=%0.2f", ...
%!                  a_bC(1), a_bC(2))})
%! title ("Three population samples from different Gamma distributions")
%! hold off

## Test output
%!shared x
%! x = [1.2 1.6 1.7 1.8 1.9 2.0 2.2 2.6 3.0 3.5 4.0 4.8 5.6 6.6 7.6];
%!test
%! [paramhat, paramci] = gamfit (x);
%! assert (paramhat, [3.4248, 0.9752], 1e-4);
%! assert (paramci, [1.7287, 0.4670; 6.7852, 2.0366], 1e-4);
%!test
%! [paramhat, paramci] = gamfit (x, 0.01);
%! assert (paramhat, [3.4248, 0.9752], 1e-4);
%! assert (paramci, [1.3945, 0.3705; 8.4113, 2.5668], 1e-4);
%!test
%! freq = [1 1 1 1 2 1 1 1 1 2 1 1 1 1 2];
%! [paramhat, paramci] = gamfit (x, [], [], freq);
%! assert (paramhat, [3.3025, 1.0615], 1e-4);
%! assert (paramci, [1.7710, 0.5415; 6.1584, 2.0806], 1e-4);
%!test
%! [paramhat, paramci] = gamfit (x, [], [], [1:15]);
%! assert (paramhat, [4.4484, 0.9689], 1e-4);
%! assert (paramci, [3.4848, 0.7482; 5.6785, 1.2546], 1e-4);
%!test
%! [paramhat, paramci] = gamfit (x, 0.01, [], [1:15]);
%! assert (paramhat, [4.4484, 0.9689], 1e-4);
%! assert (paramci, [3.2275, 0.6899; 6.1312, 1.3608], 1e-4);
%!test
%! cens = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
%! [paramhat, paramci] = gamfit (x, [], cens, [1:15]);
%! assert (paramhat, [4.7537, 0.9308], 1e-4);
%! assert (paramci, [3.7123, 0.7162; 6.0872, 1.2097], 1e-4);
%!test
%! cens = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
%! freq = [1 1 1 1 2 1 1 1 1 2 1 1 1 1 2];
%! [paramhat, paramci] = gamfit (x, [], cens, freq);
%! assert (paramhat, [3.4736, 1.0847], 1e-4);
%! assert (paramci, [1.8286, 0.5359; 6.5982, 2.1956], 1e-4);

## Test edge cases
%!test
%! [paramhat, paramci] = gamfit ([1 1 1 1 1 1]);
%! assert (paramhat, [Inf, 0]);
%! assert (paramci, [Inf, 0; Inf, 0]);
%!test
%! [paramhat, paramci] = gamfit ([1 1 1 1 1 1], [], [1 1 1 1 1 1]);
%! assert (paramhat, [NaN, NaN]);
%! assert (paramci, [NaN, NaN; NaN, NaN]);
%!test
%! [paramhat, paramci] = gamfit ([1 1 1 1 1 1], [], [], [1 1 1 1 1 1]);
%! assert (paramhat, [Inf, 0]);
%! assert (paramci, [Inf, 0; Inf, 0]);

## Test class of input preserved
%!assert (class (gamfit (single (x))), "single")

## Test input validation
%!error<gamfit: X must be a vector.> gamfit (ones (2))
%!error<gamfit: wrong value for ALPHA.> gamfit (x, 1)
%!error<gamfit: wrong value for ALPHA.> gamfit (x, -1)
%!error<gamfit: wrong value for ALPHA.> gamfit (x, {0.05})
%!error<gamfit: wrong value for ALPHA.> gamfit (x, "a")
%!error<gamfit: wrong value for ALPHA.> gamfit (x, i)
%!error<gamfit: wrong value for ALPHA.> gamfit (x, [0.01 0.02])
%!error<gamfit: X and FREQ vectors mismatch.>
%! gamfit ([1 2 3], 0.05, [], [1 5])
%!error<gamfit: FREQ must not contain negative values.>
%! gamfit ([1 2 3], 0.05, [], [1 5 -1])
%!error<gamfit: 'options' 5th argument must be a structure> ...
%! gamfit ([1:10], 0.05, [], [], 5)
%!error<gamfit: X cannot contain negative values.> gamfit ([1 2 3 -4])
%!error<gamfit: X must contain positive values when censored.> ...
%! gamfit ([1 2 0], [], [1 0 0])
