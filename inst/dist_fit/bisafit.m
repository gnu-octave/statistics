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
## @deftypefn  {statistics} {@var{paramhat} =} bisafit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} bisafit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} bisafit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} bisafit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} bisafit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} bisafit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate mean and confidence intervals for the Birnbaum-Saunders distribution.
##
## @code{@var{muhat} = bisafit (@var{x})} returns maximum likelihood estimates
## of the parameters of the Birnbaum-Saunders distribution given the data in
## @var{x}.  @qcode{@var{paramhat}(1)} is the scale parameter, @var{beta}, and
## @qcode{@var{paramhat}(2)} is the shape parameter, @var{gamma}.
##
## @code{[@var{paramhat}, @var{paramci}] = bisafit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = bisafit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = bisafit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = bisafit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = bisafit (@dots{}, @var{options})} specifies control
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
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bisacdf, bisainv, bisapdf, bisarnd, bisalike, bisastat}
## @end deftypefn

function [paramhat, paramci] = bisafit (x, alpha, censor, freq, options)

  ## Check input arguments
  if (! isvector (x))
    error ("bisafit: X must be a vector.");
  elseif (any (x <= 0))
    error ("bisafit: X must contain only positive values.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("bisafit: Wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("bisafit: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("bisafit: X and FREQ vector mismatch.");
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
      error (strcat (["bisafit: 'options' 5th argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  endif

  ## Starting points as suggested by Birnbaum and Saunders
  x_uncensored = x(censor==0);
  xubar = mean (x_uncensored);
  xuinv = mean (1 ./ x_uncensored);
  beta = sqrt (xubar ./ xuinv);
  gamma = 2 .* sqrt (sqrt (xubar .* xuinv) - 1);
  x0 = [beta, gamma];

  ## Minimize negative log-likelihood to estimate parameters
  f = @(params) bisalike (params, x, censor, freq);
  [paramhat, ~, err, output] = fminsearch (f, x0, options);
  ## Force positive parameter values
  paramhat = abs (paramhat);

  ## Handle errors
  if (err == 0)
    if (output.funcCount >= options.MaxFunEvals)
      warning ("bisafit: maximum number of function evaluations are exceeded.");
    elseif (output.iterations >= options.MaxIter)
      warning ("bisafit: maximum number of iterations are exceeded.");
    endif
  elseif (err < 0)
    error ("bisafit: NoSolution.");
  endif

  ## Compute CIs using a log normal approximation for parameters.
  if (nargout > 1)
    ## Compute asymptotic covariance
    [~, acov] = bisalike (paramhat, x, censor, freq);
    ## Get standard errors
    stderr = sqrt (diag (acov))';
    stderr = stderr ./ paramhat;
    ## Apply log transform
    phatlog = log (paramhat);
    ## Compute normal quantiles
    z = norminv (alpha / 2);
    ## Compute CI
    paramci = [phatlog; phatlog] + [stderr; stderr] .* [z, z; -z, -z];
    ## Inverse log transform
    paramci = exp (paramci);
 endif

endfunction

%!demo
%! ## Sample 3 populations from different Birnbaum-Saunders distibutions
%! rand ("seed", 5);    # for reproducibility
%! r1 = bisarnd (1, 0.5, 2000, 1);
%! rand ("seed", 2);    # for reproducibility
%! r2 = bisarnd (2, 0.3, 2000, 1);
%! rand ("seed", 7);    # for reproducibility
%! r3 = bisarnd (4, 0.5, 2000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 80, 4.2);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 1.1]);
%! xlim ([0, 8]);
%! hold on
%!
%! ## Estimate their α and β parameters
%! beta_gammaA = bisafit (r(:,1));
%! beta_gammaB = bisafit (r(:,2));
%! beta_gammaC = bisafit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0:0.1:8];
%! y = bisapdf (x, beta_gammaA(1), beta_gammaA(2));
%! plot (x, y, "-pr");
%! y = bisapdf (x, beta_gammaB(1), beta_gammaB(2));
%! plot (x, y, "-sg");
%! y = bisapdf (x, beta_gammaC(1), beta_gammaC(2));
%! plot (x, y, "-^c");
%! hold off
%! legend ({"Normalized HIST of sample 1 with β=1 and γ=0.5", ...
%!          "Normalized HIST of sample 2 with β=2 and γ=0.3", ...
%!          "Normalized HIST of sample 3 with β=4 and γ=0.5", ...
%!          sprintf("PDF for sample 1 with estimated β=%0.2f and γ=%0.2f", ...
%!                  beta_gammaA(1), beta_gammaA(2)), ...
%!          sprintf("PDF for sample 2 with estimated β=%0.2f and γ=%0.2f", ...
%!                  beta_gammaB(1), beta_gammaB(2)), ...
%!          sprintf("PDF for sample 3 with estimated β=%0.2f and γ=%0.2f", ...
%!                  beta_gammaC(1), beta_gammaC(2))})
%! title ("Three population samples from different Birnbaum-Saunders distibutions")
%! hold off

## Test output
%!test
%! paramhat = bisafit ([1:50]);
%! paramhat_out = [16.2649, 1.0156];
%! assert (paramhat, paramhat_out, 1e-4);
%!test
%! paramhat = bisafit ([1:5]);
%! paramhat_out = [2.5585, 0.5839];
%! assert (paramhat, paramhat_out, 1e-4);

## Test input validation
%!error<bisafit: X must be a vector.> bisafit (ones (2,5));
%!error<bisafit: X must contain only positive values.> bisafit ([-1 2 3 4]);
%!error<bisafit: Wrong value for ALPHA.> bisafit ([1, 2, 3, 4, 5], 1.2);
%!error<bisafit: Wrong value for ALPHA.> bisafit ([1, 2, 3, 4, 5], 0);
%!error<bisafit: Wrong value for ALPHA.> bisafit ([1, 2, 3, 4, 5], "alpha");
%!error<bisafit: X and CENSOR vector mismatch.> ...
%! bisafit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<bisafit: X and CENSOR vector mismatch.> ...
%! bisafit ([1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<bisafit: X and FREQ vector mismatch.> ...
%! bisafit ([1, 2, 3, 4, 5], 0.05, zeros (1,5), [1 1 0]);
%!error<bisafit: X and FREQ vector mismatch.> ...
%! bisafit ([1, 2, 3, 4, 5], [], [], [1 1 0 1 1]');
%!error<bisafit: 'options' 5th argument must be a structure> ...
%! bisafit ([1, 2, 3, 4, 5], 0.05, [], [], 2);
