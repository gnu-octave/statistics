## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{paramhat} =} lognfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} lognfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} lognfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} lognfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} lognfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} lognfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for the log-normal distribution.
##
## @code{@var{paramhat} = lognfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the log-normal distribution given the data in
## vector @var{x}.  @qcode{@var{paramhat}([1, 2])} corresponds to the mean and
## standard deviation, respectively, of the associated normal distribution.
##
## @code{[@var{paramhat}, @var{paramci}] = lognfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @itemize
## @item
## @var{alpha} is a scalar value in the range @math{(0,1)} specifying the
## confidence level for the confidence intervals calculated as
## @qcode{100 * (1 - @var{alpha})} percent.  By default, the optional argument
## @var{alpha} is 0.05 corresponding to 95% confidence intervals.  Pass in
## @qcode{[]} for @var{alpha} to use the default values.
##
## @item
## @var{censor} is a logical vector of the same length as @var{x} specifying
## whether each value in @var{x} is right-censored or not.  1 indicates
## observations that are right-censored and 0 indicates observations that are
## fully observed.  With censoring, @var{muhat} and @var{sigmahat} are the
## maximum likelihood estimates (MLEs).  If empty, the default is an array of
## 0s, meaning that all observations are fully observed.
##
## @item
## @var{freq} is a vector of the same length as @var{x} and it typically
## contains non-negative integer counts of the corresponding elements in
## @var{x}, but it may contain any non-integer non-negative values.  By default,
## or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @item
## @var{options} is a structure with the control parameters for
## @code{fminsearch} which is used internally to compute MLEs for censored data.
## By default, it uses the following options:
## @itemize
## @item @qcode{@var{options}.Display = "off"}
## @item @qcode{@var{options}.MaxFunEvals = 400}
## @item @qcode{@var{options}.MaxIter = 200}
## @item @qcode{@var{options}.TolX = 1e-6}
## @end itemize
## @end itemize
##
## With no censor, the estimate of the standard deviation,
## @qcode{@var{paramhat}(2)}, is the square root of the unbiased estimate of the
## variance of @qcode{log (@var{x})}.  With censored data, the maximum
## likelihood estimate is returned.
##
## @seealso{logncdf, logninv, lognpdf, lognrnd, lognlike, lognstat}
## @end deftypefn

function [paramhat, paramci] = lognfit (x, alpha, censor, freq, options)

  ## Check X for valid data
  if (! isvector (x) || ! isnumeric (x) || any (x <= 0))
    error ("lognfit: X must be a numeric vector of positive values.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("lognfit: Wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = [];
  elseif (! isequal (size (x), size (censor)))
    error ("lognfit: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = [];
  elseif (! isequal (size(x), size(freq)))
    error ("lognfit: X and FREQ vector mismatch.");
  endif

  ## Check options structure or add defaults
  if (nargin > 4 && ! isempty (options))
    if (! isstruct (options) || ! isfield (options, "Display") ||
        ! isfield (options, "MaxFunEvals") || ! isfield (options, "MaxIter")
                                           || ! isfield (options, "TolX"))
      error (strcat (["lognfit: 'options' 5th argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  else
    options = [];
  endif

  ## Fit a normal distribution to the logged data
  if (nargout <= 1)
    [muhat, sigmahat] = normfit (log (x), alpha, censor, freq, options);
    paramhat = [muhat, sigmahat];
  else
    [muhat, sigmahat, muci, sigmaci] = normfit (log (x), alpha, ...
                                                censor, freq, options);
    paramhat = [muhat, sigmahat];
    paramci = [muci, sigmaci];
  endif

endfunction


%!demo
%! ## Sample 3 populations from 3 different log-normal distibutions
%! randn ("seed", 1);    # for reproducibility
%! r1 = lognrnd(0, 0.25, 1000, 1);
%! randn ("seed", 2);    # for reproducibility
%! r2 = lognrnd(0, 0.5, 1000, 1);
%! randn ("seed", 3);    # for reproducibility
%! r3 = lognrnd(0, 1, 1000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 30, 2);
%! h = findobj(gca,'Type','patch');
%! set(h(1),'facecolor',"c");
%! set(h(2),'facecolor',"g");
%! set(h(3),'facecolor',"r");
%! hold on
%!
%! ## Estimate their mu and sigma parameters
%! mu_sigmaA = lognfit (r(:,1));
%! mu_sigmaB = lognfit (r(:,2));
%! mu_sigmaC = lognfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0:0.1:6];
%! y = lognpdf (x, mu_sigmaA(1), mu_sigmaA(2));
%! plot (x, y, "-pr");
%! y = lognpdf (x, mu_sigmaB(1), mu_sigmaB(2));
%! plot (x, y, "-sg");
%! y = lognpdf (x, mu_sigmaC(1), mu_sigmaC(2));
%! plot (x, y, "-^c");
%! ylim ([0, 2])
%! xlim ([0, 6])
%! hold off
%! legend ({"Normalized HIST of sample 1 with mu=0, σ=0.25", ...
%!          "Normalized HIST of sample 2 with mu=0, σ=0.5", ...
%!          "Normalized HIST of sample 3 with mu=0, σ=1", ...
%!          sprintf("PDF for sample 1 with estimated mu=%0.2f and σ=%0.2f", ...
%!                  mu_sigmaA(1), mu_sigmaA(2)), ...
%!          sprintf("PDF for sample 2 with estimated mu=%0.2f and σ=%0.2f", ...
%!                  mu_sigmaB(1), mu_sigmaB(2)), ...
%!          sprintf("PDF for sample 3 with estimated mu=%0.2f and σ=%0.2f", ...
%!                  mu_sigmaC(1), mu_sigmaC(2))}, "location", "northeast")
%! title ("Three population samples from different log-normal distibutions")
%! hold off

## Test output
%!test
%! x = lognrnd (3, 5, [1000, 1]);
%! [paramhat, paramci] = lognfit (x, 0.01);
%! assert (paramci(1,1) < 3);
%! assert (paramci(1,2) > 3);
%! assert (paramci(2,1) < 5);
%! assert (paramci(2,2) > 5);

## Test input validation
%!error<lognfit: X must be a numeric vector of positive values.> ...
%! lognfit (ones (20,3))
%!error<lognfit: X must be a numeric vector of positive values.> ...
%! lognfit ({1, 2, 3, 4, 5})
%!error<lognfit: X must be a numeric vector of positive values.> ...
%! lognfit ([-1, 2, 3, 4, 5])
%!error<lognfit: Wrong value for ALPHA.> ...
%! lognfit (ones (20,1), 0)
%!error<lognfit: Wrong value for ALPHA.> ...
%! lognfit (ones (20,1), -0.3)
%!error<lognfit: Wrong value for ALPHA.> ...
%! lognfit (ones (20,1), 1.2)
%!error<lognfit: Wrong value for ALPHA.> ...
%! lognfit (ones (20,1), [0.05,  0.1])
%!error<lognfit: Wrong value for ALPHA.> ...
%! lognfit (ones (20,1), 0.02+i)
%!error<lognfit: X and CENSOR vector mismatch.> ...
%! lognfit (ones (20,1), [], zeros(15,1))
%!error<lognfit: X and FREQ vector mismatch.> ...
%! lognfit (ones (20,1), [], zeros(20,1), ones(25,1))
%!error<lognfit: > lognfit (ones (20,1), [], zeros(20,1), ones(20,1), "options")
