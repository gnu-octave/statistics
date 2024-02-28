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
## @deftypefn  {statistics} {@var{paramhat} =} betafit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} betafit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} betafit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} betafit (@var{x}, @var{alpha}, @var{freq})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} betafit (@var{x}, @var{alpha}, @var{options})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} betafit (@var{x}, @var{alpha}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for the Beta distribution.
##
## @code{@var{paramhat} = betafit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the Beta distribution given the data in vector
## @var{x}.  @qcode{@var{paramhat}([1, 2])} corresponds to the @math{α} and
## @math{β} shape parameters, respectively.  Missing values, @qcode{NaNs}, are
## ignored.
##
## @code{[@var{paramhat}, @var{paramci}] = betafit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = betafit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals of the estimated
## parameter.  By default, the optional argument @var{alpha} is 0.05
## corresponding to 95% confidence intervals.
##
## @code{[@dots{}] = betafit (@var{params}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## must contain non-negative integer frequencies for the corresponding elements
## in @var{x}.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
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
## The Beta distribution is defined on the open interval @math{(0,1)}.  However,
## @code{betafit} can also compute the unbounded beta likelihood function for
## data that include exact zeros or ones.  In such cases, zeros and ones are
## treated as if they were values that have been left-censored at
## @qcode{sqrt (realmin)} or right-censored at @qcode{1 - eps/2}, respectively.
##
## Further information about the Beta distribution can be found at
## @url{https://en.wikipedia.org/wiki/Beta_distribution}
##
## @seealso{betacdf, betainv, betapdf, betarnd, betalike, betastat}
## @end deftypefn

function [paramhat, paramci] = betafit (x, alpha, varargin)

  ## Check X for being a vector
  if (isempty (x))
    phat = nan (1, 2, class (x));
    pci = nan (2, 2, class (x));
    return
  elseif (! isvector (x) || ! isreal (x))
    error ("betafit: X must be a vector of real values.");
  endif

  ## Check that X contains values in the range [0,1]
  if (any (x < 0) || any (x > 1))
    error ("betafit: X must be in the range [0,1].");
  endif

  ## Check X being a constant vector
  if (min (x) == max(x))
    error ("betafit: X must contain distinct values.");
  endif

  ## Check ALPHA
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("betafit: wrong value for ALPHA.");
    endif
  endif

  ## Add defaults
  freq = ones (size (x));
  options.Display = "off";
  options.MaxFunEvals = 400;
  options.MaxIter = 200;
  options.TolX = 1e-6;

  ## Check extra arguments for FREQ vector and/or 'options' structure
  if (nargin > 2)
    if (numel (varargin) == 1 && isstruct (varargin{1}))
      options = varargin{1};
    elseif (numel (varargin) == 1 && isnumeric (varargin{1}))
      freq = varargin{1};
    elseif (numel (varargin) == 2)
      freq = varargin{1};
      options = varargin{2};
    endif
    if (isempty (freq))
      freq = ones (size (x));
    endif
    ## Check for valid freq vector
    if (! isequal (size (x), size (freq)))
      error ("betafit: X and FREQ vectors mismatch.");
    elseif (any (freq < 0))
      error ("betafit: FREQ must not contain negative values.");
    elseif (any (fix (freq) != freq))
      error ("betafit: FREQ must contain integer values.");
    endif
    ## Check for valid options structure
    if (! isstruct (options) || ! isfield (options, "Display") ||
        ! isfield (options, "MaxFunEvals") || ! isfield (options, "MaxIter")
                                           || ! isfield (options, "TolX"))
      error (strcat (["betafit: 'options' argument must be a"], ...
                     [" structure with 'Display', 'MaxFunEvals',"], ...
                     [" 'MaxIter', and 'TolX' fields present."]));
    endif
  endif

  ## Remove missing values
  remove = isnan (x) | isnan (freq);
  x(remove) = [];
  freq(remove) = [];

  ## Expand frequency
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  ## Estimate initial parameters
  numx = length (x);
  tmp1 = prod ((1 - x) .^ (1 / numx));
  tmp2 = prod (x .^ (1 / numx));
  tmp3 = (1 - tmp1 - tmp2);
  ahat = 0.5 * (1 - tmp1) / tmp3;
  bhat = 0.5 * (1 - tmp2) / tmp3;
  init = log ([ahat, bhat]);

  ## Add tolerance for boundary conditions
  x_lo = sqrt (realmin (class (x)));
  x_hi = 1 - eps (class (x)) / 2;

  ## All values are strictly within the interval (0,1)
  if (all (x > x_lo) && all (x < x_hi))
    sumlogx = sum (log (x));
    sumlog1px = sum (log1p (-x));
    paramhat = fminsearch (@cont_negloglike, init, options);
    paramhat = exp (paramhat);

  ## Find boundary elements and process them separately
  else
    num0 = sum (x < x_lo);
    num1 = sum (x > x_hi);
    x_ct = x (x > x_lo & x < x_hi);
    numx = length (x_ct);
    sumlogx = sum (log (x_ct));
    sumlog1px = sum (log1p (-x_ct));
    paramhat = fminsearch (@mixed_negloglike, init, options);
    paramhat = exp (paramhat);
  endif

  ## Compute confidence intervals
  if (nargout == 2)
    [~, acov] = betalike (paramhat,x);
    logphat = log (paramhat);
    serrlog = sqrt (diag (acov))' ./ paramhat;
    p_int = [alpha/2; 1-alpha/2];
    paramci = exp (norminv ([p_int p_int], ...
                            [logphat; logphat], [serrlog; serrlog]));
  endif

  ## Continuous Negative log-likelihood function
  function nll = cont_negloglike (params)
    params = exp (params);
    nll = numx * betaln (params(1), params(2)) - (params(1) - 1) ...
               * sumlogx - (params(2) - 1) * sumlog1px;
  endfunction

  ## Unbounded Negative log-likelihood function
  function nll = mixed_negloglike (params)
    params = exp (params);
    nll = numx * betaln (params(1), params(2)) - (params(1) - 1) ...
               * sumlogx - (params(2) - 1) * sumlog1px;
    ## Handle zeros
    if (num0 > 0)
      nll = nll - num0 * log (betainc (x_lo, params(1), params(2), "lower"));
    endif
    ## Handle ones
    if (num1 > 0)
      nll = nll - num1 * log (betainc (x_hi, params(1), params(2), "upper"));
    endif
  endfunction

endfunction

%!demo
%! ## Sample 2 populations from different Beta distibutions
%! randg ("seed", 1);   # for reproducibility
%! r1 = betarnd (2, 5, 500, 1);
%! randg ("seed", 2);   # for reproducibility
%! r2 = betarnd (2, 2, 500, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 12, 15);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! hold on
%!
%! ## Estimate their shape parameters
%! a_b_A = betafit (r(:,1));
%! a_b_B = betafit (r(:,2));
%!
%! ## Plot their estimated PDFs
%! x = [min(r(:)):0.01:max(r(:))];
%! y = betapdf (x, a_b_A(1), a_b_A(2));
%! plot (x, y, "-pr");
%! y = betapdf (x, a_b_B(1), a_b_B(2));
%! plot (x, y, "-sg");
%! ylim ([0, 4])
%! legend ({"Normalized HIST of sample 1 with α=2 and β=5", ...
%!          "Normalized HIST of sample 2 with α=2 and β=2", ...
%!          sprintf("PDF for sample 1 with estimated α=%0.2f and β=%0.2f", ...
%!                  a_b_A(1), a_b_A(2)), ...
%!          sprintf("PDF for sample 2 with estimated α=%0.2f and β=%0.2f", ...
%!                  a_b_B(1), a_b_B(2))})
%! title ("Two population samples from different Beta distibutions")
%! hold off

## Test output
%!test
%! x = 0.01:0.02:0.99;
%! [paramhat, paramci] = betafit (x);
%! paramhat_out = [1.0199, 1.0199];
%! paramci_out = [0.6947, 0.6947; 1.4974, 1.4974];
%! assert (paramhat, paramhat_out, 1e-4);
%! assert (paramci, paramci_out, 1e-4);
%!test
%! x = 0.01:0.02:0.99;
%! [paramhat, paramci] = betafit (x, 0.01);
%! paramci_out = [0.6157, 0.6157; 1.6895, 1.6895];
%! assert (paramci, paramci_out, 1e-4);
%!test
%! x = 0.00:0.02:1;
%! [paramhat, paramci] = betafit (x);
%! paramhat_out = [0.0875, 0.1913];
%! paramci_out = [0.0822, 0.1490; 0.0931, 0.2455];
%! assert (paramhat, paramhat_out, 1e-4);
%! assert (paramci, paramci_out, 1e-4);

## Test input validation
%!error<betafit: X must be a vector of real values.> betafit ([0.2, 0.5+i]);
%!error<betafit: X must be a vector of real values.> betafit (ones (2,2) * 0.5);
%!error<betafit: X must be in the range> betafit ([0.5, 1.2]);
%!error<betafit: X must contain distinct values.> betafit ([0.1, 0.1]);
%!error<betafit: wrong value for ALPHA.> betafit ([0.01:0.1:0.99], 1.2);
%!error<betafit: X and FREQ vectors mismatch.> ...
%! betafit ([0.01:0.01:0.05], 0.05, [1, 2, 3, 2]);
%!error<betafit: FREQ must not contain negative values.> ...
%! betafit ([0.01:0.01:0.05], 0.05, [1, 2, 3, 2, -1]);
%!error<betafit: FREQ must contain integer values.> ...
%! betafit ([0.01:0.01:0.05], 0.05, [1, 2, 3, 2, 1.5]);
%!error<betafit: 'options' argument must be a structure> ...
%! betafit ([0.01:0.01:0.05], 0.05, struct ("option", 234));
%!error<betafit: 'options' argument must be a structure> ...
%! betafit ([0.01:0.01:0.05], 0.05, ones (1,5), struct ("option", 234));
