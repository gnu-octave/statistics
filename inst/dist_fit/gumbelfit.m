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
## @deftypefn  {statistics} {@var{paramhat} =} gumbelfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gumbelfit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} gumbelfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} gumbelfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} gumbelfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@dots{}] =} gumbelfit (@var{x}, @var{alpha}, @var{censor}, @var{freq}, @var{options})
##
## Estimate parameters and confidence intervals for Gumbel distribution.
##
## @code{@var{paramhat} = gumbelfit (@var{x})} returns the maximum likelihood
## estimates of the parameters of the Gumbel distribution (also known as
## the extreme value or the type I generalized extreme value distribution) given
## in @var{x}.  @qcode{@var{paramhat}(1)} is the location parameter, @var{mu},
## and @qcode{@var{paramhat}(2)} is the scale parameter, @var{beta}.
##
## @code{[@var{paramhat}, @var{paramci}] = gumbelfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimates.
##
## @code{[@dots{}] = gumbelfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals for the
## parameter estimates.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = gumbelfit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = gumbelfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @code{[@dots{}] = gumbelfit (@dots{}, @var{options})} specifies control
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
## The Gumbel distribution is used to model the distribution of the maximum (or
## the minimum) of a number of samples of various distributions.  This version
## is suitable for modeling maxima.  For modeling minima, use the alternative
## extreme value fitting function, @code{evfit}.
##
## Further information about the Gumbel distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gumbel_distribution}
##
## @seealso{gumbelcdf, gumbelinv, gumbelpdf, gumbelrnd, gumbellike, gumbelstat,
## evfit}
## @end deftypefn

function [paramhat, paramci] = gumbelfit (x, alpha, censor, freq, options)

  ## Check X for being a double precision vector
  if (! isvector (x) || ! isa (x, "double"))
    error ("gumbelfit: X must be a double-precision vector.");
  endif

  ## Check that X does not contain missing values (NaNs)
  if (any (isnan (x)))
    error ("gumbelfit: X must NOT contain missing values (NaNs).");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  else
    if (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
      error ("gumbelfit: wrong value for ALPHA.");
    endif
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("gumbelfit: X and CENSOR vectors mismatch.");
  endif

  ## Parse FREQ argument or add default
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("gumbelfit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("gumbelfit: FREQ must not contain negative values.");
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
      error (strcat ("gumbelfit: 'options' 5th argument must be a", ...
                     " structure with 'Display', 'MaxFunEvals',", ...
                     " 'MaxIter', and 'TolX' fields present."));
    endif
  endif

  ## Remove zeros and NaNs from frequency vector (if necessary)
  if (! all (freq == 1))
    remove = freq == 0 | isnan (freq);
    x(remove) = [];
    censor(remove) = [];
    freq(remove) = [];
  endif

  ## If X is a column vector, make X, CENSOR, and FREQ row vectors
  if (size (x, 1) > 1)
    x = x(:)';
    censor = censor(:)';
    freq = freq(:)';
  endif

  ## Call evfit to do the actual computation on the negative X
  try
    [paramhat, paramci] = evfit (-x, alpha, censor, freq, options);
  catch
    error ("gumbelfit: no solution for maximum likelihood estimates.");
  end_try_catch

  ## Flip sign on estimated parameter MU
  paramhat(1) = -paramhat(1);

  ## Flip sign on confidence intervals of parameter MU
  paramci(:,1) = -flip (paramci(:,1));

endfunction

%!demo
%! ## Sample 3 populations from different Gumbel distributions
%! rand ("seed", 1);    # for reproducibility
%! r1 = gumbelrnd (2, 5, 400, 1);
%! rand ("seed", 11);    # for reproducibility
%! r2 = gumbelrnd (-5, 3, 400, 1);
%! rand ("seed", 16);    # for reproducibility
%! r3 = gumbelrnd (14, 8, 400, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 25, 0.32);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! ylim ([0, 0.28])
%! xlim ([-11, 50]);
%! hold on
%!
%! ## Estimate their MU and BETA parameters
%! mu_betaA = gumbelfit (r(:,1));
%! mu_betaB = gumbelfit (r(:,2));
%! mu_betaC = gumbelfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [min(r(:)):max(r(:))];
%! y = gumbelpdf (x, mu_betaA(1), mu_betaA(2));
%! plot (x, y, "-pr");
%! y = gumbelpdf (x, mu_betaB(1), mu_betaB(2));
%! plot (x, y, "-sg");
%! y = gumbelpdf (x, mu_betaC(1), mu_betaC(2));
%! plot (x, y, "-^c");
%! legend ({"Normalized HIST of sample 1 with μ=2 and β=5", ...
%!          "Normalized HIST of sample 2 with μ=-5 and β=3", ...
%!          "Normalized HIST of sample 3 with μ=14 and β=8", ...
%!          sprintf("PDF for sample 1 with estimated μ=%0.2f and β=%0.2f", ...
%!                  mu_betaA(1), mu_betaA(2)), ...
%!          sprintf("PDF for sample 2 with estimated μ=%0.2f and β=%0.2f", ...
%!                  mu_betaB(1), mu_betaB(2)), ...
%!          sprintf("PDF for sample 3 with estimated μ=%0.2f and β=%0.2f", ...
%!                  mu_betaC(1), mu_betaC(2))})
%! title ("Three population samples from different Gumbel distributions")
%! hold off

## Test output
%!test
%! x = 1:50;
%! [paramhat, paramci] = gumbelfit (x);
%! paramhat_out = [18.3188, 13.0509];
%! paramci_out = [14.4882, 10.5294; 22.1495, 16.1763];
%! assert (paramhat, paramhat_out, 1e-4);
%! assert (paramci, paramci_out, 1e-4);
%!test
%! x = 1:50;
%! [paramhat, paramci] = gumbelfit (x, 0.01);
%! paramci_out = [13.2845, 9.8426; 23.3532, 17.3051];
%! assert (paramci, paramci_out, 1e-4);

## Test input validation
%!error<gumbelfit: X must be a double-precision vector.> gumbelfit (ones (2,5));
%!error<gumbelfit: X must be a double-precision vector.> ...
%! gumbelfit (single (ones (1,5)));
%!error<gumbelfit: X must NOT contain missing values> ...
%! gumbelfit ([1, 2, 3, 4, NaN]);
%!error<gumbelfit: wrong value for ALPHA.> gumbelfit ([1, 2, 3, 4, 5], 1.2);
%!error<gumbelfit: X and CENSOR vectors mismatch.> ...
%! gumbelfit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<gumbelfit: X and FREQ vectors mismatch.> ...
%! gumbelfit ([1, 2, 3, 4, 5], 0.05, [], [1 1 0]);
%!error<gamfit: FREQ must not contain negative values.>
%! gamfit ([1, 2, 3], 0.05, [], [1 5 -1])
%!error<gumbelfit: 'options' 5th argument> ...
%! gumbelfit ([1, 2, 3, 4, 5], 0.05, [], [], 2);
