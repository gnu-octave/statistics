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
## @deftypefn  {statistics} {@var{pshat} =} binofit (@var{x}, @var{n})
## @deftypefnx {statistics} {[@var{pshat}, @var{psci}] =} binofit (@var{x}, @var{n})
## @deftypefnx {statistics} {[@var{pshat}, @var{psci}] =} binofit (@var{x}, @var{n}, @var{alpha})
##
## Estimate parameter and confidence intervals for the binomial distribution.
##
## @code{@var{pshat} = binofit (@var{x}, @var{n})} returns the maximum
## likelihood estimate (MLE) of the probability of success for the binomial
## distribution.  @var{x} and @var{n} are scalars containing the number of
## successes and the number of trials, respectively.  If @var{x} and @var{n} are
## vectors, @code{binofit} returns a vector of estimates whose @math{i}-th
## element is the parameter estimate for @var{x}(i) and @var{n}(i).  A scalar
## value for @var{x} or @var{n} is expanded to the same size as the other input.
##
## @code{[@var{pshat}, @var{psci}] = binofit (@var{x}, @var{n}, @var{alpha})}
## also returns the @qcode{100 * (1 - @var{alpha})} percent confidence intervals
## of the estimated parameter.  By default, the optional argument @var{alpha}
## is 0.05 corresponding to 95% confidence intervals.
##
## @code{binofit} treats a vector @var{x} as a collection of measurements from
## separate samples, and returns a vector of estimates.  If you want to treat
## @var{x} as a single sample and compute a single parameter estimate and
## confidence interval, use @qcode{binofit (sum (@var{x}), sum (@var{n}))} when
## @var{n} is a vector, and
## @qcode{binofit (sum (@var{x}), @var{n} * length (@var{x}))} when @var{n} is a
## scalar.
##
## Further information about the binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Binomial_distribution}
##
## @seealso{binocdf, binoinv, binopdf, binornd, binolike, binostat}
## @end deftypefn

function [pshat, psci] = binofit (x, n, alpha)

  ## Check input arguments
  if (nargin < 2)
    error ("binofit: too few input arguments.");
  endif
  if (any (x < 0))
    error ("binofit: X cannot have negative values.");
  endif
  if (! isvector (x))
    error ("binofit: X must be a vector.");
  endif
  if (any (n < 0) || any (n != round (n)) || any (isinf (n)))
    error ("binofit: N must be a non-negative integer.");
  endif
  if (any (x > n))
    error ("binofit: N cannot be greater than X.");
  endif

  if (nargin < 3 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("binofit: Wrong value for ALPHA.");
  endif

  ## Compute pshat
  pshat = x ./ n;

  ## Compute lower confidence interval
  nu1 = 2 * x;
  nu2 = 2 * (n - x + 1);
  F   = finv (alpha / 2, nu1, nu2);
  lb  = (nu1 .* F) ./ (nu2 + nu1 .* F);
  x0  = find (x == 0);
  if (! isempty (x0))
    lb(x0) = 0;
  endif

  ## Compute upper confidence interval
  nu1 = 2 * (x + 1);
  nu2 = 2 * (n - x);
  F   = finv (1 - alpha / 2, nu1, nu2);
  ub  = (nu1 .* F) ./ (nu2 + nu1 .* F);
  xn  = find (x == n);
  if (! isempty (xn))
    ub(xn) = 1;
  endif

  psci = [lb(:), ub(:)];

endfunction

%!demo
%! ## Sample 2 populations from 2 different binomial distibutions
%! rand ("seed", 1);    # for reproducibility
%! r1 = binornd (50, 0.15, 1000, 1);
%! rand ("seed", 2);    # for reproducibility
%! r2 = binornd (100, 0.5, 1000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 23, 0.35);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! hold on
%!
%! ## Estimate their probability of success
%! pshatA = binofit (r(:,1), 50);
%! pshatB = binofit (r(:,2), 100);
%!
%! ## Plot their estimated PDFs
%! x = [min(r(:,1)):max(r(:,1))];
%! y = binopdf (x, 50, mean (pshatA));
%! plot (x, y, "-pg");
%! x = [min(r(:,2)):max(r(:,2))];
%! y = binopdf (x, 100, mean (pshatB));
%! plot (x, y, "-sc");
%! ylim ([0, 0.2])
%! legend ({"Normalized HIST of sample 1 with ps=0.15", ...
%!          "Normalized HIST of sample 2 with ps=0.50", ...
%!          sprintf("PDF for sample 1 with estimated ps=%0.2f", ...
%!                  mean (pshatA)), ...
%!          sprintf("PDF for sample 2 with estimated ps=%0.2f", ...
%!                  mean (pshatB))})
%! title ("Two population samples from different binomial distibutions")
%! hold off

## Test output
%!test
%! x = 0:3;
%! [pshat, psci] = binofit (x, 3);
%! assert (pshat, [0, 0.3333, 0.6667, 1], 1e-4);
%! assert (psci(1,:), [0, 0.7076], 1e-4);
%! assert (psci(2,:), [0.0084, 0.9057], 1e-4);
%! assert (psci(3,:), [0.0943, 0.9916], 1e-4);
%! assert (psci(4,:), [0.2924, 1.0000], 1e-4);

## Test input validation
%!error<binofit: too few input arguments.> binofit ([1 2 3 4])
%!error<binofit: X cannot have negative values.> binofit (-1, [1 2 3 3])
%!error<binofit: N must be a non-negative integer.> binofit (1, [1 2 -1 3])
%!error<binofit: Wrong value for ALPHA.> binofit (1, [1 2 3], 0)
%!error<binofit: Wrong value for ALPHA.> binofit (1, [1 2 3], 1.2)
%!error<binofit: Wrong value for ALPHA.> binofit (1, [1 2 3], [0.02 0.05])
