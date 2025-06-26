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
## @deftypefn  {statistics} {@var{lambdahat} =} poissfit (@var{x})
## @deftypefnx {statistics} {[@var{lambdahat}, @var{lambdaci}] =} poissfit (@var{x})
## @deftypefnx {statistics} {[@var{lambdahat}, @var{lambdaci}] =} poissfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{lambdahat}, @var{lambdaci}] =} poissfit (@var{x}, @var{alpha}, @var{freq})
##
## Estimate parameter and confidence intervals for the Poisson distribution.
##
## @code{@var{lambdahat} = poissfit (@var{x})} returns the maximum likelihood
## estimate of the rate parameter, @var{lambda}, of the Poisson distribution
## given the data in @var{x}.  @var{x} must be a vector of non-negative values.
##
## @code{[@var{lambdahat}, @var{lambdaci}] = poissfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimate.
##
## @code{[@var{lambdahat}, @var{lambdaci}] = poissfit (@var{x}, @var{alpha})}
## also returns the @qcode{100 * (1 - @var{alpha})} percent confidence intervals
## of the estimated parameter.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = poissfit (@var{x}, @var{alpha}, @var{freq})} accepts a
## frequency vector or matrix, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}.  @var{freq} cannot contain negative values.
##
## Further information about the Poisson distribution can be found at
## @url{https://en.wikipedia.org/wiki/Poisson_distribution}
##
## @seealso{poisscdf, poissinv, poisspdf, poissrnd, poisslike, poisstat}
## @end deftypefn

function [lambdahat, lambdaci] = poissfit (x, alpha, freq)

  ## Check input arguments
  if (any (x < 0))
    error ("poissfit: X cannot have negative values.");
  endif

  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("poissfit: wrong value for ALPHA.");
  endif

  if (nargin < 3 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("poissfit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("poissfit: FREQ must not contain negative values.");
  endif

  if (isvector (x))
    x = x(:);
    freq = freq(:);
  endif

  ## Compute lambdahat
  n = sum (freq, 1);
  lambdahat = double (sum (x .* freq) ./ n);

  ## Compute confidence intervals
  lambdasum = n .* lambdahat;
  ## Select elements for exact method or normal approximation
  k = (lambdasum < 100);
  if (any (k))    # exact method
    lb(k) = chi2inv (alpha / 2, 2 * lambdasum(k)) / 2;
    ub(k) = chi2inv (1 - alpha / 2, 2 * (lambdasum(k) + 1)) / 2;
  endif
  k = ! k;
  if (any (k))    # normal approximation
      lb(k) = norminv (alpha / 2, lambdasum(k), sqrt (lambdasum(k)));
      ub(k) = norminv (1 - alpha / 2, lambdasum(k), sqrt (lambdasum(k)));
  endif

  lambdaci = [lb; ub] / n;

endfunction

%!demo
%! ## Sample 3 populations from 3 different Poisson distributions
%! randp ("seed", 2);    # for reproducibility
%! r1 = poissrnd (1, 1000, 1);
%! randp ("seed", 2);    # for reproducibility
%! r2 = poissrnd (4, 1000, 1);
%! randp ("seed", 3);    # for reproducibility
%! r3 = poissrnd (10, 1000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0:20], 1);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! hold on
%!
%! ## Estimate their lambda parameter
%! lambdahat = poissfit (r);
%!
%! ## Plot their estimated PDFs
%! x = [0:20];
%! y = poisspdf (x, lambdahat(1));
%! plot (x, y, "-pr");
%! y = poisspdf (x, lambdahat(2));
%! plot (x, y, "-sg");
%! y = poisspdf (x, lambdahat(3));
%! plot (x, y, "-^c");
%! xlim ([0, 20])
%! ylim ([0, 0.4])
%! legend ({"Normalized HIST of sample 1 with λ=1", ...
%!          "Normalized HIST of sample 2 with λ=4", ...
%!          "Normalized HIST of sample 3 with λ=10", ...
%!          sprintf("PDF for sample 1 with estimated λ=%0.2f", ...
%!                  lambdahat(1)), ...
%!          sprintf("PDF for sample 2 with estimated λ=%0.2f", ...
%!                  lambdahat(2)), ...
%!          sprintf("PDF for sample 3 with estimated λ=%0.2f", ...
%!                  lambdahat(3))})
%! title ("Three population samples from different Poisson distributions")
%! hold off

## Test output
%!test
%! x = [1 3 2 4 5 4 3 4];
%! [lhat, lci] = poissfit (x);
%! assert (lhat, 3.25)
%! assert (lci, [2.123007901949543; 4.762003010390628], 1e-14)
%!test
%! x = [1 3 2 4 5 4 3 4];
%! [lhat, lci] = poissfit (x, 0.01);
%! assert (lhat, 3.25)
%! assert (lci, [1.842572740234582; 5.281369033298528], 1e-14)
%!test
%! x = [1 2 3 4 5];
%! f = [1 1 2 3 1];
%! [lhat, lci] = poissfit (x, [], f);
%! assert (lhat, 3.25)
%! assert (lci, [2.123007901949543; 4.762003010390628], 1e-14)
%!test
%! x = [1 2 3 4 5];
%! f = [1 1 2 3 1];
%! [lhat, lci] = poissfit (x, 0.01, f);
%! assert (lhat, 3.25)
%! assert (lci, [1.842572740234582; 5.281369033298528], 1e-14)

## Test input validation
%!error<poissfit: X cannot have negative values.> poissfit ([1 2 -1 3])
%!error<poissfit: wrong value for ALPHA.> poissfit ([1 2 3], 0)
%!error<poissfit: wrong value for ALPHA.> poissfit ([1 2 3], 1.2)
%!error<poissfit: wrong value for ALPHA.> poissfit ([1 2 3], [0.02 0.05])
%!error<poissfit: X and FREQ vectors mismatch.>
%! poissfit ([1 2 3], [], [1 5])
%!error<poissfit: FREQ must not contain negative values.>
%! poissfit ([1 2 3], [], [1 5 -1])
