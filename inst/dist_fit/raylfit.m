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
## @deftypefn  {statistics} {@var{sigmaA} =} raylfit (@var{x})
## @deftypefnx {statistics} {[@var{sigmaA}, @var{sigmaci}] =} raylfit (@var{x})
## @deftypefnx {statistics} {[@var{sigmaA}, @var{sigmaci}] =} raylfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{sigmaA}, @var{sigmaci}] =} raylfit (@var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@var{sigmaA}, @var{sigmaci}] =} raylfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})
##
## Estimate parameter and confidence intervals for the Rayleigh distribution.
##
## @code{@var{sigmaA} = raylfit (@var{x})} returns the maximum likelihood
## estimate of the rate parameter, @var{lambda}, of the Rayleigh distribution
## given the data in @var{x}.  @var{x} must be a vector of non-negative values.
##
## @code{[@var{sigmaA}, @var{sigmaci}] = raylfit (@var{x})} returns the 95%
## confidence intervals for the parameter estimate.
##
## @code{[@var{sigmaA}, @var{sigmaci}] = raylfit (@var{x}, @var{alpha})}
## also returns the @qcode{100 * (1 - @var{alpha})} percent confidence intervals
## of the estimated parameter.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = raylfit (@var{x}, @var{alpha}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = raylfit (@var{x}, @var{alpha}, @var{censor}, @var{freq})}
## accepts a frequency vector or matrix, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}.  @var{freq} cannot contain negative values.
##
## Further information about the Rayleigh distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rayleigh_distribution}
##
## @seealso{raylcdf, raylinv, raylpdf, raylrnd, rayllike, raylstat}
## @end deftypefn

function [sigmaA, sigmaci] = raylfit (x, alpha, censor, freq)

  ## Check input arguments
  if (any (x < 0))
    error ("raylfit: X cannot have negative values.");
  endif
  if (! isvector (x))
    error ("raylfit: X must be a vector.");
  endif

  ## Check alpha
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("raylfit: wrong value for ALPHA.");
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("raylfit: X and CENSOR vectors mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("raylfit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("raylfit: FREQ must not contain negative values.");
  endif

  ## Remove any censored data
  censored = censor == 1;
  freq(censored) = [];
  x(censored) = [];

  ## Expand frequency vector (if necessary)
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  ## Compute sigmaA
  sigmaA = sqrt (0.5 * mean (x .^ 2));

  ## Compute confidence intervals (based on chi-squared)
  if (nargout > 1)
    sx = 2 * numel (x);
    ci = [1-alpha/2; alpha/2];
    sigmaci = sqrt (sx * sigmaA .^ 2 ./ chi2inv (ci, sx));
  endif

endfunction

%!demo
%! ## Sample 3 populations from 3 different Rayleigh distibutions
%! rand ("seed", 2);    # for reproducibility
%! r1 = raylrnd (1, 1000, 1);
%! rand ("seed", 2);    # for reproducibility
%! r2 = raylrnd (2, 1000, 1);
%! rand ("seed", 3);    # for reproducibility
%! r3 = raylrnd (4, 1000, 1);
%! r = [r1, r2, r3];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0.5:0.5:10.5], 2);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! set (h(3), "facecolor", "r");
%! hold on
%!
%! ## Estimate their lambda parameter
%! sigmaA = raylfit (r(:,1));
%! sigmaB = raylfit (r(:,2));
%! sigmaC = raylfit (r(:,3));
%!
%! ## Plot their estimated PDFs
%! x = [0:0.1:10];
%! y = raylpdf (x, sigmaA);
%! plot (x, y, "-pr");
%! y = raylpdf (x, sigmaB);
%! plot (x, y, "-sg");
%! y = raylpdf (x, sigmaC);
%! plot (x, y, "-^c");
%! xlim ([0, 10])
%! ylim ([0, 0.7])
%! legend ({"Normalized HIST of sample 1 with σ=1", ...
%!          "Normalized HIST of sample 2 with σ=2", ...
%!          "Normalized HIST of sample 3 with σ=4", ...
%!          sprintf("PDF for sample 1 with estimated σ=%0.2f", ...
%!                  sigmaA), ...
%!          sprintf("PDF for sample 2 with estimated σ=%0.2f", ...
%!                  sigmaB), ...
%!          sprintf("PDF for sample 3 with estimated σ=%0.2f", ...
%!                  sigmaC)})
%! title ("Three population samples from different Rayleigh distibutions")
%! hold off

## Test output
%!test
%! x = [1 3 2 4 5 4 3 4];
%! [shat, sci] = raylfit (x);
%! assert (shat, 2.4495, 1e-4)
%! assert (sci, [1.8243; 3.7279], 1e-4)
%!test
%! x = [1 3 2 4 5 4 3 4];
%! [shat, sci] = raylfit (x, 0.01);
%! assert (shat, 2.4495, 1e-4)
%! assert (sci, [1.6738; 4.3208], 1e-4)
%!test
%! x = [1 2 3 4 5];
%! f = [1 1 2 3 1];
%! [shat, sci] = raylfit (x, [], [], f);
%! assert (shat, 2.4495, 1e-4)
%! assert (sci, [1.8243; 3.7279], 1e-4)
%!test
%! x = [1 2 3 4 5];
%! f = [1 1 2 3 1];
%! [shat, sci] = raylfit (x, 0.01, [], f);
%! assert (shat, 2.4495, 1e-4)
%! assert (sci, [1.6738; 4.3208], 1e-4)
%!test
%! x = [1 2 3 4 5 6];
%! c = [0 0 0 0 0 1];
%! f = [1 1 2 3 1 1];
%! [shat, sci] = raylfit (x, 0.01, c, f);
%! assert (shat, 2.4495, 1e-4)
%! assert (sci, [1.6738; 4.3208], 1e-4)

## Test input validation
%!error<raylfit: X must be a vector.> raylfit (ones (2,5));
%!error<raylfit: X cannot have negative values.> raylfit ([1 2 -1 3])
%!error<raylfit: wrong value for ALPHA.> raylfit ([1 2 3], 0)
%!error<raylfit: wrong value for ALPHA.> raylfit ([1 2 3], 1.2)
%!error<raylfit: wrong value for ALPHA.> raylfit ([1 2 3], [0.02 0.05])
%!error<raylfit: X and CENSOR vectors mismatch.> ...
%! raylfit ([1, 2, 3, 4, 5], 0.05, [1 1 0]);
%!error<raylfit: X and CENSOR vectors mismatch.> ...
%! raylfit ([1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<raylfit: X and FREQ vectors mismatch.> ...
%! raylfit ([1, 2, 3, 4, 5], 0.05, zeros (1,5), [1 1 0]);
%!error<raylfit: X and FREQ vectors mismatch.> ...
%! raylfit ([1, 2, 3, 4, 5], [], [], [1 1 0 1 1]');
%!error<raylfit: X and FREQ vectors mismatch.>
%! raylfit ([1 2 3], [], [], [1 5])
%!error<raylfit: FREQ must not contain negative values.>
%! raylfit ([1 2 3], [], [], [1 5 -1])
