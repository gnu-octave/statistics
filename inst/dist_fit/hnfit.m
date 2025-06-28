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
## @deftypefn  {statistics} {[@var{paramhat}, @var{paramci}] =} hnfit (@var{x}, @var{mu})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} hnfit (@var{x}, @var{mu}, @var{alpha})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} hnfit (@var{x}, @var{mu}, @var{alpha}, @var{freq})
##
## Estimate parameters and confidence intervals for the half-normal distribution.
##
## @code{@var{paramhat} = hnfit (@var{x}, @var{mu})} returns the maximum
## likelihood estimates of the parameters of the half-normal distribution given
## the data in vector @var{x} and the location parameter @var{mu}.
## @qcode{@var{paramhat}(1)} is the location parameter, @var{mu}, and
## @qcode{@var{paramhat}(2)} is the scale parameter, @var{sigma}.  Although
## @var{mu} is returned in the estimated @var{paramhat}, @code{hnfit} does not
## estimate the location parameter @var{mu}, and it must be assumed to be known,
## given as a fixed parameter in input argument @var{mu}.
##
## @code{[@var{paramhat}, @var{paramci}] = hnfit (@var{x}, @var{mu})} returns
## the 95% confidence intervals for the estimated scale parameter @var{sigma}.
## The first colummn of @var{paramci} includes the location parameter @var{mu}
## without any confidence bounds.
##
## @code{[@dots{}] = hnfit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals of the estimated
## scale parameter.  By default, the optional argument @var{alpha} is 0.05
## corresponding to 95% confidence intervals.
##
## @code{[@dots{}] = hnfit (@var{params}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## must contain non-negative integer frequencies for the corresponding elements
## in @var{x}.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## The half-normal CDF is only defined for @qcode{@var{x} >= @var{mu}}.
##
## Further information about the half-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Half-normal_distribution}
##
## @seealso{hncdf, hninv, hnpdf, hnrnd, hnlike, hnstat}
## @end deftypefn

function [paramhat, paramci] = hnfit (x, mu, alpha, freq)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("hnfit: function called with too few input arguments.");
  endif

  ## Check X for being a vector
  if (isempty (x))
    phat = nan (1, 2, class (x));
    pci = nan (2, 2, class (x));
    return
  elseif (! isvector (x) || ! isreal (x))
    error ("hnfit: X must be a vector of real values.");
  endif

  ## Check for MU being a scalar real value
  if (! isscalar (mu) || ! isreal (mu))
    error ("hnfit: MU must be a real scalar value.");
  endif

  ## Check X >= MU
  if (any (x < mu))
    error ("hnfit: X cannot contain values less than MU.");
  endif

  ## Parse ALPHA argument or add default
  if (nargin < 3 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("hnfit: wrong value for ALPHA.");
  endif

  ## Parse FREQ argument or add default
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("hnfit: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("hnfit: FREQ must not contain negative values.");
  endif

  ## Expand frequency vector (if necessary)
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  ## Estimate parameters
  sz = numel (x);
  x = x - mu;
  sigmahat = sqrt (sum (x .* x) ./ sz);
  paramhat = [mu, sigmahat];

  ## Compute confidence intervals
  if (nargout == 2)
    chi2cr = chi2inv ([alpha/2, 1-alpha/2], sz);
    shatlo = sigmahat * sqrt (sz / chi2inv (1 - alpha / 2, sz));
    shathi = sigmahat * sqrt (sz / chi2inv (alpha / 2, sz));
    paramci = [mu, shatlo; mu, shathi];
  endif

endfunction

%!demo
%! ## Sample 2 populations from different half-normal distributions
%! rand ("seed", 1);   # for reproducibility
%! r1 = hnrnd (0, 5, 5000, 1);
%! rand ("seed", 2);   # for reproducibility
%! r2 = hnrnd (0, 2, 5000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, [0.5:20], 1);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! hold on
%!
%! ## Estimate their shape parameters
%! mu_sigmaA = hnfit (r(:,1), 0);
%! mu_sigmaB = hnfit (r(:,2), 0);
%!
%! ## Plot their estimated PDFs
%! x = [0:0.2:10];
%! y = hnpdf (x, mu_sigmaA(1), mu_sigmaA(2));
%! plot (x, y, "-pr");
%! y = hnpdf (x, mu_sigmaB(1), mu_sigmaB(2));
%! plot (x, y, "-sg");
%! xlim ([0, 10])
%! ylim ([0, 0.5])
%! legend ({"Normalized HIST of sample 1 with μ=0 and σ=5", ...
%!          "Normalized HIST of sample 2 with μ=0 and σ=2", ...
%!          sprintf("PDF for sample 1 with estimated μ=%0.2f and σ=%0.2f", ...
%!                  mu_sigmaA(1), mu_sigmaA(2)), ...
%!          sprintf("PDF for sample 2 with estimated μ=%0.2f and σ=%0.2f", ...
%!                  mu_sigmaB(1), mu_sigmaB(2))})
%! title ("Two population samples from different half-normal distributions")
%! hold off

## Test output
%!test
%! x = 1:20;
%! [paramhat, paramci] = hnfit (x, 0);
%! assert (paramhat, [0, 11.9791], 1e-4);
%! assert (paramci, [0, 9.1648; 0, 17.2987], 1e-4);
%!test
%! x = 1:20;
%! [paramhat, paramci] = hnfit (x, 0, 0.01);
%! assert (paramci, [0, 8.4709; 0, 19.6487], 1e-4);

## Test input validation
%!error<hnfit: function called with too few input arguments.> hnfit ()
%!error<hnfit: function called with too few input arguments.> hnfit (1)
%!error<hnfit: X must be a vector of real values.> hnfit ([0.2, 0.5+i], 0);
%!error<hnfit: X must be a vector of real values.> hnfit (ones (2,2) * 0.5, 0);
%!error<hnfit: MU must be a real scalar value.> ...
%! hnfit ([0.5, 1.2], [0, 1]);
%!error<hnfit: MU must be a real scalar value.> ...
%! hnfit ([0.5, 1.2], 5+i);
%!error<hnfit: X cannot contain values less than MU.> ...
%! hnfit ([1:5], 2);
%!error<hnfit: wrong value for ALPHA.> hnfit ([0.01:0.1:0.99], 0, 1.2);
%!error<hnfit: wrong value for ALPHA.> hnfit ([0.01:0.1:0.99], 0, i);
%!error<hnfit: wrong value for ALPHA.> hnfit ([0.01:0.1:0.99], 0, -1);
%!error<hnfit: wrong value for ALPHA.> hnfit ([0.01:0.1:0.99], 0, [0.05, 0.01]);
%!error<hnfit: X and FREQ vectors mismatch.>
%! hnfit ([1 2 3], 0, [], [1 5])
%!error<hnfit: FREQ must not contain negative values.>
%! hnfit ([1 2 3], 0, [], [1 5 -1])
