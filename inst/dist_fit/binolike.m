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
## @deftypefn  {statistics} {@var{nlogL} =} binolike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} binolike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} binolike (@var{params}, @var{x}, @var{freq})
##
## Negative log-likelihood for the binomial distribution.
##
## @code{@var{nlogL} = binolike (@var{params}, @var{x})} returns the negative
## log likelihood of the binomial distribution with (1) parameter @var{n} and
## (2) parameter @var{ps}, given in the two-element vector @var{params}, where
## @var{n} is the number of trials and @var{ps} is the probability of success,
## given the number of successes in @var{x}.  Unlike @code{binofit}, which
## handles each element in @var{x} independently, @code{binolike} returns the
## negative log likelihood of the entire vector @var{x}.
##
## @code{[@var{nlogL}, @var{acov}] = binolike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{params} are their asymptotic variances.
##
## @code{[@dots{}] = binolike (@var{params}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## typically contains integer frequencies for the corresponding elements in
## @var{x}, but it can contain any non-integer non-negative values.  By default,
## or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Binomial_distribution}
##
## @seealso{binocdf, binoinv, binopdf, binornd, binofit, binostat}
## @end deftypefn

function [nlogL, acov] = binolike (params, x, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("binolike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("binolike: X must be a vector.");
  endif

  if (length (params) != 2)
    error ("binolike: PARAMS must be a two-element vector.");
  endif

  if (params(1) < 0 || params(1) != round (params(1)) || isinf (params(1)))
    error (strcat (["binolike: number of trials, PARAMS(1), must be a"], ...
                   [" finite non-negative integer."]));
  endif

  if (params(2) < 0 || params(2) > 1)
    error (strcat (["binolike: probability of success, PARAMS(2), must be"], ...
                   [" in the range [0,1]."]));
  endif

  ## Parse FREQ argument or add default
  if (nargin < 3 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("binolike: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("binolike: FREQ must not contain negative values.");
  endif

  ## Expand frequency vector (if necessary)
  if (! all (freq == 1))
    ## Remove NaNs and zeros
    remove = isnan (freq) | freq == 0;
    x(remove) = [];
    freq(remove) = [];
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  if (any (x < 0))
    error ("binolike: X cannot have negative values.");
  endif

  if (any (x > params(1)))
    error (strcat (["binolike: number of successes, X, must be at least"], ...
                   [" as large as the number of trials, N."]));
  endif

  ## Compute negative log-likelihood and asymptotic covariance
  n = params(1);
  ps = params(2);
  numx = length (x);
  nlogL = -sum (log (binopdf (x, n, ps)));
  tmp = ps * (1 - ps) / (n * numx);
  acov = [0, 0; 0, tmp];

endfunction

## Test output
%!assert (binolike ([3, 0.333], [0:3]), 6.8302, 1e-4)
%!assert (binolike ([3, 0.333], 0), 1.2149, 1e-4)
%!assert (binolike ([3, 0.333], 1), 0.8109, 1e-4)
%!assert (binolike ([3, 0.333], 2), 1.5056, 1e-4)
%!assert (binolike ([3, 0.333], 3), 3.2988, 1e-4)
%!test
%! [nlogL, acov] = binolike ([3, 0.333], 3);
%! assert (acov(4), 0.0740, 1e-4)

## Test input validation
%!error<binolike: function called with too few input arguments.> binolike (3.25)
%!error<binolike: X must be a vector.> binolike ([5, 0.2], ones (2))
%!error<binolike: PARAMS must be a two-element vector.> ...
%! binolike ([1, 0.2, 3], [1, 3, 5, 7])
%!error<binolike: number of trials,> binolike ([1.5, 0.2], 1)
%!error<binolike: number of trials,> binolike ([-1, 0.2], 1)
%!error<binolike: number of trials,> binolike ([Inf, 0.2], 1)
%!error<binolike: probability of success,> binolike ([5, 1.2], [3, 5])
%!error<binolike: probability of success,> binolike ([5, -0.2], [3, 5])
%!error<binolike: X and FREQ vectors mismatch.> ...
%! binolike ([5, 0.5], ones (10, 1), ones (8,1))
%!error<binolike: FREQ must not contain negative values.> ...
%! binolike ([5, 0.5], ones (1, 8), [1 1 1 1 1 1 1 -1])
%!error<binolike: X cannot have negative values.> binolike ([5, 0.2], [-1, 3])
%!error<binolike: number of successes,> binolike ([5, 0.2], [3, 5, 7])
