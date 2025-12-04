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
## @deftypefn  {statistics} {@var{nlogL} =} nbinlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{avar}] =} nbinlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{avar}] =} nbinlike (@var{params}, @var{x}, @var{freq})
##
## Negative log-likelihood for the negative binomial distribution.
##
## @code{@var{nlogL} = nbinlike (@var{params}, @var{x})} returns the negative
## log likelihood of the negative binomial distribution with (1) parameter
## @var{r} and (2) parameter @var{ps}, given in the two-element vector
## @var{params}, where @var{r} is the number of successes until the experiment
## is stopped and @var{ps} is the probability of success in each experiment,
## given the number of failures in @var{x}.
##
## @code{[@var{nlogL}, @var{avar}] = nbinlike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{avar}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{params} are their asymptotic variances.
##
## @code{[@dots{}] = nbinlike (@var{params}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## must contain non-negative integer frequencies for the corresponding elements
## in @var{x}.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## When @var{r} is an integer, the negative binomial distribution is also known
## as the Pascal distribution and it models the number of failures in @var{x}
## before a specified number of successes is reached in a series of independent,
## identical trials.  Its parameters are the probability of success in a single
## trial, @var{ps}, and the number of successes, @var{r}.  A special case of the
## negative binomial distribution, when @qcode{@var{r} = 1}, is the geometric
## distribution, which models the number of failures before the first success.
##
## @var{r} can also have non-integer positive values, in which form the negative
## binomial distribution, also known as the Polya distribution, has no
## interpretation in terms of repeated trials, but, like the Poisson
## distribution, it is useful in modeling count data.  The negative binomial
## distribution is more general than the Poisson distribution because it has a
## variance that is greater than its mean, making it suitable for count data
## that do not meet the assumptions of the Poisson distribution.  In the limit,
## as @var{r} increases to infinity, the negative binomial distribution
## approaches the Poisson distribution.
##
## Further information about the negative binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Negative_binomial_distribution}
##
## @seealso{nbincdf, nbininv, nbinpdf, nbinrnd, nbinfit, nbinstat}
## @end deftypefn

function [nlogL, avar] = nbinlike (params, x, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("nbinlike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("nbinlike: X must be a vector.");
  endif

  if (any (x < 0))
    error ("nbinlike: X cannot have negative values.");
  endif

  if (any (x != fix (x)))
    error ("nbinlike: number of failures, X, must be integers.");
  endif

  if (length (params) != 2)
    error ("nbinlike: PARAMS must be a two-element vector.");
  endif

  if (params(1) <= 0)
    error (strcat ("nbinlike: number of successes, PARAMS(1), must be", ...
                   " a real positive value."));
  endif

  if (params(2) < 0 || params(2) > 1)
    error (strcat ("nbinlike: probability of success, PARAMS(2), must be", ...
                   " in the range [0,1]."));
  endif

  if (nargin < 3 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("nbinlike: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("nbinlike: FREQ must not contain negative values.");
  elseif (any (fix (freq) != freq))
    error ("nbinlike: FREQ must contain integer values.");
  endif

  ## Expand frequency
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  ## Compute negative log-likelihood and asymptotic variance
  r = params(1);
  ps = params(2);
  nx = numel (x);
  glnr = gammaln (r + x) - gammaln (x + 1) - gammaln (r);
  sumx = sum (x);
  nlogL = -(sum (glnr) + nx * r * log (ps)) - sumx * log (1 - ps);

  if (nargout == 2)
    dL11 = sum (psi (1, r + x) - psi (1, r));
    dL12 = nx ./ ps;
    dL22 = -nx .*r ./ ps .^ 2 - sumx ./ (1 - ps) .^ 2;
    nH = -[dL11, dL12; dL12, dL22];
    if (any (isnan (nH(:))))
        avar = [NaN, NaN; NaN, NaN];
    else
        avar = inv (nH);
    end
  endif

endfunction

## Test output
%!assert (nbinlike ([2.42086, 0.0867043], [1:50]), 205.5942, 1e-4)
%!assert (nbinlike ([3.58823, 0.254697], [1:20]), 63.6435, 1e-4)
%!assert (nbinlike ([8.80671, 0.615565], [1:10]), 24.7410, 1e-4)
%!assert (nbinlike ([22.1756, 0.831306], [1:8]), 17.9528, 1e-4)
%!assert (nbinlike ([22.1756, 0.831306], [1:9], [ones(1,8), 0]), 17.9528, 1e-4)

## Test input validation
%!error<nbinlike: function called with too few input arguments.> nbinlike (3.25)
%!error<nbinlike: X must be a vector.> nbinlike ([5, 0.2], ones (2))
%!error<nbinlike: X cannot have negative values.> nbinlike ([5, 0.2], [-1, 3])
%!error<nbinlike: PARAMS must be a two-element vector.> ...
%! nbinlike ([1, 0.2, 3], [1, 3, 5, 7])
%!error<nbinlike: number of successes,> nbinlike ([-5, 0.2], [1:15])
%!error<nbinlike: number of successes,> nbinlike ([0, 0.2], [1:15])
%!error<nbinlike: probability of success,> nbinlike ([5, 1.2], [3, 5])
%!error<nbinlike: probability of success,> nbinlike ([5, -0.2], [3, 5])
%!error<nbinlike: X and FREQ vectors mismatch.> ...
%! nbinlike ([5, 0.2], ones (10, 1), ones (8,1))
%!error<nbinlike: FREQ must not contain negative values.> ...
%! nbinlike ([5, 0.2], ones (1, 8), [1 1 1 1 1 1 1 -1])
%!error<nbinlike: FREQ must contain integer values.> ...
%! nbinlike ([5, 0.2], ones (1, 8), [1 1 1 1 1 1 1 1.5])
