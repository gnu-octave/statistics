## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {statistics} {@var{lp} =} discrimlogp (@var{X}, @var{Mu}, @var{Sigma}, @var{LogDetSigma}, @var{Prior}, @var{DiscrimType})
##
## Log unconditional probability density of every row of @var{X} under a
## discriminant, returned as an @math{Nx1} vector.
##
## The density is @math{P(x) = sum_k P(k) P(x|k)} over the classes, with
## @math{P(k)} the prior and @math{P(x|k)} the multivariate normal density of
## class @math{k}.  The class log densities are formed from the squared
## Mahalanobis distances and @var{LogDetSigma}, so the value follows the
## covariance the model reports rather than a second one derived here, and
## they are summed by their maximum so a well separated observation is a
## density rather than a sum of underflowed zeros.
##
## @end deftypefn

function lp = discrimlogp (X, Mu, Sigma, LogDetSigma, Prior, DiscrimType)

  K = rows (Mu);
  P = columns (X);
  M = discrimmahal (X, Mu, Sigma, DiscrimType);

  logdet = LogDetSigma;
  if (isscalar (logdet))
    logdet = repmat (logdet, K, 1);
  endif

  L = zeros (rows (X), K);
  for k = 1:K
    L(:, k) = -0.5 * (M(:, k) + logdet(k) + P * log (2 * pi)) ...
              + log (Prior(k));
  endfor

  mx = max (L, [], 2);
  lp = mx + log (sum (exp (L - mx), 2));

endfunction
