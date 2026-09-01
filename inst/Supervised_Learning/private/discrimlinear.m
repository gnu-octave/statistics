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
## @deftypefn {Private Function} {[@var{Z}, @var{b}, @var{DeltaPredictor}] =} discrimlinear (@var{Mu}, @var{Sigma}, @var{Prior}, @var{DiscrimType}, @var{Delta})
##
## Per-class linear coefficients of a discriminant, after @var{Delta}.
##
## A linear discriminant scores class @math{k} as @math{Z(k,:) * x' + b(k)},
## which is the Mahalanobis form with the terms common to every class dropped.
## The coefficients are taken about the prior-weighted mean of the class means,
##
## @example
## @var{Z} = (@var{Mu} - @var{Prior} * @var{Mu}) / @var{Sigma}
## @end example
##
## @var{Delta} eliminates a predictor by zeroing coefficients, and it compares
## against the @strong{standardized} coefficient, @code{Z(k,j) * s(j)} with
## @math{s} the within-class standard deviations.  Scaling matters: a threshold
## on the raw coefficients would depend on the units each predictor is measured
## in, so a centimetre and a metre would eliminate different predictors.
##
## @var{DeltaPredictor} is the @var{Delta} at which each predictor drops out
## altogether, @code{max (abs (Z(:,j)) * s(j))} over the classes.
##
## The quadratic family has no linear coefficients to eliminate and does not
## reach here.
##
## @end deftypefn

function [Z, b, DeltaPredictor] = discrimlinear (Mu, Sigma, Prior, ...
                                                 DiscrimType, Delta)

  K = rows (Mu);
  mubar = Prior(:)' * Mu;
  SigmaInv = discriminv (Sigma, DiscrimType);
  Z = (Mu - mubar) * SigmaInv(:,:,1);

  if (strncmp (DiscrimType, 'diag', 4))
    sd = sqrt (Sigma(1,:,1));
  else
    sd = sqrt (diag (Sigma(:,:,1)))';
  endif
  DeltaPredictor = max (abs (Z) .* sd, [], 1);
  if (K == 1)
    DeltaPredictor = abs (Z) .* sd;
  endif

  if (Delta > 0)
    Z(abs (Z) .* sd <= Delta) = 0;
  endif

  ## The constant a shrunken centroid implies.  At Delta of 0 this is the
  ## ordinary linear discriminant constant, the terms it drops being the same
  ## for every class and so lost in the normalization.
  b = zeros (1, K);
  for k = 1:K
    b(k) = -0.5 * (Z(k,:) * (Mu(k,:) + mubar)') + log (Prior(k));
  endfor

endfunction
