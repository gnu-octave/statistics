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
## @deftypefn {statistics} {@var{M} =} discrimmahal (@var{X}, @var{Mu}, @var{Sigma}, @var{DiscrimType})
##
## Squared Mahalanobis distance from every row of @var{X} to every class mean
## of a discriminant, returned as an @math{NxK} matrix.
##
## @var{M}(i,j) is @math{(x_i - mu_j)' * inv (S_j) * (x_i - mu_j)}, where
## @math{S_j} is the class covariance for a quadratic type and the one shared
## covariance for a linear one.  The inverse comes from @code{discriminv}, so
## the regularized covariance the model reports is the one measured against,
## and the diagonal and pseudo types keep their own conventions.
##
## @end deftypefn

function M = discrimmahal (X, Mu, Sigma, DiscrimType)

  K = rows (Mu);
  SigmaInv = discriminv (Sigma, DiscrimType);
  nInv = size (SigmaInv, 3);
  M = zeros (rows (X), K);
  for k = 1:K
    Si = SigmaInv(:,:,min (k, nInv));
    Z = X - Mu(k, :);
    M(:, k) = sum ((Z * Si) .* Z, 2);
  endfor

endfunction
