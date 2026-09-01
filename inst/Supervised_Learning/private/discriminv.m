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
## @deftypefn {Private Function} {@var{SigmaInv} =} discriminv (@var{Sigma}, @var{DiscrimType})
##
## Inverse within-class covariance of a discriminant, one @math{PxP} slice per
## class, whatever shape @var{Sigma} is stored in.
##
## The diagonal types invert elementwise, leaving a zero where a predictor has
## no variance, so that predictor contributes nothing rather than an infinity.
## The pseudo types take a pseudo-inverse in correlation space, which is the
## same convention their @code{LogDetSigma} uses.
##
## @end deftypefn

function SigmaInv = discriminv (Sigma, DiscrimType)

  K = size (Sigma, 3);

  if (strncmp (DiscrimType, 'diag', 4))
    P = columns (Sigma);
    SigmaInv = zeros (P, P, K);
    for k = 1:K
      d = Sigma(1,:,k);
      v = zeros (1, P);
      v(d > 0) = 1 ./ d(d > 0);
      SigmaInv(:,:,k) = diag (v);
    endfor

  elseif (strncmp (DiscrimType, 'pseudo', 6))
    P = columns (Sigma);
    SigmaInv = zeros (P, P, K);
    for k = 1:K
      SigmaInv(:,:,k) = pseudoinverse (Sigma(:,:,k));
    endfor

  else
    SigmaInv = zeros (size (Sigma));
    for k = 1:K
      SigmaInv(:,:,k) = inv (Sigma(:,:,k));
    endfor
  endif

endfunction

## Pseudo-inverse taken on the correlation matrix and scaled back, the same
## convention the pseudo log determinant uses.  A predictor with no variance is
## left out and comes back as a zero row and column, which is how the oracle
## reports its coefficient.
function Si = pseudoinverse (S)

  P = columns (S);
  Si = zeros (P, P);
  d = diag (S)';
  pos = d > 0;
  if (! any (pos))
    return;
  endif
  dp = d(pos);
  R = S(pos, pos) ./ sqrt (dp' * dp);
  Si(pos, pos) = pinv ((R + R') / 2) ./ sqrt (dp' * dp);

endfunction
