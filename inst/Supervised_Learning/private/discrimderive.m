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
## @deftypefn {statistics} {[@var{Sigma}, @var{LogDetSigma}, @var{Gamma}] =} discrimderive (@var{Base}, @var{DiscrimType}, @var{Gamma}, @var{MinGamma})
##
## Derive a discriminant's covariance from the covariance its fit estimated.
##
## @var{Base} is the unregularized within-class covariance, @math{PxP} for the
## linear family and @math{PxPxK} for the quadratic one.  It is the state a fit
## keeps, and every type in the same family is derived from it, which is what
## makes @code{DiscrimType} assignable after the fit.
##
## @var{Sigma} takes the shape the type calls for: @math{PxP} or @math{PxPxK}
## for the plain and pseudo types, @math{1xP} or @math{1xPxK} for the diagonal
## ones.  @var{LogDetSigma} is a scalar for the linear family and @math{Kx1} for
## the quadratic one.
##
## @var{Gamma} is returned as well as taken, because the type and the
## regularization are one state: a diagonal type is @math{Gamma = 1}, a pseudo
## type is @math{Gamma = 0}, and a plain linear type on a singular covariance is
## raised to @var{MinGamma}, which describes the data and is measured once by
## @code{discrimmingamma}.
##
## @end deftypefn

function [Sigma, LogDetSigma, Gamma] = discrimderive (Base, DiscrimType, ...
                                                      Gamma, MinGamma)

  P = columns (Base);
  K = size (Base, 3);
  LogDetSigma = zeros (K, 1);


  ## Every type takes its log determinant in correlation space, which is not a
  ## nicety: on a rank-3-of-4 iris fixture the direct log (det (Sigma)) gives
  ## -41.2988 where the oracle reports -41.3112, the correlation form matching
  ## to every digit.  Splitting the scale off leaves a matrix of ones on the
  ## diagonal, which is far better conditioned to take a determinant of.

  if (strncmp (DiscrimType, 'diag', 4))
    ## A diagonal type is Gamma = 1: regularizing all the way to the diagonal
    ## and naming the diagonal type are the same state, in both directions.
    ## Its correlation matrix is the identity, so only the scale term remains.
    Gamma = 1;
    Sigma = zeros (1, P, K);
    for k = 1:K
      d = diag (Base(:,:,k))';
      Sigma(1,:,k) = d;
      LogDetSigma(k) = sum (log (d(d > 0)));
    endfor

  elseif (strncmp (DiscrimType, 'pseudo', 6))
    ## The pseudo types never regularize.  They invert whatever rank is there,
    ## so a singular covariance is a case they answer rather than refuse, and
    ## the determinant runs over the directions that carry variance.
    Gamma = 0;
    Sigma = Base;
    for k = 1:K
      [d, R] = correlate (Base(:,:,k));
      ev = eig (R)';
      tol = max (size (R)) * eps (max (ev));
      LogDetSigma(k) = sum (log (d)) + sum (log (ev(ev > tol)));
    endfor

  else
    Sigma = zeros (P, P, K);
    ## MinGamma is the least regularization that makes the correlation matrix
    ## invertible, and a plain linear type is raised to it rather than failing.
    ## The quadratic family is not: it refuses a singular class covariance
    ## outright, so a bump there would only mask the error, and Gamma on a
    ## quadratic discriminant may be nothing but 0 or 1.
    lin = (K == 1);
    for k = 1:K
      S = Base(:,:,k);
      D = diag (diag (S));
      if (lin)
        g = max (Gamma, MinGamma);
      else
        g = Gamma;
      endif
      Sigma(:,:,k) = S * (1 - g) + D * g;
      [d, R] = correlate (Sigma(:,:,k));
      LogDetSigma(k) = sum (log (d)) + log (det (R));
    endfor
    if (lin)
      Gamma = max (Gamma, MinGamma);
    endif
  endif

endfunction

## Split a covariance into its scale and its correlation, keeping only the
## predictors that carry variance.  A predictor with none contributes nothing
## to the determinant rather than driving it to -Inf.
function [d, R] = correlate (S)

  d = diag (S)';
  pos = d > 0;
  d = d(pos);
  if (isempty (d))
    R = [];
    return;
  endif
  R = S(pos, pos) ./ sqrt (d' * d);
  R = (R + R') / 2;

endfunction
