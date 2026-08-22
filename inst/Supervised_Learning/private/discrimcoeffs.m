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
## @deftypefn {statistics} {@var{Coeffs} =} discrimcoeffs (@var{Mu}, @var{Sigma}, @var{LogDetSigma}, @var{Prior}, @var{DiscrimType}, @var{Delta}, @var{ClassNames})
##
## Pairwise boundary coefficients of a discriminant, a @math{KxK} structure
## whose @math{(i,j)} entry separates class @math{i} from class @math{j}.
##
## The linear family reports @code{Const} and @code{Linear}; the quadratic one
## adds @code{Quadratic}, which follows @var{Sigma}'s shape and so is
## @math{1xP} for @qcode{'diagQuadratic'}.  The diagonal entries carry the two
## class names and nothing else.
##
## @end deftypefn

function Coeffs = discrimcoeffs (Mu, Sigma, LogDetSigma, Prior, ...
                                 DiscrimType, Delta, ClassNames)

  K = rows (Mu);
  [~, fam] = discrimcanon (DiscrimType);
  isquad = strcmp (fam, 'quadratic');
  isdiag = strncmp (DiscrimType, 'diag', 4);
  SigmaInv = discriminv (Sigma, DiscrimType);
  if (! isquad)
    [Z, bk] = discrimlinear (Mu, Sigma, Prior, DiscrimType, Delta);
  endif

  Coeffs = struct ();
  for i = 1:K
    for j = 1:K
      Coeffs(i,j).DiscrimType = '';
      Coeffs(i,j).Const = [];
      Coeffs(i,j).Linear = [];
      if (isquad)
        Coeffs(i,j).Quadratic = [];
      endif
      Coeffs(i,j).Class1 = ClassNames(i,:);
      Coeffs(i,j).Class2 = ClassNames(j,:);
      if (i == j)
        continue;
      endif

      Coeffs(i,j).DiscrimType = DiscrimType;
      if (isquad)
        Si = SigmaInv(:,:,i);
        Sj = SigmaInv(:,:,j);
        L = Si * Mu(i,:)' - Sj * Mu(j,:)';
        C = -0.5 * (Mu(i,:) * Si * Mu(i,:)' - Mu(j,:) * Sj * Mu(j,:)') ...
            - 0.5 * (LogDetSigma(i) - LogDetSigma(j)) ...
            + log (Prior(i) / Prior(j));
        Q = -0.5 * (Si - Sj);
        if (isdiag)
          ## A diagonal model reports the diagonal, matching Sigma's own shape.
          Q = diag (Q)';
        endif
        Coeffs(i,j).Quadratic = Q;
      else
        ## Taken from the per-class coefficients, so a predictor Delta has
        ## eliminated is eliminated from every pair that names it.
        L = (Z(i,:) - Z(j,:))';
        C = bk(i) - bk(j);
      endif
      Coeffs(i,j).Linear = L;
      Coeffs(i,j).Const = C;
    endfor
  endfor

endfunction
