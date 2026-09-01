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
## @deftypefn {Private Function} {@var{g} =} discrimmingamma (@var{S})
##
## Least regularization that leaves a covariance invertible.
##
## @var{g} is the smallest multiple of @code{eps} for which
## @code{(1 - @var{g}) * R + @var{g} * I} is invertible, where @var{R} is
## @var{S}'s correlation matrix, and 0 when @var{R} already is.  Multiples of
## @code{eps} are the quantization the oracle's own readings show, 0 on a well
## conditioned fit and 8 and 10 times @code{eps} on two singular ones.
##
## It describes the data rather than any one discriminant type, so it is taken
## once from the pooled within-class covariance and reported by every type,
## whether or not that type regularizes.
##
## @end deftypefn

function g = discrimmingamma (S)

  g = 0;
  d = diag (S)';
  pos = d > 0;
  if (! any (pos))
    return;
  endif
  dp = d(pos);
  R = S(pos, pos) ./ sqrt (dp' * dp);
  R = (R + R') / 2;
  n = rows (R);
  I = eye (n);
  if (rcond (R) > n * eps)
    return;
  endif
  for m = 1:64
    gm = m * eps;
    if (rcond (R * (1 - gm) + I * gm) > n * eps)
      g = gm;
      return;
    endif
  endfor

endfunction
