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
## @deftypefn {Private Function} {@var{ZtZ} =} choose_storage (@var{ZtZ})
##
## Choose sparse or full storage for the random-effects cross product
## @var{ZtZ}, which @code{factor_K} follows.
##
## @code{K = I + L'*ZtZ*L} has the non-zero pattern of @var{ZtZ} whatever the
## covariance parameters are, so the choice is made once, from the fill of its
## Cholesky factor estimated under an approximate minimum degree order.  Sparse
## is faster up to about a fifth of the triangle filled, as with nested or
## few-level crossed terms; crossed terms with many levels each fill in, and
## @var{ZtZ} is returned full.
##
## This helper is shared by @code{__lmefit__}, @code{__glmefit__} and
## @code{__lme_dfsatt__}.
##
## @end deftypefn

function ZtZ = choose_storage (ZtZ)
  q = columns (ZtZ);
  P = spones (ZtZ) + speye (q);
  o = amd (P);
  if (sum (symbfact (P(o,o))) > 0.2 * q * (q + 1) / 2)
    ZtZ = full (ZtZ);
  endif
endfunction
