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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{Rk}, @var{flag}, @var{pk}] =} factor_K (@var{L}, @var{ZtZ})
##
## Cholesky factor of the mixed-model matrix @code{K = I + L'*ZtZ*L}.
##
## @var{Rk} is taken in the row order @var{pk}, so that @code{Rk'*Rk} is
## @code{K(pk,pk)}, and @var{flag} is nonzero when @var{K} is not positive
## definite.  @var{K} is factorised sparse, under a fill-reducing order, when
## @var{ZtZ} is sparse, and full otherwise, with @var{pk} the identity order.
## @var{K} is symmetric in exact arithmetic but not bitwise, the two products
## being separate calls, and @code{chol} reads one triangle, so it is
## symmetrised first.
##
## This helper is shared by @code{__lmefit__}, @code{__glmefit__} and
## @code{__lme_dfsatt__}.
##
## @end deftypefn

function [Rk, flag, pk] = factor_K (L, ZtZ)
  K = speye (columns (ZtZ)) + L' * ZtZ * L;
  K = (K + K') / 2;
  if (issparse (K))
    [Rk, flag, pk] = chol (K, "vector");
  else
    [Rk, flag] = chol (K);
    pk = 1:columns (K);
  endif
endfunction
