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
## @deftypefn {Private Function} {@var{L} =} tril_from_theta (@var{th}, @var{q})
##
## Lower-triangular @var{q}-by-@var{q} matrix from its @var{q}*(@var{q}+1)/2
## entries @var{th}, taken in column-major order.
##
## This helper is shared by @code{build_Lfull} and @code{__lme_dfsatt__}.
##
## @end deftypefn

function L = tril_from_theta (th, q)
  L = zeros (q, q);
  idx = 0;
  for j = 1:q
    for i = j:q
      idx += 1;
      L(i, j) = th(idx);
    endfor
  endfor
endfunction
