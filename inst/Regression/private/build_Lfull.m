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
## @deftypefn {Private Function} {@var{Lf} =} build_Lfull (@var{theta}, @var{qk}, @var{nlev})
##
## Block-diagonal relative covariance factor of a mixed model.
##
## @var{theta} holds, term after term, the q_k*(q_k+1)/2 column-major
## lower-triangle entries of the factor L_k of each term's relative covariance.
## @var{Lf} is the sparse block-diagonal matrix with L_k repeated over the
## @var{nlev}(k) levels of term k, held sparse so that products against it stay
## inside the block structure.
##
## This helper is shared by @code{__lmefit__} and @code{__glmefit__}.
##
## @end deftypefn

function Lf = build_Lfull (theta, qk, nlev)
  blocks = cell (1, numel (qk));
  off = 0;
  for k = 1:numel (qk)
    m = qk(k)*(qk(k)+1)/2;
    Lk = sparse (tril_from_theta (theta(off+(1:m)), qk(k)));
    blocks{k} = kron (speye (nlev(k)), Lk);
    off += m;
  endfor
  Lf = blkdiag (blocks{:});
endfunction
