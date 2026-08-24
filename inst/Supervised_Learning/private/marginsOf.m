## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{m} =} marginsOf (@var{s}, @var{gY}, @var{L})
##
## Classification margin of every observation, over one or more score
## matrices.
##
## @var{s} is @math{NxK} when @var{L} is one and @math{NxKxL} otherwise,
## @var{gY} the index into the class names of each observation's true class.
## @var{m} is @math{NxL}: the score of the true class less the best score
## among the others, which for two classes is simply the other one.
##
## @end deftypefn

function m = marginsOf (s, gY, L)

  n = rows (s);
  m = zeros (n, L);
  rowidx = (1:n)';

  for k = 1:L
    if (L == 1)
      sk = s;
    else
      sk = s(:,:,k);
    endif
    idx = sub2ind (size (sk), rowidx, gY);
    strue = sk(idx);
    so = sk;
    so(idx) = -Inf;
    m(:,k) = strue - max (so, [], 2);
  endfor

endfunction
