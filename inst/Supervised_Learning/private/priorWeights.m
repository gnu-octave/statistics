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
## @deftypefn {Private Function} priorWeights (@var{Prior}, @var{gY}, @var{n})
##
## Spread each class prior evenly over the observations of that class.
##
## Every observation of class @var{k} carries @qcode{Prior(k) / n_k}, so the
## observations of a class always sum to its prior.  @var{gY} holds the class
## index of each of the @var{n} observations.
##
## @end deftypefn

function w = priorWeights (Prior, gY, n)

  w = zeros (n, 1);
  for k = 1:numel (Prior)
    idx = (gY == k);
    nk = sum (idx);
    if (nk > 0)
      w(idx) = Prior(k) / nk;
    endif
  endfor

endfunction
