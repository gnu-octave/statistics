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
## @deftypefn {Private Function} {@var{m} =} weightedMedian (@var{Y}, @var{W})
##
## The weighted median of @var{Y}, the value at which half the weight lies
## on each side.
##
## @end deftypefn

## The weighted median, the value at which half the weight lies on each side.
## It is what a fit of the epsilon-insensitive loss starts its intercept
## from, that loss being an absolute deviation once the band is crossed.
function m = weightedMedian (Y, W)

  [ys, k] = sort (Y);
  ws = W(k);
  c = cumsum (ws) / sum (ws);
  j = find (c >= 0.5, 1);
  m = ys(j);

endfunction
