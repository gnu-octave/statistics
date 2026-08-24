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
##

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{w} =} priorNormalize (@var{w}, @var{gY}, @var{Prior})
##
## Scale observation weights so that each class carries its prior.
##
## Within a class the weights keep their relative sizes; across classes each
## class's weights sum to that class's prior.  A loss computed on weights that
## were not normalized this way would let an over-represented class speak for
## more than its prior says it should.
##
## A class with no weight at all is left alone rather than divided by zero.
##
## @seealso{ClassificationNaiveBayes, edgeWeights}
## @end deftypefn

function w = priorNormalize (w, gY, Prior)

  w = w(:);
  for i = 1:numel (Prior)
    idx = gY == i;
    sw = sum (w(idx));
    if (sw > 0)
      w(idx) = w(idx) * Prior(i) / sw;
    endif
  endfor
  sw = sum (w);
  if (sw > 0)
    w = w / sw;
  endif

endfunction
