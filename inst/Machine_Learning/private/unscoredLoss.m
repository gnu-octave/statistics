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
## @deftypefn {Private Function} {@var{L} =} unscoredLoss (@var{LossFun}, @var{S}, @var{gY}, @var{w}, @var{Cost}, @var{K})
##
## Classification loss of scores that may hold rows without a score.
##
## A row of @var{S} holding @code{NaN} has no score but still has a label, the
## class @var{K}, so it counts under @qcode{'classiferror'},
## @qcode{'classifcost'} and @qcode{'mincost'}, and makes any other loss, which
## reads the scores, @code{NaN}.  The remaining arguments are those of
## @code{classificationLoss}.
##
## @end deftypefn

function L = unscoredLoss (LossFun, S, gY, w, Cost, K)

  none = any (isnan (S), 2);
  if (any (none))
    if (! any (strcmp (LossFun, {'classiferror', 'classifcost', 'mincost'})))
      L = NaN;
      return;
    endif
    S(none,:) = 0;
    S(none,K) = 1;
  endif
  L = classificationLoss (LossFun, S, gY, w, Cost);

endfunction
