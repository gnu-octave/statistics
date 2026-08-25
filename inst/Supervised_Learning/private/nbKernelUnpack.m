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
## @deftypefn {Private Function} {@var{DP} =} nbKernelUnpack (@var{DP}, @var{D}, @var{K}, @var{S}, @var{W})
##
## Rebuild the fitted kernel densities @code{nbKernelPack} wrote out as samples.
##
## @var{DP} is a @code{DistributionParameters} cell carrying a sample wherever
## @var{D} names a @qcode{'kernel'} predictor, and @var{K}, @var{S} and @var{W}
## are the model's @code{Kernel}, @code{Support} and @code{Width}.  Each cell is
## refitted through @code{prob.KernelDistribution.fit} with the bandwidth the
## model recorded, so the density comes back identical rather than merely close:
## nothing is re-estimated, since the bandwidth is passed rather than chosen.
##
## A class with no observations of a predictor keeps its empty cell.
##
## @end deftypefn

function DP = nbKernelUnpack (DP, D, K, S, W)

  for j = find (strcmp (D(:)', 'kernel'))
    for i = 1:rows (DP)
      if (! isempty (DP{i,j}) && ! isobject (DP{i,j}))
        DP{i,j} = prob.KernelDistribution.fit (DP{i,j}, K{j}, S{j}, W(i,j));
      endif
    endfor
  endfor

endfunction
