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
## @deftypefn {Private Function} {@var{DP} =} nbKernelPack (@var{DP}, @var{D})
##
## Replace every fitted kernel density in @var{DP} by the sample it was fitted
## to, so that a naive Bayes model can be written to a file.
##
## @var{DP} is a @code{DistributionParameters} cell, class by predictor, and
## @var{D} the matching @code{DistributionNames}.  A @qcode{'kernel'} predictor
## holds a @code{prob.KernelDistribution} object, which Octave's @code{save}
## cannot serialize: it warns and writes a struct that will not load back as an
## object.  Every other distribution already holds plain numeric data and is
## left alone.
##
## The sample is enough to rebuild the density exactly, because the kernel, the
## support and the bandwidth actually used are stored beside it on the model.
## @code{nbKernelUnpack} does the rebuilding.
##
## @end deftypefn

function DP = nbKernelPack (DP, D)

  for j = find (strcmp (D(:)', 'kernel'))
    for i = 1:rows (DP)
      if (isobject (DP{i,j}))
        DP{i,j} = DP{i,j}.InputData.data;
      endif
    endfor
  endfor

endfunction
