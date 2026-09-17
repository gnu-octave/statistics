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
## @deftypefn {Private Function} {@var{out} =} bagOutlier (@var{P}, @var{g})
##
## Outlier measure of each observation from a proximity matrix.
##
## @var{P} is an @math{NxN} proximity matrix and @var{g} a class index for
## each observation, or empty to treat every observation as one class.  The
## raw measure of an observation is the size of its class divided by the sum
## of its squared proximities to the observations of the class, itself
## included.  @var{out} is the absolute deviation of the raw measure from the
## median of its class, divided by the median absolute deviation of the
## class.  When that deviation is zero, @var{out} is the raw measure itself,
## and in a class of one or two observations it is zero.
##
## @end deftypefn

function out = bagOutlier (P, g)

  N = rows (P);
  if (isempty (g))
    g = ones (N, 1);
  endif
  g = g(:);
  P2 = P .^ 2;
  out = zeros (N, 1);
  for k = unique (g)'
    idx = g == k;
    n = sum (idx);
    if (n <= 2)
      continue;
    endif
    raw = n ./ sum (P2(idx,idx), 2);
    dev = abs (raw - median (raw));
    mdev = median (dev);
    if (mdev == 0)
      out(idx) = raw;
    else
      out(idx) = dev / mdev;
    endif
  endfor

endfunction
