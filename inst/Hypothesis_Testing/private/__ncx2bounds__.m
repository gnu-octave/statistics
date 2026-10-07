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
## @deftypefn  {Private Function} {@var{lambda} =} __ncx2bounds__ (@var{chisq}, @var{df}, @var{alpha})
##
## Confidence bounds for the noncentrality of a chi-square statistic.
##
## @var{lambda} holds the noncentralities at which @var{chisq} falls at the
## upper and at the lower @math{@var{alpha}/2} point of the noncentral
## chi-square distribution with @var{df} degrees of freedom; a bound is zero
## where the central distribution already puts @var{chisq} beyond that point.
##
## @end deftypefn

function lambda = __ncx2bounds__ (chisq, df, alpha)

  target = [1 - alpha / 2, alpha / 2];
  lambda = [0, 0];
  for j = 1:2
    if (chi2cdf (chisq, df) > target(j))
      hi = max (chisq, 1);
      while (ncx2cdf (chisq, df, hi) > target(j))
        hi *= 2;
      endwhile
      lambda(j) = fzero (@(l) ncx2cdf (chisq, df, l) - target(j), [0, hi]);
    endif
  endfor

endfunction
