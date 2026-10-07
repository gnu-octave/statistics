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
## @deftypefn  {Private Function} {[@var{lambda}, @var{lambdahat}] =} __ncfbounds__ (@var{F}, @var{df1}, @var{df2}, @var{alpha})
##
## Confidence bounds for the noncentrality of an F statistic, and its unbiased
## estimate.
##
## @var{lambda} holds the noncentralities at which @var{F} falls at the upper
## and at the lower @math{@var{alpha}/2} point of the noncentral F
## distribution with @var{df1} and @var{df2} degrees of freedom; a bound is
## zero where the central distribution already puts @var{F} beyond that point.
## @code{ncfcdf} loses accuracy at large noncentralities, so where the upper
## bound would pass @math{10^5} both bounds are @qcode{NaN}.
##
## @var{lambdahat} is @math{max (F df1 (df2 - 2) / df2 - df1, 0)}, unbiased
## for the noncentrality before the truncation at zero, and @qcode{NaN} where
## @math{df2 <= 2} leaves the mean of @var{F} undefined.
##
## @end deftypefn

function [lambda, lambdahat] = __ncfbounds__ (F, df1, df2, alpha)

  lmax = 1e5;
  if (df2 > 2)
    lambdahat = max (F * df1 * (df2 - 2) / df2 - df1, 0);
  else
    lambdahat = NaN;
  endif

  lambda = [NaN, NaN];
  if (! isfinite (F) || lambdahat > lmax)
    return;
  endif

  target = [1 - alpha / 2, alpha / 2];
  bounds = [0, 0];
  for j = 1:2
    if (fcdf (F, df1, df2) > target(j))
      hi = max (F * df1, 1);
      while (ncfcdf (F, df1, df2, hi) > target(j))
        hi *= 2;
        if (hi > 2 * lmax)
          return;
        endif
      endwhile
      bounds(j) = fzero (@(l) ncfcdf (F, df1, df2, l) - target(j), [0, hi]);
      if (bounds(j) > lmax)
        return;
      endif
    endif
  endfor
  lambda = bounds;

endfunction
