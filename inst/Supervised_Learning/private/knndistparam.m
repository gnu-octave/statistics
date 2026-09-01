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

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{p} =} knndistparam (@var{Distance}, @var{X}, @var{standardized})
##
## Distance parameter a nearest-neighbour metric implies.
##
## Three metrics carry one: @qcode{'minkowski'} the exponent, defaulting to 2;
## @qcode{'seuclidean'} a per-predictor scale, the standard deviations of
## @var{X}, or ones where the predictors were standardized already; and
## @qcode{'mahalanobis'} the covariance of @var{X}.  Every other metric has
## none and @var{p} is empty.
##
## It is recomputed rather than carried over whenever @code{Distance} changes,
## which is what the oracle does: a scale belonging to one metric means
## nothing under another.
##
## @end deftypefn

function p = knndistparam (Distance, X, standardized)

  p = [];
  if (! ischar (Distance))
    return;
  endif
  switch (lower (Distance))
    case 'minkowski'
      p = 2;
    case 'seuclidean'
      if (standardized)
        p = ones (1, columns (X));
      else
        p = std (X, [], 1);
      endif
    case 'mahalanobis'
      p = cov (X);
  endswitch

endfunction
