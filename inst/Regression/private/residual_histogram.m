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
## @deftypefn {Private Function} {@var{h} =} residual_histogram (@var{r})
##
## Draw the residuals @var{r} with @code{hist} as a histogram of probability
## density, in bins from the @qcode{'scott'} method of @code{histcounts}, and
## return its patch @var{h}.  @qcode{NaN} residuals are left out.
##
## @end deftypefn

function h = residual_histogram (r)

  r = r(! isnan (r));
  [~, edges] = histcounts (r, 'BinMethod', 'scott');
  width = edges(2) - edges(1);
  centers = edges(1:end-1) + width / 2;
  ## HIST reads a scalar as a number of bins
  if (isscalar (centers))
    centers(2) = centers + width;
  endif

  color = [0.1490, 0.5490, 0.8660];
  hist (r, centers, 1 / width, 'FaceColor', color);
  ax = gca ();
  children = get (ax, 'children');
  h = children(1);
  xlabel (ax, 'Residuals');
  ylabel (ax, 'Probability density');
  title (ax, 'Histogram of residuals');

endfunction
