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
## @deftypefn  {Private Function} {@var{h} =} residual_histogram (@var{ax}, @var{r})
## @deftypefnx {Private Function} {@var{h} =} residual_histogram (@var{ax}, @var{r}, @var{props})
##
## Draw the residuals @var{r} into the axes @var{ax} as a histogram of
## probability density, in bins from the @qcode{'scott'} method of
## @code{histcounts}, and return its patch @var{h}.  @qcode{NaN} residuals are
## left out.  @var{props} is a cell array of patch properties; without it the
## bars take the default blue face.
##
## @end deftypefn

function h = residual_histogram (ax, r, props)

  if (nargin < 3)
    props = {'FaceColor', [0.1490, 0.5490, 0.8660]};
  endif

  r = r(! isnan (r));
  [counts, edges] = histcounts (r, 'BinMethod', 'scott');
  dens = counts / (numel (r) * (edges(2) - edges(1)));
  nb = numel (counts);
  left = edges(1:end-1);
  right = edges(2:end);
  Xp = [left; left; right; right];
  Yp = [zeros(1, nb); dens; dens; zeros(1, nb)];
  ## A new plot replaces the old one, unless 'hold' is on
  ax = newplot (ax);
  h = patch (ax, 'XData', Xp, 'YData', Yp, props{:});
  xlabel (ax, 'Residuals');
  ylabel (ax, 'Probability density');
  title (ax, 'Histogram of residuals');

endfunction
