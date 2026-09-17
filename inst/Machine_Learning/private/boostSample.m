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
## @deftypefn {Private Function} {[@var{idx}, @var{sw}, @var{cnt}] =} boostSample (@var{w}, @var{m}, @var{replace})
## The rows a resampling boosting method grows one learner on.
##
## With @var{replace}, @var{m} rows are drawn with replacement in proportion to
## the weights @var{w}, each draw carrying the same sample weight; without it,
## @var{m} rows are drawn uniformly without replacement, each carrying its
## weight from @var{w}.  @var{idx} lists the rows drawn, @var{sw} their sample
## weights, which sum to one, and @var{cnt} how many times each row was drawn.
## These are the draws MATLAB R2024a makes.
## @end deftypefn

function [idx, sw, cnt] = boostSample (w, m, replace)

  w = w(:);
  n = numel (w);
  if (replace)
    cw = [0; cumsum(w)];
    cw /= cw(end);
    idx = lookup (cw, rand (m, 1));
    sw = ones (m, 1) / m;
  else
    order = randperm (n);
    idx = sort (order(1:m));
    idx = idx(:);
    sw = w(idx);
    if (sum (sw) > 0)
      sw /= sum (sw);
    else
      sw = ones (m, 1) / m;
    endif
  endif
  cnt = accumarray (idx, 1, [n, 1]);

endfunction
