## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{yFit} =} kfoldResponse (@var{Trained}, @var{Partition}, @var{X}, @var{L})
##
## Assemble the out-of-fold predictions of a cross-validated regression
## model.
##
## Each observation is predicted by the fold that held it out, and an
## observation no fold held out comes back @qcode{NaN} rather than
## predicted.  Under a holdout partition that is most of them.  @var{yFit}
## is @math{NxL}, one column per regularization strength.
##
## The predictions are untransformed: the fold models carry no response
## transform and the parent applies its own once.
##
## @end deftypefn

function yFit = kfoldResponse (Trained, Partition, X, L)

  yFit = nan (rows (X), L);

  for k = 1:Partition.NumTestSets
    idx = test (Partition, k);
    if (! any (idx))
      continue;
    endif
    yFit(idx,:) = predict (Trained{k}, X(idx,:));
  endfor

endfunction
