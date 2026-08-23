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
## @deftypefn {Private Function} {@var{sets} =} foldSets (@var{Partition}, @var{Folds}, @var{Mode}, @var{n})
##
## The sets of observations a @code{kfold} method reports over.
##
## Under @qcode{'individual'} there is one set per named fold, holding the
## observations that fold held out.  Under @qcode{'average'} there is one
## set holding the observations of all the named folds together.
##
## The distinction matters, and it is not the one the name suggests:
## averaging pools the observations rather than averaging the per-fold
## values, so folds of unequal size do not weigh equally.  Measured on
## R2024a: a four-fold cross-validation of carsmall, whose folds hold 24,
## 23, 23 and 23 observations, reports 64.0992586069075 where the mean of
## its own four fold losses is 64.2539.  The first is the mean squared error
## over all ninety-three out-of-fold predictions, and that is what MATLAB
## returns.
##
## @end deftypefn

function sets = foldSets (Partition, Folds, Mode, n)

  if (strcmp (Mode, 'individual'))
    sets = cell (numel (Folds), 1);
    for i = 1:numel (Folds)
      sets{i} = test (Partition, Folds(i));
    endfor
  else
    pooled = false (n, 1);
    for i = 1:numel (Folds)
      pooled = pooled | test (Partition, Folds(i));
    endfor
    sets = {pooled};
  endif

endfunction
