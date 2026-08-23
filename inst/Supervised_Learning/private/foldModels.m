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
## @deftypefn {Private Function} {@var{Trained} =} foldModels (@var{learner}, @var{F}, @var{Partition}, @var{args})
##
## Fit one model per fold of a partition.
##
## @var{learner} names the class to fit, @var{F} is the frame
## @code{classFrame} returned (or, for a regression model, a structure
## carrying @qcode{X}, @qcode{Y} and @qcode{Weights} alone), @var{Partition}
## a @code{cvpartition} over the retained observations, and @var{args} the
## Name-Value pairs each fold is fitted with.
##
## Every fold is given the parent's class names, prior and cost rather than
## being left to re-derive them from its own rows.  Measured on R2024a: the
## prior of a fold of an unbalanced problem is the parent's, not the fold's
## own frequencies.  A fold that happened to hold one class only would
## otherwise renumber the classes, and every score assembled from it would
## be inverted without anything raising.
##
## What is @emph{not} passed down is any option that resolves to
## @qcode{'auto'}.  Each fold works those out from its own training rows, so
## a fold of a hundred observations partitioned five ways gets
## @qcode{Lambda} of one eightieth rather than one hundredth, and its own
## @qcode{Epsilon}.  Both were measured rather than assumed.
##
## @end deftypefn

function Trained = foldModels (learner, F, Partition, args)

  K = Partition.NumTestSets;
  Trained = cell (K, 1);

  for k = 1:K
    idx = training (Partition, k);
    fargs = args;
    if (! isempty (F.Weights))
      fargs = [fargs, {'Weights', F.Weights(idx)}];
    endif
    Trained{k} = feval (learner, F.X(idx,:), F.Y(idx,:), fargs{:});
  endfor

endfunction
