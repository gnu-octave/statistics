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
## @deftypefn {Private Function} {[@var{labels}, @var{scores}] =} kfoldScores (@var{Trained}, @var{Partition}, @var{X}, @var{Y}, @var{ClassNames}, @var{L})
##
## Assemble the out-of-fold labels and scores of a cross-validated binary
## classifier.
##
## Each observation is predicted by the fold that held it out, and an
## observation no fold held out is left missing rather than given a class.
## Under a holdout partition that is most of them.  @var{Y} is the response,
## used only for its type, so that the labels come back in it.
##
## @var{L} is the number of regularization strengths the fold models carry,
## which is one for a kernel model and one per value of @qcode{'Lambda'} for
## a linear one.  @var{scores} is @math{Nx2xL}, or @math{Nx2} when @var{L}
## is one.
##
## The scores are whatever the fold models report, transform included: a
## logistic fold returns posteriors.  MATLAB leaves the transform with the
## folds and has the cross-validated model report none of its own, so a
## parent that transformed again would apply it twice.
##
## @end deftypefn

function [labels, scores] = kfoldScores (Trained, Partition, X, Y, ...
                                         ClassNames, L)

  n = rows (X);
  K = Partition.NumTestSets;

  if (iscellstr (Y))
    labels = repmat ({''}, n, L);
  elseif (islogical (Y))
    labels = false (n, L);
  elseif (ischar (Y))
    labels = repmat (' ', n, columns (ClassNames));
  else
    labels = nan (n, L);
  endif
  scores = nan (n, 2, L);

  for k = 1:K
    idx = test (Partition, k);
    if (! any (idx))
      continue;
    endif
    [lab, sc] = predict (Trained{k}, X(idx,:));
    if (ischar (Y))
      labels(idx,:) = lab;
    else
      labels(idx,:) = lab;
    endif
    if (L == 1)
      scores(idx,:,1) = sc;
    else
      scores(idx,:,:) = sc;
    endif
  endfor

  if (L == 1)
    scores = scores(:,:,1);
  endif

endfunction
