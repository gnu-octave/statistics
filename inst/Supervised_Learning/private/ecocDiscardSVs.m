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
## @deftypefn {Private Function} {[@var{learners}, @var{n}] =} ecocDiscardSVs (@var{learners})
##
## Discard the support vectors of every binary learner that has any.
##
## @var{learners} is the cell of binary learners of an error correcting code.
## Only a support vector machine on a linear kernel can give its vectors up,
## the linear model standing in for them exactly; any other learner is left
## as it is rather than refused, a code being free to mix them.
##
## @var{n} is how many learners were collapsed, which is zero when none of
## them was a linear support vector machine and is what the caller warns on.
##
## @seealso{CompactClassificationECOC, ClassificationECOC}
## @end deftypefn

function [learners, n] = ecocDiscardSVs (learners)

  n = 0;
  svm = {'ClassificationSVM', 'CompactClassificationSVM'};
  for j = 1:numel (learners)
    L = learners{j};
    if (! any (strcmp (class (L), svm)))
      continue;
    endif
    if (! strcmpi (L.KernelParameters.Function, 'linear'))
      continue;
    endif
    learners{j} = discardSupportVectors (L);
    n++;
  endfor

endfunction
