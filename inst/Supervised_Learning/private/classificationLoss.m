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
## @deftypefn {Private Function} {@var{l} =} classificationLoss (@var{LossFun}, @var{s}, @var{gY}, @var{w}, @var{Cost})
##
## One of MATLAB's classification losses, over an @math{NxK} score matrix.
##
## @var{s} holds the scores after @qcode{ScoreTransform}, @var{gY} the index
## into the class names of each observation's true class, @var{w} weights
## summing to one, and @var{Cost} the misclassification cost matrix.
##
## The rival of the true class is the best of the others, which for two
## classes is simply the other one, so a binary caller is unaffected.
##
## Every loss but the two cost-based ones is a function of the score the
## model gives the @emph{true} class, and not of the margin.  Measured on
## R2024a across both learners: with an untransformed support vector machine
## the two readings differ by a factor of two, and with a logit transformed
## logistic regression they differ by an affine map, and the true-class score
## is the one that reproduces all six losses in both cases.
##
## @end deftypefn

## The classification losses MATLAB offers, all of them functions of the
## score the model gives the true class.
function l = classificationLoss (LossFun, s, gY, w, Cost)

  n = rows (s);
  idx = sub2ind (size (s), (1:n)', gY);
  strue = s(idx);
  ## The best score among the classes that are not the true one
  so = s;
  so(idx) = -Inf;
  [sother, kother] = max (so, [], 2);

  switch (LossFun)

    case 'classiferror'
      l = sum (w .* (strue <= sother));

    case 'hinge'
      l = sum (w .* max (0, 1 - strue));

    case 'quadratic'
      l = sum (w .* ((1 - strue) .^ 2));

    case 'logit'
      l = sum (w .* log (1 + exp (-strue)));

    case 'binodeviance'
      l = sum (w .* log (1 + exp (-2 * strue)));

    case 'exponential'
      l = sum (w .* exp (-strue));

    case 'mincost'
      expected = s * Cost;
      [~, k] = min (expected, [], 2);
      l = sum (w .* Cost(sub2ind (size (Cost), gY, k)));

    case 'classifcost'
      ## The class the model would pick, the true one keeping a tie
      k = gY;
      worse = sother > strue;
      k(worse) = kother(worse);
      l = sum (w .* Cost(sub2ind (size (Cost), gY, k)));

  endswitch

endfunction
