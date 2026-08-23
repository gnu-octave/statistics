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
## @deftypefn {Private Function} {@var{P} =} interactionPairs (@var{M})
## The predictor index pairs of the two-way terms in a term matrix.
##
## @var{M} is the term matrix a GAM fits, one row per term and one column per
## predictor, true wherever the term multiplies that predictor.  @var{P} is
## the @math{Kx2} matrix of predictor indices of the rows naming exactly two
## predictors, in the order those rows appear, and @code{zeros (0, 2)} when
## there are none.
##
## A row naming one predictor is a main effect and a row naming three or more
## is a higher-order term.  Neither has a two-column form, so neither appears
## here; the term matrix itself remains the complete record.
## @end deftypefn

function P = interactionPairs (M)

  if (isempty (M))
    P = zeros (0, 2);
    return;
  endif

  M = logical (M);
  keep = find (sum (M, 2) == 2);
  P = zeros (numel (keep), 2);
  for i = 1:numel (keep)
    P(i,:) = find (M(keep(i),:));
  endfor

endfunction
