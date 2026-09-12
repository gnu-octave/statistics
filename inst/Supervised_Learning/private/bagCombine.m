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
## @deftypefn {Private Function} {[@var{A}, @var{SD}, @var{none}] =} bagCombine (@var{P}, @var{use}, @var{tw}, @var{mode})
##
## Weighted mean and standard deviation of a bagged ensemble's tree outputs.
##
## @var{P} is the @math{NxKxT} array @code{bagTreeOutputs} returns, @var{use}
## an @math{NxT} logical matrix of which tree may answer for which
## observation, and @var{tw} a weight for each tree.  @var{mode} is
## @qcode{'ensemble'}, averaging over every tree at once, @qcode{'cumulative'},
## averaging over the first @math{t} trees for each @math{t}, or
## @qcode{'individual'}, taking each tree on its own.
##
## @var{A} and @var{SD} are @math{NxKxL}, with @math{L} one for
## @qcode{'ensemble'} and @math{T} otherwise.  @var{SD} is the weighted
## population standard deviation over the trees used.  @var{none} is an
## @math{NxL} logical matrix marking the observations no tree could answer
## for, whose entries of @var{A} are not numbers.
##
## @end deftypefn

function [A, SD, none] = bagCombine (P, use, tw, mode)

  [N, K, T] = size (P);
  W3 = reshape (double (use) .* tw, N, 1, T);
  switch (mode)
    case 'ensemble'
      den = sum (W3, 3);
      num = sum (P .* W3, 3);
      sq = sum ((P .^ 2) .* W3, 3);
    case 'cumulative'
      den = cumsum (W3, 3);
      num = cumsum (P .* W3, 3);
      sq = cumsum ((P .^ 2) .* W3, 3);
    otherwise
      den = W3;
      num = P .* W3;
      sq = (P .^ 2) .* W3;
  endswitch
  none = reshape (den == 0, N, size (den, 3));
  A = num ./ den;
  SD = sqrt (max (sq ./ den - A .^ 2, 0));
  SD(isnan (A)) = NaN;

endfunction
