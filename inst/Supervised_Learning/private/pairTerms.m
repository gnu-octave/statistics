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
## @deftypefn {Private Function} {@var{T} =} pairTerms (@var{p})
## Every pairwise interaction of @var{p} predictors, as a logical matrix with
## one row per term and one column per predictor, two of them true.
## @end deftypefn

function T = pairTerms (p)

  if (p < 2)
    T = false (0, p);
    return;
  endif
  pairs = nchoosek (1:p, 2);
  n = rows (pairs);
  T = false (n, p);
  T(sub2ind ([n, p], (1:n)', pairs(:,1))) = true;
  T(sub2ind ([n, p], (1:n)', pairs(:,2))) = true;

endfunction
