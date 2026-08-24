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
##

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{L} =} nbLevels (@var{x})
##
## The distinct levels of one categorical predictor, in ascending order.
##
## The order is the one every per-level quantity is reported in, so it is
## fixed here once rather than recovered wherever a level is looked up.
##
## @seealso{ClassificationNaiveBayes}
## @end deftypefn

function L = nbLevels (x)

  L = unique (x(:));
  L = L(! isnan (L));

endfunction
