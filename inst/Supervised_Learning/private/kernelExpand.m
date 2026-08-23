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
## @deftypefn {Private Function} {@var{T} =} kernelExpand (@var{X}, @var{B})
##
## Map predictor data into the expanded feature space of a random basis.
##
## @var{X} is an @math{NxP} numeric matrix and @var{B} a basis structure from
## @code{kernelBasis}.  @var{T} is @math{NxM}, and the inner product of two
## of its rows approximates the Gaussian kernel of the two rows of @var{X}
## they came from.
##
## @end deftypefn

function T = kernelExpand (X, B)

  c = sqrt (2 / B.NumDimensions);
  Z = X * B.Weights;
  T = c * [cos(Z), sin(Z)];
  if (! isempty (B.OddOffset))
    T = [T, c * cos(X * B.OddWeights + B.OddOffset)];
  endif

endfunction
