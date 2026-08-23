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
## @deftypefn {Private Function} {@var{sigma} =} autoKernelScale (@var{X})
##
## The kernel scale a kernel learner resolves @qcode{'auto'} to.
##
## @end deftypefn

## The kernel scale MATLAB calls 'auto'.  MathWorks documents a heuristic
## that subsamples and warns that its estimates vary from one call to the
## next, so no particular number is the parity target.  This one is the
## median distance between observations, which is the usual choice and is
## reproducible: a kernel that wide puts about half the pairs inside one
## bandwidth of each other.  Beyond a thousand observations the median is
## taken over a random thousand of them, the full matrix being quadratic.
function sigma = autoKernelScale (X)

  n = rows (X);
  if (n > 1000)
    X = X(randperm (n, 1000), :);
    n = 1000;
  endif
  if (n < 2)
    sigma = 1;
    return;
  endif

  d = zeros (n * (n - 1) / 2, 1);
  at = 0;
  for i = 1:n-1
    k = n - i;
    d(at+1:at+k) = sqrt (sum ((X(i+1:end,:) - X(i,:)) .^ 2, 2));
    at += k;
  endfor
  sigma = median (d);
  if (! (sigma > 0))
    sigma = 1;
  endif

endfunction
