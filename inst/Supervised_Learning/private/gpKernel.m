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
## @deftypefn {Private Function} gpKernel (@var{A}, @var{B}, @var{name}, @var{theta})
##
## Evaluate a Gaussian process covariance function between two sets of points.
##
## @var{A} is @math{NxD} and @var{B} is @math{MxD}, and the returned matrix is
## @math{NxM} holding the covariance of every pair.  @var{name} is one of the
## ten built-in covariance functions, spelled as the model stores it, and
## @var{theta} is its parameter vector.
##
## The isotropic kernels take @qcode{[SigmaL; SigmaF]}, the rational quadratic
## @qcode{[SigmaL; AlphaRQ; SigmaF]}, the automatic relevance determination
## kernels one length scale per predictor followed by @qcode{SigmaF}, and the
## ARD rational quadratic those length scales, @qcode{AlphaRQ} and
## @qcode{SigmaF}.  Every kernel is stationary, so the covariance depends on
## the two points only through the distance between them, scaled by the length
## scale of each predictor.
##
## @end deftypefn

function K = gpKernel (A, B, name, theta)

  d = columns (A);
  ard = (numel (name) > 3 && strncmpi (name, 'ARD', 3));

  ## Split the parameter vector and scale each predictor by its length scale,
  ## which reduces every kernel below to a function of the scaled distance.
  if (ard)
    L = theta(1:d)(:)';
    if (strcmpi (name, 'ARDRationalQuadratic'))
      alphaRQ = theta(d+1);
      sigmaF = theta(d+2);
    else
      sigmaF = theta(d+1);
    endif
  else
    L = repmat (theta(1), 1, d);
    if (strcmpi (name, 'RationalQuadratic'))
      alphaRQ = theta(2);
      sigmaF = theta(3);
    else
      sigmaF = theta(2);
    endif
  endif

  ## Squared distances, accumulated one predictor at a time.  The expanded
  ## form, sum (A.^2) - 2*A*B' + sum (B.^2), is faster but does not return
  ## exactly zero for a point against itself: it leaves a residue of the order
  ## of eps, and the square root of 1.8e-15 is 4.2e-8, so a stationary kernel
  ## stops returning exactly SigmaF^2 on its own diagonal.  The rough kernels
  ## take the square root and inherit that as an error eight orders larger
  ## than the one it came from.
  r2 = zeros (rows (A), rows (B));
  for j = 1:d
    r2 += ((A(:,j) - B(:,j)') / L(j)) .^ 2;
  endfor

  switch (lower (strrep (lower (name), 'ard', '')))

    case 'exponential'
      K = sigmaF^2 * exp (-sqrt (r2));

    case 'squaredexponential'
      K = sigmaF^2 * exp (-0.5 * r2);

    case 'matern32'
      s = sqrt (3) * sqrt (r2);
      K = sigmaF^2 * (1 + s) .* exp (-s);

    case 'matern52'
      s = sqrt (5) * sqrt (r2);
      K = sigmaF^2 * (1 + s + (5/3) * r2) .* exp (-s);

    case 'rationalquadratic'
      K = sigmaF^2 * (1 + r2 / (2 * alphaRQ)) .^ (-alphaRQ);

    otherwise
      error ("gpKernel: unrecognized kernel function '%s'.", name);

  endswitch

endfunction
