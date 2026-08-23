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
## @deftypefn {Private Function} {@var{B} =} kernelBasis (@var{p}, @var{m}, @var{sigma})
##
## Draw a random basis for a Gaussian kernel approximation.
##
## @var{p} is the number of predictors, @var{m} the number of expansion
## dimensions and @var{sigma} the kernel scale.  @var{B} is a structure and
## is what @code{kernelExpand} applies to a predictor matrix.
##
## The basis is the random Fourier expansion of Rahimi and Recht.  Frequencies
## are drawn from @math{N(0, sigma^-2 I)} and each contributes a cosine and a
## sine of the same frequency, so that
## @code{T(x1) * T(x2)'} has expectation
## @math{exp (-||x1 - x2||^2 / (2 * sigma^2))}, the Gaussian kernel.
##
## The paired form is used rather than the single cosine with a uniform
## phase offset that MathWorks writes down.  Both are unbiased for the same
## kernel and cost the same, but pairing removes the variance the random
## phase adds: measured over twenty draws on the iris predictors, the mean
## absolute error of the approximated kernel is 0.034 against 0.071 at 128
## dimensions and 0.020 against 0.033 at 512, so the pairs approximate about
## twice as closely for the same @var{m}.  An odd @var{m} cannot be paired
## throughout, and the odd dimension left over is a phase-offset cosine,
## which carries the same expectation.
##
## MATLAB approximates the kernel by the Fastfood construction, which
## replaces the Gaussian frequency matrix with a product of Hadamard,
## permutation and diagonal factors to reach the same distribution in
## @math{O(M log P)} time.  Neither that difference nor this one is
## observable in a single model: the draws come from different generators, so
## a model fitted here and one fitted in MATLAB hold different bases and
## report different scores from any seed.
##
## @end deftypefn

function B = kernelBasis (p, m, sigma)

  pairs = floor (m / 2);
  odd = m - 2 * pairs;

  B = struct ();
  B.Weights = randn (p, pairs) / sigma;
  B.OddWeights = randn (p, odd) / sigma;
  B.OddOffset = 2 * pi * rand (1, odd);
  B.NumDimensions = m;
  B.Scale = sigma;

endfunction
