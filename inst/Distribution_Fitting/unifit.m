## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{ahat} =} unifit (@var{x})
## @deftypefnx {statistics} {[@var{ahat}, @var{bhat}] =} unifit (@var{x})
## @deftypefnx {statistics} {[@var{ahat}, @var{bhat}, @var{aci}, @var{bci}] =} unifit (@var{x})
## @deftypefnx {statistics} {[@dots{}] =} unifit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} unifit (@var{x}, @var{alpha}, @var{freq})
##
## Estimate parameters and confidence intervals for the continuous uniform
## distribution.
##
## @code{[@var{ahat}, @var{bhat}] = unifit (@var{x})} returns the maximum
## likelihood estimates of the lower and upper endpoints, @var{a} and @var{b},
## of the continuous uniform distribution given the data in @var{x}.  Each
## estimate is returned as a separate output.
##
## @var{x} may be a vector, which is fitted as a single sample, or a matrix,
## which is fitted column by column.  For a matrix of @math{n} columns
## @var{ahat} and @var{bhat} are @math{1}-by-@math{n} row vectors and @var{aci}
## and @var{bci} are @math{2}-by-@math{n}.
##
## @code{[@var{ahat}, @var{bhat}, @var{aci}, @var{bci}] = unifit (@var{x})} also
## returns the 95% confidence intervals of the two estimates, one column per
## column of @var{x}, with the lower bound in the first row and the upper bound
## in the second.  @var{ahat} is the upper bound of @var{aci} and @var{bhat} the
## lower bound of @var{bci}, since no sample can fall outside the fitted range.
##
## @code{[@dots{}] = unifit (@var{x}, @var{alpha})} also returns the
## @qcode{100 * (1 - @var{alpha})} percent confidence intervals of the
## estimated parameters.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = unifit (@var{x}, @var{alpha}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## typically contains integer frequencies for the corresponding elements in
## @var{x}, but it can contain any non-integer non-negative values.  By default,
## or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.  This third
## argument is an Octave extension; MATLAB's @code{unifit} takes two inputs at
## most, and @var{freq} is accepted for a vector @var{x} only.
##
## Further information about the continuous uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Continuous_uniform_distribution}
##
## @seealso{unifcdf, unifinv, unifpdf, unifrnd, unifstat}
## @end deftypefn

function [ahat, bhat, aci, bci] = unifit (x, alpha, freq)

  ## Check input arguments
  if (nargin < 1)
    error ("unifit: function called with too few input arguments.");
  endif

  if (! isnumeric (x) || ! isreal (x))
    error ("unifit: X must be a vector or matrix of real values.");
  endif

  ## Check ALPHA
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("unifit: wrong value for ALPHA.");
  endif

  ## Check frequency vector
  if (nargin > 2 && ! isempty (freq))
    if (! isvector (x))
      error ("unifit: FREQ is supported for a vector X only.");
    elseif (! isequal (size (x), size (freq)))
      error ("unifit: X and FREQ vector mismatch.");
    elseif (any (freq < 0))
      error ("unifit: FREQ cannot have negative values.");
    endif
    ## Expand frequency vector
    if (! all (freq == 1))
      xf = [];
      for i = 1:numel (freq)
        xf = [xf, repmat(x(i), 1, freq(i))];
      endfor
      x = xf;
    endif
  endif

  if (isempty (x))
    ahat = [];
    bhat = [];
    aci = [];
    bci = [];
    return
  endif

  ## A vector is one sample whichever way it lies; a matrix is fitted by column
  if (isvector (x))
    x = x(:);
  endif

  ## Compute A and B estimates
  ahat = min (x, [], 1);
  bhat = max (x, [], 1);

  ## Compute the confidence intervals of A and B.  No sample can fall outside
  ## the fitted range, so AHAT bounds ACI from above and BHAT bounds BCI from
  ## below.
  if (nargout > 2)
    tmp = (bhat - ahat) ./ alpha .^ (1 ./ rows (x));
    aci = [bhat - tmp; ahat];
    bci = [bhat; ahat + tmp];
  endif

endfunction

%!demo
%! ## Sample 2 populations from different continuous uniform distributions
%! rng (42);
%! r1 = unifrnd (2, 5, 2000, 1);
%! r2 = unifrnd (3, 9, 2000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 0:0.5:10, 2);
%! h = findobj (gca, 'Type', 'patch');
%! set (h(1), 'facecolor', 'c');
%! set (h(2), 'facecolor', 'g');
%! hold on
%!
%! ## Estimate their probability of success
%! a_bA = unifit (r(:,1));
%! a_bB = unifit (r(:,2));
%!
%! ## Plot their estimated PDFs
%! x = [0:10];
%! y = unifpdf (x, a_bA(1), a_bA(2));
%! plot (x, y, '-pg');
%! y = unifpdf (x, a_bB(1), a_bB(2));
%! plot (x, y, '-sc');
%! xlim ([1, 10])
%! ylim ([0, 0.5])
%! legend ({'Normalized HIST of sample 1 with a=2 and b=5', ...
%!          'Normalized HIST of sample 2 with a=3 and b=9', ...
%!          sprintf("PDF for sample 1 with estimated a=%0.2f and b=%0.2f", ...
%!                  a_bA(1), a_bA(2)), ...
%!          sprintf("PDF for sample 2 with estimated a=%0.2f and b=%0.2f", ...
%!                  a_bB(1), a_bB(2))})
%! title ('Two population samples from different continuous uniform distributions')
%! hold off

## Test output
## Values below are R2024a's, measured 2026-08-17.
%!test
%! [ahat, bhat] = unifit (0:5);
%! assert_equal (ahat, 0);
%! assert_equal (bhat, 5);
%!test
%! [ahat, bhat, aci, bci] = unifit (0:5);
%! assert_equal (aci, [-3.237744862210329; 0], 1e-12);
%! assert_equal (bci, [5; 8.237744862210329], 1e-12);
%!test
%! [~, ~, aci, bci] = unifit (0:5, 0.10);
%! assert_equal (aci, [-2.338996338110347; 0], 1e-12);
%! assert_equal (bci, [5; 7.338996338110347], 1e-12);
%!test
%! ## a column vector is the same single sample as a row
%! [ahat, bhat, aci, bci] = unifit ((0:5)');
%! assert_equal ([ahat, bhat], [0, 5]);
%! assert_equal (aci, [-3.237744862210329; 0], 1e-12);
%!test
%! ## a matrix is fitted column by column
%! [ahat, bhat, aci, bci] = unifit ([0 10; 1 11; 2 12; 3 13; 4 14; 5 15]);
%! assert_equal (ahat, [0, 10]);
%! assert_equal (bhat, [5, 15]);
%! assert_equal (aci, [-3.237744862210329, 6.762255137789671; 0, 10], 1e-12);
%! assert_equal (bci, [5, 15; 8.237744862210329, 18.237744862210327], 1e-12);
%!test
%! ## negative data is ordinary for a uniform distribution
%! [ahat, bhat, aci, bci] = unifit ([-2, -1, 0, 1, 2]);
%! assert_equal ([ahat, bhat], [-2, 2]);
%! assert_equal (aci, [-5.282256812104322; -2], 1e-12);
%! assert_equal (bci, [2; 5.282256812104322], 1e-12);
%!test
%! ## a one-element sample has no width
%! [ahat, bhat, aci, bci] = unifit (5);
%! assert_equal ([ahat, bhat], [5, 5]);
%! assert_equal (aci, [5; 5]);
%! assert_equal (bci, [5; 5]);
%!test
%! ## empty data gives empty estimates
%! [ahat, bhat, aci, bci] = unifit ([]);
%! assert_equal (isempty (ahat), true);
%! assert_equal (isempty (bhat), true);
%! assert_equal (isempty (aci), true);
%!test
%! ## the endpoints ignore NaN
%! [ahat, bhat] = unifit ([0 1 NaN 3]);
%! assert_equal ([ahat, bhat], [0, 3]);
%!test
%! ## FREQ counts repeated observations
%! [a1, b1] = unifit ([1 1 2 3]);
%! [a2, b2] = unifit ([1 2 3], [], [2 1 1]);
%! assert_equal ([a1, b1], [a2, b2]);

## Test input validation
%!error<unifit: function called with too few input arguments.> unifit ()
%!error<unifit: X must be a vector or matrix of real values.> unifit ({1, 2})
%!error<unifit: X must be a vector or matrix of real values.> unifit (1+2i)
%!error<unifit: wrong value for ALPHA.> unifit (1, 0)
%!error<unifit: wrong value for ALPHA.> unifit (1, 1.2)
%!error<unifit: wrong value for ALPHA.> unifit (1, [0.02 0.05])
%!error<unifit: FREQ is supported for a vector X only.> ...
%! unifit ([1 2; 3 4], [], [1 1; 1 1])
%!error<unifit: X and FREQ vector mismatch.> ...
%! unifit ([1.5, 0.2], [], [0, 0, 0, 0, 0])
%!error<unifit: FREQ cannot have negative values.> ...
%! unifit ([1.5, 0.2], [], [1, -1])
