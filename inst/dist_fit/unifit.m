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
## @deftypefn  {statistics} {@var{paramhat} =} unifit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} unifit (@var{x})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} unifit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{paramhat}, @var{paramci}] =} unifit (@var{x}, @var{alpha}, @var{freq})
##
## Estimate parameter and confidence intervals for the continuous uniform distribution.
##
## @code{@var{paramhat} = unifit (@var{x})} returns the maximum likelihood
## estimate (MLE) of the parameters @var{a} and @var{b} of the continuous
## uniform distribution given the data in @var{x}. @var{x} must be a vector.
##
## @code{[@var{paramhat}, @var{paramci}] = unifit (@var{x}, @var{alpha})} also
## returns the @qcode{100 * (1 - @var{alpha})} percent confidence intervals of
## the estimated parameter.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = unifit (@var{x}, @var{alpha}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## typically contains integer frequencies for the corresponding elements in
## @var{x}, but it can contain any non-integer non-negative values.  By default,
## or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the continuous uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Continuous_uniform_distribution}
##
## @seealso{unifcdf, unifinv, unifpdf, unifrnd, unifstat}
## @end deftypefn

function [paramhat, paramci] = unifit (x, alpha, freq)

  ## Check input arguments
  if (nargin < 1)
    error ("unifit: function called with too few input arguments.");
  endif

  ## Check data in X
  if (any (x < 0))
    error ("unifit: X cannot have negative values.");
  endif
  if (! isvector (x))
    error ("unifit: X must be a vector.");
  endif

  ## Check ALPHA
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("unifit: wrong value for ALPHA.");
  endif

  ## Check frequency vector
  if (nargin < 3 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("unifit: X and FREQ vector mismatch.");
  elseif (any (freq < 0))
    error ("unifit: FREQ cannot have negative values.");
  endif

  ## Expand frequency and censor vectors (if necessary)
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  ## Compute A and B estimates
  ahat = min (x);
  bhat = max (x);
  paramhat = [ahat, bhat];

  ## Compute confidence interval of A and B
  if (nargout > 1)
    tmp = (bhat - ahat) ./ alpha .^ (1 ./ numel (x));
    paramci = [bhat-tmp, ahat+tmp; ahat, bhat];
  endif

endfunction

%!demo
%! ## Sample 2 populations from different continuous uniform distributions
%! rand ("seed", 5);    # for reproducibility
%! r1 = unifrnd (2, 5, 2000, 1);
%! rand ("seed", 6);    # for reproducibility
%! r2 = unifrnd (3, 9, 2000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 0:0.5:10, 2);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! hold on
%!
%! ## Estimate their probability of success
%! a_bA = unifit (r(:,1));
%! a_bB = unifit (r(:,2));
%!
%! ## Plot their estimated PDFs
%! x = [0:10];
%! y = unifpdf (x, a_bA(1), a_bA(2));
%! plot (x, y, "-pg");
%! y = unifpdf (x, a_bB(1), a_bB(2));
%! plot (x, y, "-sc");
%! xlim ([1, 10])
%! ylim ([0, 0.5])
%! legend ({"Normalized HIST of sample 1 with a=2 and b=5", ...
%!          "Normalized HIST of sample 2 with a=3 and b=9", ...
%!          sprintf("PDF for sample 1 with estimated a=%0.2f and b=%0.2f", ...
%!                  a_bA(1), a_bA(2)), ...
%!          sprintf("PDF for sample 2 with estimated a=%0.2f and b=%0.2f", ...
%!                  a_bB(1), a_bB(2))})
%! title ("Two population samples from different continuous uniform distributions")
%! hold off

## Test output
%!test
%! x = 0:5;
%! [paramhat, paramci] = unifit (x);
%! assert (paramhat, [0, 5]);
%! assert (paramci, [-3.2377, 8.2377; 0, 5], 1e-4);
%!test
%! x = 0:5;
%! [paramhat, paramci] = unifit (x, [], [1 1 1 1 1 1]);
%! assert (paramhat, [0, 5]);
%! assert (paramci, [-3.2377, 8.2377; 0, 5], 1e-4);
%!assert (unifit ([1 1 2 3]), unifit ([1 2 3], [] ,[2 1 1]))

## Test input validation
%!error<unifit: function called with too few input arguments.> unifit ()
%!error<unifit: X cannot have negative values.> unifit (-1, [1 2 3 3])
%!error<unifit: wrong value for ALPHA.> unifit (1, 0)
%!error<unifit: wrong value for ALPHA.> unifit (1, 1.2)
%!error<unifit: wrong value for ALPHA.> unifit (1, [0.02 0.05])
%!error<unifit: X and FREQ vector mismatch.> ...
%! unifit ([1.5, 0.2], [], [0, 0, 0, 0, 0])
%!error<unifit: FREQ cannot have negative values.> ...
%! unifit ([1.5, 0.2], [], [1, -1])
%!error<unifit: X and FREQ vector mismatch.> ...
%! unifit ([1.5, 0.2], [], [1, 1, 1])
