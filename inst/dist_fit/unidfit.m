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
## @deftypefn  {statistics} {@var{Nhat} =} unidfit (@var{x})
## @deftypefnx {statistics} {[@var{Nhat}, @var{Nci}] =} unidfit (@var{x})
## @deftypefnx {statistics} {[@var{Nhat}, @var{Nci}] =} unidfit (@var{x}, @var{alpha})
## @deftypefnx {statistics} {[@var{Nhat}, @var{Nci}] =} unidfit (@var{x}, @var{alpha}, @var{freq})
##
## Estimate parameter and confidence intervals for the discrete uniform distribution.
##
## @code{@var{Nhat} = unidfit (@var{x})} returns the maximum likelihood estimate
## (MLE) of the maximum observable value for the discrete uniform distribution.
## @var{x} must be a vector.
##
## @code{[@var{Nhat}, @var{Nci}] = unidfit (@var{x}, @var{alpha})} also
## returns the @qcode{100 * (1 - @var{alpha})} percent confidence intervals of
## the estimated parameter.  By default, the optional argument @var{alpha} is
## 0.05 corresponding to 95% confidence intervals.  Pass in @qcode{[]} for
## @var{alpha} to use the default values.
##
## @code{[@dots{}] = unidfit (@var{x}, @var{alpha}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## typically contains integer frequencies for the corresponding elements in
## @var{x}, but it can contain any non-integer non-negative values.  By default,
## or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the discrete uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Discrete_uniform_distribution}
##
## @seealso{unidcdf, unidinv, unidpdf, unidrnd, unidstat}
## @end deftypefn

function [Nhat, Nci] = unidfit (x, alpha, freq)

  ## Check input arguments
  if (nargin < 1)
    error ("unidfit: function called with too few input arguments.");
  endif

  ## Check data inX
  if (any (x < 0))
    error ("unidfit: X cannot have negative values.");
  endif
  if (! isvector (x))
    error ("unidfit: X must be a vector.");
  endif

  ## Check ALPHA
  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  elseif (! isscalar (alpha) || ! isreal (alpha) || alpha <= 0 || alpha >= 1)
    error ("unidfit: wrong value for ALPHA.");
  endif

  ## Check frequency vector
  if (nargin < 3 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("unidfit: X and FREQ vector mismatch.");
  elseif (any (freq < 0))
    error ("unidfit: FREQ cannot have negative values.");
  endif

  ## Expand frequency and censor vectors (if necessary)
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
    freq = ones (size (x));
  endif

  ## Compute N estimate
  Nhat = max (x);

  ## Compute confidence interval of N
  if (nargout > 1)
    Nci = [Nhat; ceil(Nhat ./ alpha .^ (1 ./ numel (x)))];
  endif

endfunction

%!demo
%! ## Sample 2 populations from different discrete uniform distibutions
%! rand ("seed", 1);    # for reproducibility
%! r1 = unidrnd (5, 1000, 1);
%! rand ("seed", 2);    # for reproducibility
%! r2 = unidrnd (9, 1000, 1);
%! r = [r1, r2];
%!
%! ## Plot them normalized and fix their colors
%! hist (r, 0:0.5:20.5, 1);
%! h = findobj (gca, "Type", "patch");
%! set (h(1), "facecolor", "c");
%! set (h(2), "facecolor", "g");
%! hold on
%!
%! ## Estimate their probability of success
%! NhatA = unidfit (r(:,1));
%! NhatB = unidfit (r(:,2));
%!
%! ## Plot their estimated PDFs
%! x = [0:10];
%! y = unidpdf (x, NhatA);
%! plot (x, y, "-pg");
%! y = unidpdf (x, NhatB);
%! plot (x, y, "-sc");
%! xlim ([0, 10])
%! ylim ([0, 0.4])
%! legend ({"Normalized HIST of sample 1 with N=5", ...
%!          "Normalized HIST of sample 2 with N=9", ...
%!          sprintf("PDF for sample 1 with estimated N=%0.2f", NhatA), ...
%!          sprintf("PDF for sample 2 with estimated N=%0.2f", NhatB)})
%! title ("Two population samples from different discrete uniform distibutions")
%! hold off

## Test output
%!test
%! x = 0:5;
%! [Nhat, Nci] = unidfit (x);
%! assert (Nhat, 5);
%! assert (Nci, [5; 9]);
%!test
%! x = 0:5;
%! [Nhat, Nci] = unidfit (x, [], [1 1 1 1 1 1]);
%! assert (Nhat, 5);
%! assert (Nci, [5; 9]);
%!assert (unidfit ([1 1 2 3]), unidfit ([1 2 3], [] ,[2 1 1]))

## Test input validation
%!error<unidfit: function called with too few input arguments.> unidfit ()
%!error<unidfit: X cannot have negative values.> unidfit (-1, [1 2 3 3])
%!error<unidfit: wrong value for ALPHA.> unidfit (1, 0)
%!error<unidfit: wrong value for ALPHA.> unidfit (1, 1.2)
%!error<unidfit: wrong value for ALPHA.> unidfit (1, [0.02 0.05])
%!error<unidfit: X and FREQ vector mismatch.> ...
%! unidfit ([1.5, 0.2], [], [0, 0, 0, 0, 0])
%!error<unidfit: X and FREQ vector mismatch.> ...
%! unidfit ([1.5, 0.2], [], [1, 1, 1])
%!error<unidfit: FREQ cannot have negative values.> ...
%! unidfit ([1.5, 0.2], [], [1, -1])
