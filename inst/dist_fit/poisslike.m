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
## @deftypefn  {statistics} {@var{nlogL} =} poisslike (@var{lambda}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{avar}] =} poisslike (@var{lambda}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} poisslike (@var{lambda}, @var{x}, @var{freq})
##
## Negative log-likelihood for the Poisson distribution.
##
## @code{@var{nlogL} = poisslike (@var{lambda}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the Poisson
## distribution with rate parameter @var{lambda}.  @var{x} must be a vector of
## non-negative values.
##
## @code{[@var{nlogL}, @var{avar}] = poisslike (@var{lambda}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{avar}.  If the input
## rate parameter, @var{lambda}, is the maximum likelihood estimate, @var{avar}
## is its asymptotic variance.
##
## @code{[@dots{}] = poisslike (@var{lambda}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## typically contains integer frequencies for the corresponding elements in
## @var{x}, but it can contain any non-integer non-negative values.  By default,
## or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the Poisson distribution can be found at
## @url{https://en.wikipedia.org/wiki/Poisson_distribution}
##
## @seealso{poisscdf, poissinv, poisspdf, poissrnd, poissfit, poisstat}
## @end deftypefn

function [nlogL, avar] = poisslike (lambda, x, freq=[])

  ## Check input arguments
  if (nargin < 2)
    error ("poisslike: function called with too few input arguments.");
  endif

  if (! isscalar (lambda) || ! isnumeric (lambda) || lambda <= 0)
    error ("poisslike: LAMBDA must be a positive scalar.");
  endif

  if (! isvector (x) || any (x < 0))
    error ("poisslike: X must be a vector of non-negative values.");
  endif

  if (isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("poisslike: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("poisslike: FREQ must not contain negative values.");
  endif

  ## Compute negative log-likelihood and asymptotic covariance
  n = sum (freq);
  nlogL = - sum (freq .* log (poisspdf (x, lambda)));
  avar = lambda / n;

endfunction

## Test output
%!test
%! x = [1 3 2 4 5 4 3 4];
%! [nlogL, avar] = poisslike (3.25, x);
%! assert (nlogL, 13.9533, 1e-4)
%!test
%! x = [1 2 3 4 5];
%! f = [1 1 2 3 1];
%! [nlogL, avar] = poisslike (3.25, x, f);
%! assert (nlogL, 13.9533, 1e-4)

## Test input validation
%!error<poisslike: function called with too few input arguments.> poisslike (1)
%!error<poisslike: LAMBDA must be a positive scalar.> poisslike ([1 2 3], [1 2])
%!error<poisslike: X must be a vector of non-negative values.> ...
%! poisslike (3.25, ones (10, 2))
%!error<poisslike: X must be a vector of non-negative values.> ...
%! poisslike (3.25, [1 2 3 -4 5])
%!error<poisslike: X and FREQ vectors mismatch.> ...
%! poisslike (3.25, ones (10, 1), ones (8,1))
%!error<poisslike: FREQ must not contain negative values.> ...
%! poisslike (3.25, ones (1, 8), [1 1 1 1 1 1 1 -1])
