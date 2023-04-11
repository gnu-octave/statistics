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
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} poisslike (@var{lambda}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} poisslike (@var{lambda}, @var{x}, @var{freq})
##
## Negative log-likelihood for the Poisson distribution.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{lambda} is a scalar containing the rate parameter of the exponential
## distribution
## @item
## @var{x} is the vector of given values.
## @item
## @var{freq} is a vector of the same length as @var{x}, which typically
## contains integer frequencies for the corresponding elements in @var{x}.
## @var{freq} cannot contain negative values.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{nlogL} is the negative log-likelihood.
## @item
## @var{avar} is the inverse of the Fisher information matrix.
## (The Fisher information matrix is the second derivative of the negative log
## likelihood with respect to the parameter value.)
## @end itemize
##
## @seealso{poisscdf, poissinv, poisspdf, poissrnd, poissfit, poisstat}
## @end deftypefn

function [nlogL, acov] = poisslike (lambda, x, freq=[])

  ## Check input arguments
  if (nargin < 2)
    error ("poisslike: too few input arguments.");
  endif

  if (! isscalar (lambda) || ! isnumeric (lambda) || lambda <= 0)
    error ("poisslike: LAMBDA must be a positive scalar.");
  endif

  if (! isvector (x))
    error ("poisslike: X must be a vector.");
  endif

  if (isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("poisslike: FREQ must match X in size.");
  elseif (any (freq < 0))
    error ("poisslike: FREQ must not contain negative values.");
  endif

  ## Compute negative log-likelihood and asymptotic covariance
  n = sum (freq);
  nlogL = - sum (freq .* log (poisspdf (x, lambda)));
  acov = lambda / n;

endfunction

## test output
%!test
%! x = [1 3 2 4 5 4 3 4];
%! [nlogL, acov] = poisslike (3.25, x);
%! assert (nlogL, 13.9533, 1e-4)
%!test
%! x = [1 2 3 4 5];
%! f = [1 1 2 3 1];
%! [nlogL, acov] = poisslike (3.25, x, f);
%! assert (nlogL, 13.9533, 1e-4)

## test input validation
%!error<poisslike: too few input arguments.> poisslike (3.25)
%!error<poisslike: LAMBDA must be a positive scalar.> poisslike ([1 2 3], [1 2])
%!error<poisslike: X must be a vector.> poisslike (3.25, ones (10, 2))
%!error<poisslike: FREQ must match X in size.> ...
%! poisslike (3.25, ones (10, 1), ones (8,1))
%!error<poisslike: FREQ must not contain negative values.> ...
%! poisslike (3.25, ones (1, 8), [1 1 1 1 1 1 1 -1])
