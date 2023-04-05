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
## @deftypefn  {statistics} {@var{nlogL} =} poisslike (@var{lambda}, @var{data})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} poisslike (@var{lambda}, @var{data})
## @deftypefnx {statistics} {[@dots{}] =} poisslike (@var{lambda}, @var{data}, @var{freq})
##
## Negative log-likelihood for the Poisson distribution.
##
## @seealso{poisscdf, poissinv, poisspdf, poissrnd, poissfit, poisstat}
## @end deftypefn

function [nlogL, acov] = poisslike (lambda, data, freq=[])

  ## Check input arguments
  if (nargin < 2)
    error ("poisslike: too few input arguments.");
  endif

  if (! isscalar (lambda) || ! isnumeric (lambda) || lambda <= 0)
    error ("poisslike: LAMBDA must be a positive scalar.");
  endif
  if (! isvector (data))
    error ("poisslike: DATA must be a vector.");
  endif

  if (isempty (freq))
    freq = ones (size (data));
  endif

  ## Compute negative log-likelihood and asymptotic covariance
  n = sum (freq);
  nlogL = - sum (freq .* log (poisspdf (data, lambda)));
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
%!error<poisslike: DATA must be a vector.> poisslike (3.25, ones (10, 2))
