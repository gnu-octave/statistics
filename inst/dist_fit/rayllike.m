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
## @deftypefn  {statistics} {@var{nlogL} =} rayllike (@var{sigma}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} rayllike (@var{sigma}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} rayllike (@var{sigma}, @var{x}, @var{freq})
##
## Negative log-likelihood for the Rayleigh distribution.
##
## @code{@var{nlogL} = rayllike (@var{sigma}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the Rayleigh
## distribution with rate parameter @var{sigma}.  @var{x} must be a vector of
## non-negative values.
##
## @code{[@var{nlogL}, @var{acov}] = rayllike (@var{sigma}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## rate parameter, @var{sigma}, is the maximum likelihood estimate, @var{acov}
## is its asymptotic variance.
##
## @code{[@dots{}] = rayllike (@var{sigma}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## typically contains integer frequencies for the corresponding elements in
## @var{x}, but it can contain any non-integer non-negative values.  By default,
## or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the Rayleigh distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rayleigh_distribution}
##
## @seealso{raylcdf, raylinv, raylpdf, raylrnd, raylfit, raylstat}
## @end deftypefn

function [nlogL, acov] = rayllike (sigma, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("rayllike: function called with too few input arguments.");
  endif

  if (! isscalar (sigma) || ! isnumeric (sigma) || sigma <= 0)
    error ("rayllike: SIGMA must be a positive scalar.");
  endif

  if (! isvector (x) || any (x < 0))
    error ("rayllike: X must be a vector of non-negative values.");
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("rayllike: X and CENSOR vectors mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("rayllike: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("rayllike: FREQ must not contain negative values.");
  endif

  ## Compute negative log-likelihood and asymptotic covariance
  zsq = (x / sigma) .^ 2;
  logfz = -2 * log (sigma) - zsq / 2 + log (x);
  dlogfz = (zsq - 2) / sigma;
  logS = - zsq / 2;
  dlogS = - zsq * 3 / sigma ^ 2;
  logfz(censor == 1) = logS(censor == 1);
  nlogL = - sum (freq .* logfz);
  if (nargout > 1)
    d2 = (2 - 3 * zsq) / sigma ^ 2;
    d2(censor == 1) = - 3 * zsq(censor == 1) / sigma ^ 2;
    acov = - sum (freq .* d2);
  endif

endfunction

## Test output
%!test
%! x = [1 3 2 4 5 4 3 4];
%! [nlogL, acov] = rayllike (3.25, x);
%! assert (nlogL, 14.7442, 1e-4)
%!test
%! x = [1 2 3 4 5];
%! f = [1 1 2 3 1];
%! [nlogL, acov] = rayllike (3.25, x, [], f);
%! assert (nlogL, 14.7442, 1e-4)
%!test
%! x = [1 2 3 4 5 6];
%! f = [1 1 2 3 1 0];
%! [nlogL, acov] = rayllike (3.25, x, [], f);
%! assert (nlogL, 14.7442, 1e-4)
%!test
%! x = [1 2 3 4 5 6];
%! c = [0 0 0 0 0 1];
%! f = [1 1 2 3 1 0];
%! [nlogL, acov] = rayllike (3.25, x, c, f);
%! assert (nlogL, 14.7442, 1e-4)

## Test input validation
%!error<rayllike: function called with too few input arguments.> rayllike (1)
%!error<rayllike: SIGMA must be a positive scalar.> rayllike ([1 2 3], [1 2])
%!error<rayllike: X must be a vector of non-negative values.> ...
%! rayllike (3.25, ones (10, 2))
%!error<rayllike: X must be a vector of non-negative values.> ...
%! rayllike (3.25, [1 2 3 -4 5])
%!error<rayllike: X and CENSOR vectors mismatch.> ...
%! rayllike (3.25, [1, 2, 3, 4, 5], [1 1 0]);
%!error<rayllike: X and CENSOR vectors mismatch.> ...
%! rayllike (3.25, [1, 2, 3, 4, 5], [1 1 0 1 1]');
%!error<rayllike: X and FREQ vectors mismatch.> ...
%! rayllike (3.25, [1, 2, 3, 4, 5], zeros (1,5), [1 1 0]);
%!error<rayllike: X and FREQ vectors mismatch.> ...
%! rayllike (3.25, [1, 2, 3, 4, 5], [], [1 1 0 1 1]');
%!error<rayllike: FREQ must not contain negative values.> ...
%! rayllike (3.25, ones (1, 8), [], [1 1 1 1 1 1 1 -1])
