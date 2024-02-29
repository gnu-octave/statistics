## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{nlogL} =} gplike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} gplike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} gplike (@var{params}, @var{x}, @var{freq})
##
## Negative log-likelihood for the generalized Pareto distribution.
##
## @code{@var{nlogL} = gplike (@var{params}, @var{x})} returns the negative
## log-likelihood of the data in @var{x} corresponding to the generalized Pareto
## distribution with (1) shape parameter @var{k}, (2) scale parameter
## @var{sigma}, and (3) location parameter @var{theta} given in the
## three-element vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = gplike (@var{params}, @var{x})} returns
## the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{acov} are their asymptotic variances.   @var{acov}
## is based on the observed Fisher's information, not the expected information.
##
## @code{[@dots{}] = gplike (@var{params}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## typically contains integer frequencies for the corresponding elements in
## @var{x}, but it can contain any non-integer non-negative values.  By default,
## or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## When @qcode{@var{k} = 0} and @qcode{@var{mu} = 0}, the Generalized Pareto CDF
## is equivalent to the exponential distribution.  When @qcode{@var{k} > 0} and
## @code{@var{mu} = @var{k} / @var{k}} the Generalized Pareto is equivalent to
## the Pareto distribution.  The mean of the Generalized Pareto is not finite
## when @qcode{@var{k} >= 1} and the variance is not finite when
## @qcode{@var{k} >= 1/2}.  When @qcode{@var{k} >= 0}, the Generalized Pareto
## has positive density for @qcode{@var{x} > @var{mu}}, or, when
## @qcode{@var{mu} < 0}, for
## @qcode{0 <= (@var{x} - @var{mu}) / @var{sigma} <= -1 / @var{k}}.
##
## Further information about the generalized Pareto distribution can be found at
## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
##
## @seealso{gpcdf, gpinv, gppdf, gprnd, gpfit, gpstat}
## @end deftypefn

function [nlogL, acov] = gplike (params, x, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("gplike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("gplike: X must be a vector.");
  endif

  if (numel (params) != 3)
    error ("gplike: PARAMS must be a three-element vector.");
  endif

  ## Parse FREQ argument or add default
  if (nargin < 3 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("gplike: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("gplike: FREQ must not contain negative values.");
  endif

  ## Expand frequency vector (if necessary)
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  ## Get K, SIGMA, and THETA parameters
  k = params(1);
  sigma = params(2);
  theta = params(3);
  ## Get sample size and sigma x
  sz = numel (x);
  z = (x - theta) ./ sigma;

  ## For K > 0
  if (abs (k) > eps)
    if (k > 0 || max (z) < -1 / k)
      sumLn = sum (log1p (k .* z));
      nlogL = sz * log (sigma) + (1 + 1 / k) .* sumLn;
      if (nargout > 1)
        z_kz = z ./ (1 + k .* z);
        sumv = sum (z_kz);
        sumvsq = sum (z_kz .^ 2);
        nH11 = 2 * sumLn ./ k ^ 3 - ...
               2 * sumv ./ k ^ 2 - (1 + 1 / k) .* sumvsq;
        nH12 = (-sumv + (k + 1) .* sumvsq) ./ sigma;
        nH22 = (-sz + 2 * (k + 1) .* sumv - ...
                k * (k + 1) .* sumvsq) ./ sigma ^ 2;
        acov = [nH22, -nH12; -nH12, nH11] / (nH11 * nH22 - nH12 * nH12);
        acov = [acov, [0; 0]; 0, 0, 0];
      endif
    else
      ## The support of the GP when k<0 is 0 < y < abs(sigma/k)
      nlogL = Inf;
      if (nargout > 1)
        acov = [NaN, NaN, 0; NaN, NaN, 0; 0, 0, 0];
      endif
    endif
  else # For k = 0
    nlogL = sz * log (sigma) + sum (z);
    if (nargout > 1)
      sumz = sum (z);
      sumzsq = sum (z .^ 2);
      sumzcb = sum (z .^ 3);
      nH11 = (2 / 3) * sumzcb - sumzsq;
      nH12 = (-sumz + sumzsq) ./ sigma;
      nH22 = (-sz + 2 * sumz) ./ sigma ^ 2;
      acov = [nH22, -nH12; -nH12, nH11] / (nH11 * nH22 - nH12 * nH12);
      acov = [acov, [0; 0]; 0, 0, 0];
    endif
  endif
endfunction

## Test output
%!test
%! k = 0.8937; sigma = 1.3230; theta = 1;
%! x = [2.2196, 11.9301, 4.3673, 1.0949, 6.5626, ...
%!      1.2109, 1.8576, 1.0039, 12.7917, 2.2590];
%! [nlogL, acov] = gplike ([k, sigma, theta], x);
%! assert (nlogL, 21.736, 1e-3);
%! assert (acov, [0.7249, -0.7351, 0; -0.7351, 1.3040, 0; 0, 0, 0], 1e-4);
%!assert (gplike ([2, 3, 0], 4), 3.047536764863501, 1e-14)
%!assert (gplike ([2, 3, 4], 8), 3.047536764863501, 1e-14)
%!assert (gplike ([1, 2, 0], 4), 2.890371757896165, 1e-14)
%!assert (gplike ([1, 2, 4], 8), 2.890371757896165, 1e-14)
%!assert (gplike ([2, 3, 0], [1:10]), 32.57864322725392, 1e-14)
%!assert (gplike ([2, 3, 2], [1:10] + 2), 32.57864322725392, 1e-14)
%!assert (gplike ([2, 3, 0], [1:10], ones (1,10)), 32.57864322725392, 1e-14)
%!assert (gplike ([1, 2, 0], [1:10]), 31.65666282460443, 1e-14)
%!assert (gplike ([1, 2, 3], [1:10] + 3), 31.65666282460443, 1e-14)
%!assert (gplike ([1, 2, 0], [1:10], ones (1,10)), 31.65666282460443, 1e-14)
%!assert (gplike ([1, NaN, 0], [1:10]), NaN)

## Test input validation
%!error<gplike: function called with too few input arguments.> gplike ()
%!error<gplike: function called with too few input arguments.> gplike (1)
%!error<gplike: X must be a vector.> gplike ([1, 2, 0], [])
%!error<gplike: X must be a vector.> gplike ([1, 2, 0], ones (2))
%!error<gplike: PARAMS must be a three-element vector.> gplike (2, [1:10])
%!error<gplike: PARAMS must be a three-element vector.> gplike ([2, 3], [1:10])
%!error<gplike: X and FREQ vectors mismatch.> ...
%! gplike ([1, 2, 0], ones (10, 1), ones (8,1))
%!error<gplike: FREQ must not contain negative values.> ...
%! gplike ([1, 2, 0], ones (1, 8), [1 1 1 1 1 1 1 -1])
