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
##
## Negative log-likelihood for the generalized Pareto distribution.
##
## @code{@var{nlogL} = gplike (@var{params}, @var{x})} returns the negative
## log-likelihood of the data in @var{x} corresponding to the generalized Pareto
## distribution with (1) shape parameter @var{k} and (2) scale parameter
## @var{sigma} given in the two-element vector @var{params}.  @code{gplike}
## does not allow a location parameter and it must be assumed known, and
## subtracted from @var{x} before calling @code{gplike}.
##
## @code{[@var{nlogL}, @var{acov}] = gplike (@var{params}, @var{x})} returns
## the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{acov} are their asymptotic variances.   @var{acov}
## is based on the observed Fisher's information, not the expected information.
##
## When @qcode{@var{k} = 0} and @qcode{@var{mu} = 0}, the Generalized Pareto CDF
## is equivalent to the exponential distribution.  When @qcode{@var{k} > 0} and
## @code{@var{mu} = @var{k} / @var{k}} the Generalized Pareto is equivalent to
## the Pareto distribution.  The mean of the Generalized Pareto is not finite
## when @qcode{@var{k} >= 1} and the variance is not finite when
## @qcode{@var{k} >= 1/2}.  When @qcode{@var{k} >= 0}, the Generalized Pareto
## has positive density for @qcode{@var{x} > @var{mu}}, or, when
## @qcode{@var{mu} < 0},for
## @qcode{0 <= (@var{x} - @var{mu}) / @var{sigma} <= -1 / @var{k}}.
##
## Further information about the generalized Pareto distribution can be found at
## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
##
## @seealso{gpcdf, gpinv, gppdf, gprnd, gpfit, gpstat}
## @end deftypefn

function [nlogL, acov] = gplike (params, x)

  ## Check input arguments
  if (nargin < 2)
    error ("gplike: function called with too few input arguments.");
  endif
  if (! isvector (x))
    error ("gplike: X must be a vector.");
  endif
  if (numel (params) != 2)
    error ("gplike: PARAMS must be a two-element vector.");
  endif
  ## Get SHAPE and SCALE parameters
  shape = params(1);
  scale = params(2);
  ## Get sample size and scale x
  sz = numel (x);
  z = x ./ scale;

  ## For SHAPE > 0
  if (abs (shape) > eps)
    if (shape > 0 || max (z) < -1 / shape)
      sumLn = sum (log1p (shape .* z));
      nlogL = sz * log (scale) + (1 + 1 / shape) .* sumLn;
      if (nargout > 1)
        v = z ./ (1 + shape .* z);
        sumv = sum (v);
        sumvsq = sum (v .^ 2);
        nH11 = 2 * sumLn ./ shape ^ 3 - ...
               2 * sumv ./ shape ^ 2 - (1 + 1 / shape) .* sumvsq;
        nH12 = (-sumv + (shape + 1) .* sumvsq) ./ scale;
        nH22 = (-sz + 2 * (shape + 1) .* sumv - ...
                shape * (shape + 1) .* sumvsq) ./ scale ^ 2;
        acov = [nH22, -nH12; -nH12, nH11] / (nH11 * nH22 - nH12 * nH12);
      endif
    else
      ## The support of the GP when shape<0 is 0 < y < abs(scale/shape)
      nlogL = Inf;
      if (nargout > 1)
        acov = [NaN NaN; NaN NaN];
      endif
    endif
  else # For shape = 0
    ## Handle limit explicitly to prevent (1/0) * log(1) == Inf*0 == NaN.
    nlogL = sz*log(scale) + sum(z);
    if (nargout > 1)
      sumz = sum (z);
      sumzsq = sum (z .^ 2);
      sumzcb = sum (z .^ 3);
      nH11 = (2 / 3) * sumzcb - sumzsq;
      nH12 = (-sumz + sumzsq) ./ scale;
      nH22 = (-sz + 2 * sumz) ./ scale ^ 2;
      acov = [nH22, -nH12; -nH12, nH11] / (nH11 * nH22 - nH12 * nH12);
    endif
  endif
endfunction

## Test output
%!assert (gplike ([2, 3], 4), 3.047536764863501, 1e-14)
%!assert (gplike ([1, 2], 4), 2.890371757896165, 1e-14)
%!assert (gplike ([2, 3], [1:10]), 32.57864322725392, 1e-14)
%!assert (gplike ([1, 2], [1:10]), 31.65666282460443, 1e-14)
%!assert (gplike ([1, NaN], [1:10]), NaN)

## Test input validation
%!error<gplike: function called with too few input arguments.> gplike ()
%!error<gplike: function called with too few input arguments.> gplike (1)
%!error<gplike: X must be a vector.> gplike ([1, 2], [])
%!error<gplike: X must be a vector.> gplike ([1, 2], ones (2))
%!error<gplike: PARAMS must be a two-element vector.> gplike (2, [1:10])
