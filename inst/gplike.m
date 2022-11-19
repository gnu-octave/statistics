## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{nlogL} = gplike (@var{params}, @var{data})
## @deftypefnx {Function File} [@var{nlogL}, @var{avar}] = gplike (@var{params}, @var{data})
##
## Negative log-likelihood for the generalized Pareto distribution.
##
## @code{@var{nlogL} = gplike (@var{params}, @var{data})} returns the negative
## of the log-likelihood for the two-parameter generalized Pareto distribution,
## evaluated at @code{@var{params(1)} = SHAPE} and
## @code{@var{params(2)} = SCALE} given @var{data}.  @code{gplike} does not
## allow a LOCATION parameter.  @var{nlogL} is a scalar.
##
## @code{[@var{nlogL}, @var{avar}] = gplike (@var{params}, @var{data})} returns
## the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{acov} are their asymptotic variances.   @var{acov}
## is based on the observed Fisher's information, not the expected information.
##
## When @var{shape} = 0 and @var{location} = 0, the generalized Pareto
## distribution is equivalent to the exponential distribution.  When
## @code{@var{shape} > 0} and @code{@var{location} = @var{scale} / @var{shape}},
## the generalized Pareto distribution is equivalent to the Pareto distribution.
## The mean of the generalized Pareto distribution is not finite when
## @code{@var{shape} >= 1}, and the variance is not finite when
## @code{@var{shape} >= 1/2}.  When @code{@var{shape} >= 0}, the generalized
## Pareto distribution has positive density for @code{@var{x} > @var{location}},
## or, when @code{@var{shape} < 0}, for
## @code{0 <= (@var{x} -  @var{location}) / @var{scale} <= -1 / @var{shape}}.
##
## @seealso{gpcdf, gpinv, gppdf, gprnd, gpfit, gpstat}
## @end deftypefn

function [nlogL, acov] = gplike (params, data)

  ## Check input arguments
  if (nargin < 2)
    error ("gplike: too few input arguments.");
  endif
  if (! isvector (data))
    error ("gplike: DATA must be a vector.");
  endif
  if (numel (params) != 2)
    error ("gplike: PARAMS must be a two-element vector.");
  endif
  ## Get SHAPE and SCALE parameters
  shape = params(1);
  scale = params(2);
  ## Get sample size and scale data
  sz = numel (data);
  z = data ./ scale;

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
      % The support of the GP when shape<0 is 0 < y < abs(scale/shape)
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

## Test input validation
%!error<gplike: too few input arguments.> gplike ()
%!error<gplike: too few input arguments.> gplike (1)
%!error<gplike: DATA must be a vector.> gplike ([1, 2], [])
%!error<gplike: DATA must be a vector.> gplike ([1, 2], ones (2))
%!error<gplike: PARAMS must be a two-element vector.> gplike (2, [1:10])

## Test results against MATLAB output
%!assert (gplike ([2, 3], 4), 3.047536764863501, 1e-14)
%!assert (gplike ([1, 2], 4), 2.890371757896165, 1e-14)
%!assert (gplike ([2, 3], [1:10]), 32.57864322725392, 1e-14)
%!assert (gplike ([1, 2], [1:10]), 31.65666282460443, 1e-14)
%!assert (gplike ([1, NaN], [1:10]), NaN)
