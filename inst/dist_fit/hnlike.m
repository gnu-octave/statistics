## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{nlogL} =} hnlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} hnlike (@var{params}, @var{x})
##
## Negative log-likelihood for the half-normal distribution.
##
## @code{@var{nlogL} = hnlike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the half-normal
## distribution with (1) location parameter @var{mu} and (2) scale parameter
## @var{sigma} given in the two-element vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = hnlike (@var{params}, @var{x})} returns
## the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{params} are their asymptotic variances.
##
## The half-normal CDF is only defined for @qcode{@var{x} >= @var{mu}}.
##
## Further information about the half-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Half-normal_distribution}
##
## @seealso{hncdf, hninv, hnpdf, hnrnd, hnfit}
## @end deftypefn

function [nlogL, acov] = hnlike (params, x)

  ## Check input arguments and add defaults
  if (nargin < 2)
    error ("hnlike: function called with too few input arguments.");
  endif
  if (numel (params) != 2)
    error ("hnlike: wrong parameters length.");
  endif

  ## Check X for being a vector
  if (isempty (x))
    phat = nan (1, 2, class (x));
    pci = nan (2, 2, class (x));
    return
  elseif (! isvector (x) || ! isreal (x))
    error ("hnlike: X must be a vector of real values.");
  endif

  ## Get MU and SIGMA parameters
  mu = params(1);
  sigma = params(2);

  ## Force X to column vector
  x = x(:);

  ## Return NaN for out of range parameters or data.
  sigma(sigma <= 0) = NaN;
  x(x < mu) = NaN;
  z = (x - mu) ./ sigma;

  ## Sum up the individual log-likelihood terms
  nlogL = -sum (-0.5 .* z .* z - log (sqrt (pi ./ 2) .* sigma));

  ## Compute asymptotic covariance (if requested)
  if (nargout == 2)
    nH = -sum (1 - 3 .* z .* z);
    avar =  (sigma .^ 2) ./ nH;
    acov = [0, 0; 0, avar];
  endif

endfunction

## Test output
%!test
%! x = 1:20;
%! paramhat = hnfit (x, 0);
%! [nlogL, acov] = hnlike (paramhat, x);
%! assert (nlogL, 64.179177404891300, 1e-14);

## Test input validation
%!error<hnlike: function called with too few input arguments.> ...
%! hnlike ([12, 15]);
%!error<hnlike: wrong parameters length.> hnlike ([12, 15, 3], [1:50]);
%!error<hnlike: wrong parameters length.> hnlike ([3], [1:50]);
%!error<hnlike: X must be a vector of real values.> ...
%! hnlike ([0, 3], ones (2));
%!error<hnlike: X must be a vector of real values.> ...
%! hnlike ([0, 3], [1, 2, 3, 4, 5+i]);
