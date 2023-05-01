## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{nlogL} =} evlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{avar}] =} evlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} evlike (@var{params}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} evlike (@var{params}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the extreme value distribution.
##
## @code{@var{nlogL} = evlike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the extreme value
## distribution (also known as the Gumbel or the type I generalized extreme
## value distribution) with (1) location parameter @var{mu} and (2) scale
## parameter @var{sigma} given in the two-element vector @var{paramhat}.
##
## @code{[@var{nlogL}, @var{acov}] = evlike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{acov} are their asymptotic variances.
##
## @code{[@dots{}] = evlike (@var{params}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = evlike (@var{params}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## The Gumbel distribution is used to model the distribution of the maximum (or
## the minimum) of a number of samples of various distributions.  This version
## is suitable for modeling minima.  For modeling maxima, use the alternative
## Gumbel likelihood function, @code{gumbellike}.
##
## Further information about the Gumbel distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gumbel_distribution}
##
## @seealso{evcdf, evinv, evpdf, evrnd, evfit, evstat, gumbellike}
## @end deftypefn

function [nlogL, avar] = evlike (params, x, censor, freq)

  ## Check input arguments and add defaults
  if (nargin < 2)
    error ("evlike: too few input arguments.");
  endif
  if (numel (params) != 2)
    error ("evlike: wrong parameters length.");
  endif
  if (! isvector (x))
    error ("evlike: X must be a vector.");
  endif
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("evlike: X and CENSOR vectors mismatch.");
  endif
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (isequal (size (x), size (freq)))
    nulls = find (freq == 0);
    if (numel (nulls) > 0)
      x(nulls) = [];
      censor(nulls) = [];
      freq(nulls) = [];
    endif
  else
    error ("evlike: X and FREQ vectors mismatch.");
  endif

  ## Get mu and sigma values
  mu = params(1);
  sigma = params(2);

  ## sigma must be positive, otherwise make it NaN
  if (sigma <= 0)
    sigma = NaN;
  endif

  ## Compute the individual log-likelihood terms.  Force a log(0)==-Inf for
  ## x from extreme right tail, instead of getting exp(Inf-Inf)==NaN.
  z = (x - mu) ./ sigma;
  expz = exp (z);
  L = (z - log (sigma)) .* (1 - censor) - expz;
  L(z == Inf) = -Inf;

  ## Neg-log-like is the sum of the individual contributions
  nlogL = -sum (freq .* L);

  ## Compute the negative hessian at the parameter values.
  ## Invert to get the observed information matrix.
  if (nargout == 2)
    unc = (1-censor);
    nH11 = sum(freq .* expz);
    nH12 = sum(freq .* ((z + 1) .* expz - unc));
    nH22 = sum(freq .* (z .* (z+2) .* expz - ((2 .* z + 1) .* unc)));
    avar =  (sigma .^ 2) * ...
            [nH22 -nH12; -nH12 nH11] / (nH11 * nH22 - nH12 * nH12);
  endif
endfunction

## Test output
%!test
%! x = 1:50;
%! [nlogL, avar] = evlike ([2.3, 1.2], x);
%! avar_out = [-1.2778e-13, 3.1859e-15; 3.1859e-15, -7.9430e-17];
%! assert (nlogL, 3.242264755689906e+17, 1e-14);
%! assert (avar, avar_out, 1e-3);
%!test
%! x = 1:50;
%! [nlogL, avar] = evlike ([2.3, 1.2], x * 0.5);
%! avar_out = [-7.6094e-05, 3.9819e-06; 3.9819e-06, -2.0836e-07];
%! assert (nlogL, 481898704.0472211, 1e-6);
%! assert (avar, avar_out, 1e-3);
%!test
%! x = 1:50;
%! [nlogL, avar] = evlike ([21, 15], x);
%! avar_out = [11.73913876598908, -5.9546128523121216; ...
%!             -5.954612852312121, 3.708060045170236];
%! assert (nlogL, 223.7612479380652, 1e-13);
%! assert (avar, avar_out, 1e-14);

## Test input validation
%!error<evlike: too few input arguments.> evlike ([12, 15]);
%!error<evlike: wrong parameters length.> evlike ([12, 15, 3], [1:50]);
%!error<evlike: X must be a vector.> evlike ([12, 3], ones (10, 2));
%!error<evlike: X and CENSOR> evlike ([12, 15], [1:50], [1, 2, 3]);
%!error<evlike: X and FREQ> evlike ([12, 15], [1:50], [], [1, 2, 3]);
