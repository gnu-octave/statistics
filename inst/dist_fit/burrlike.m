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
## @deftypefn  {statistics} {@var{nlogL} =} burrlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} burrlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} burrlike (@var{params}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} burrlike (@var{params}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the Burr type XII distribution.
##
## @code{@var{nlogL} = burrlike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the Burr type XII
## distribution with (1) scale parameter @var{lambda}, (2) first shape parameter
## @var{c}, and (3) second shape parameter @var{k} given in the three-element
## vector @var{paramhat}.
##
## @code{[@var{nlogL}, @var{acov}] = burrlike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{acov} are their asymptotic variances.
##
## @code{[@dots{}] = burrlike (@var{params}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = burrlike (@var{params}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the Burr type XII distribution can be found at
## @url{https://en.wikipedia.org/wiki/Burr_distribution}
##
## @seealso{burrcdf, burrinv, burrpdf, burrrnd, burrfit, burrstat}
## @end deftypefn

function [nlogL, acov] = burrlike (params, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("burrlike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("burrlike: X must be a vector.");
  endif

  if (any (x < 0))
    error ("burrlike: X cannot have negative values.");
  endif

  if (length (params) != 3)
    error ("burrlike: PARAMS must be a three-element vector.");
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("burrlike: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("burrlike: X and FREQ vector mismatch.");
  endif

  ## Get parameters
  lambda = params(1);
  c = params(2);
  k = params(3);

  ## Precalculate some values
  xl = x ./ lambda;
  log_xl = log (xl);
  l1_xlc = log1p (xl .^ c);

  ## Avoid realmax overflow by approximation
  is_inf = isinf (l1_xlc);
  l1_xlc(is_inf) = c .* log (xl(is_inf));

  ## Force censoring vector into logical
  notc = ! censor;
  cens = ! notc;

  ## Compute neg-loglikelihood
  likeL = zeros (size (x));
  likeL(notc) = (c - 1) .* log_xl(notc) - (k + 1) .* l1_xlc(notc);
  likeL(cens) = -k .* l1_xlc(cens);
  nlogL = sum (freq (notc)) * log (lambda / k / c) - sum (freq .* likeL);

  ## Compute asymptotic covariance
  if (nargout > 1)
    ## Preallocate variables
    nH = zeros (3);
    d2V1 = zeros (size (x));
    d2V2 = d2V1;

    ## Precalculate some more values
    xlc       = xl .^ c;
    log_xl    = log (xl);
    xlc1      = (1 + xlc);
    xlc1sq    = xlc1.^2;
    invxlc1sq = (1 + 1./xlc).^2;

    ## Find realmax overflow
    is_inf = isinf (xlc);
    is_fin = ! is_inf;

    ## Compute each element of the negative Hessian
    d2V1(is_fin) = -((1 + c) ./ xlc(is_fin) + 1) ./ invxlc1sq(is_fin);
    d2V1(is_inf) = -1;
    d2V2(notc) = d2V1(notc) .* (k + 1) + 1;
    d2V2(cens) = d2V1(cens) .* k;
    nH(1,1) = c ./ lambda .^ 2 .* sum (freq .* d2V2);

    d2V1(is_fin) = xlc(is_fin) .* (c .* log_xl(is_fin) + xlc1(is_fin)) ...
                               ./ xlc1sq(is_fin);
    d2V1(is_inf) = 1;
    d2V2(notc) = (k + 1) .* d2V1(notc) - 1;
    d2V2(cens) = k .* d2V1(cens);
    nH(1,2) = sum (freq ./ lambda .* d2V2);
    nH(2,1) = nH(1,2);

    d2V1(is_fin) = xlc(is_fin) .* log_xl(is_fin) .^ 2 ./ xlc1sq(is_fin);
    d2V1(is_inf) = 0;
    d2V2(notc) = d2V1(notc) .* k + d2V1(notc);
    d2V2(cens) = d2V1(cens) .* k;
    nH(2,2) = -(sum (freq (notc))) ./ c .^ 2 - sum (freq .* d2V2);

    d2V1(is_fin) = xlc(is_fin) ./ xlc1(is_fin);
    d2V1(is_inf) = 1;
    nH(1,3) = (c ./ lambda) .* sum (freq .* d2V1);
    nH(3,1) = nH(1,3);
    nH(2,3) = -sum (freq .* d2V1 .* log_xl);
    nH(3,2) = nH(2,3);

    nH(3,3) = -(sum (freq (notc))) ./ k .^ 2;
    nH = -nH;
    ## Check negative Hessian is positive definite
    [R, p] = chol (nH);
    if (p > 0)
      warning ("burrlike: non positive definite Hessian matrix.");
      acov = NaN (3);
      return
    endif
    ## ACOV estimate is the negative inverse of the Hessian.
    Rinv = inv (R);
    acov = Rinv * Rinv;
  endif

endfunction


## Test output

## Test input validation
%!error<burrlike: function called with too few input arguments.> burrlike (3.25)
%!error<burrlike: X must be a vector.> burrlike ([1, 2, 3], ones (2))
%!error<burrlike: X cannot have negative values.> burrlike ([1, 2, 3], [-1, 3])
%!error<burrlike: PARAMS must be a three-element vector.> ...
%! burrlike ([1, 2], [1, 3, 5, 7])
%!error<burrlike: PARAMS must be a three-element vector.> ...
%! burrlike ([1, 2, 3, 4], [1, 3, 5, 7])
%!error<burrlike: X and CENSOR vector mismatch.> ...
%! burrlike ([1, 2, 3], [1:5], [0, 0, 0])
%!error<burrlike: X and FREQ vector mismatch.> ...
%! burrlike ([1, 2, 3], [1:5], [0, 0, 0, 0, 0], [1, 1, 1])
%!error<burrlike: X and FREQ vector mismatch.> ...
%! burrlike ([1, 2, 3], [1:5], [], [1, 1, 1])
