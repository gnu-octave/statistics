## Copyright (C) 2012 Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{nlogL} =} gevlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} gevlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} gevlike (@var{params}, @var{x}, @var{freq})
##
## Negative log-likelihood for the generalized extreme value (GEV) distribution.
##
## @code{@var{nlogL} = gevlike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the GEV distribution
## with (1) shape parameter @var{k}, (2) scale parameter @var{sigma}, and (3)
## location parameter @var{mu} given in the three-element vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = gevlike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{acov} are their asymptotic variances.
##
## @code{[@dots{}] = gevlike (@var{params}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq}
## must contain non-negative integer frequencies for the corresponding elements
## in @var{x}.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## When @qcode{@var{k} < 0}, the GEV is the type III extreme value distribution.
## When @qcode{@var{k} > 0}, the GEV distribution is the type II, or Frechet,
## extreme value distribution.  If @var{W} has a Weibull distribution as
## computed by the @code{wblcdf} function, then @qcode{-@var{W}} has a type III
## extreme value distribution and @qcode{1/@var{W}} has a type II extreme value
## distribution.  In the limit as @var{k} approaches @qcode{0}, the GEV is the
## mirror image of the type I extreme value distribution as computed by the
## @code{evcdf} function.
##
## The mean of the GEV distribution is not finite when @qcode{@var{k} >= 1}, and
## the variance is not finite when @qcode{@var{k} >= 1/2}.  The GEV distribution
## has positive density only for values of @var{x} such that
## @qcode{@var{k} * (@var{x} - @var{mu}) / @var{sigma} > -1}.
##
## Further information about the generalized extreme value distribution can be
## found at
## @url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
##
## @seealso{gevcdf, gevinv, gevpdf, gevrnd, gevfit, gevstat}
## @end deftypefn

function [nlogL, acov] = gevlike (params, x, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("gevlike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("gevlike: X must be a vector.");
  endif

  if (length (params) != 3)
    error ("gevlike: PARAMS must be a three-element vector.");
  endif

  if (nargin < 3 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("gevlike: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("gevlike: FREQ must not contain negative values.");
  elseif (any (fix (freq) != freq))
    error ("gevlike: FREQ must contain integer values.");
  endif

  ## Expand frequency
  if (! all (freq == 1))
    xf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
    endfor
    x = xf;
  endif

  k = params(1);
  sigma = params(2);
  mu = params(3);

  ## Calculate negative log likelihood
  [nll, k_terms] = gevnll (x, k, sigma, mu);
  nlogL = sum (nll(:));

  ## Optionally calculate the first and second derivatives of the negative log
  ## likelihood with respect to parameters
  if (nargout > 1)
  	[Grad, kk_terms] = gevgrad (x, k, sigma, mu, k_terms);
    FIM = gevfim (x, k, sigma, mu, k_terms, kk_terms);
    acov = inv (FIM);
  endif

endfunction

## Internal function to calculate negative log likelihood for gevlike
function [nlogL, k_terms] = gevnll (x, k, sigma, mu)

  k_terms = [];
  a = (x - mu) ./ sigma;

  if (all (k == 0))
    nlogL = exp(-a) + a + log(sigma);
  else
    aa = k .* a;
    ## Use a series expansion to find the log likelihood more accurately
    ## when k is small
    if (min (abs (aa)) < 1E-3 && max (abs (aa)) < 0.5)
      k_terms = 1;
      sgn = 1;
      i = 0;
      while 1
        sgn = -sgn;
        i++;
        newterm = (sgn  / (i + 1)) * (aa .^ i);
        k_terms = k_terms + newterm;
        if (max (abs (newterm)) <= eps)
          break
        endif
      endwhile
      nlogL = exp (-a .* k_terms) + a .* (k + 1) .* k_terms + log (sigma);
    else
      b = 1 + aa;
      nlogL = b .^ (-1 ./ k) + (1 + 1 ./ k) .* log (b) + log (sigma);
      nlogL(b <= 0) = Inf;
    endif
  endif

endfunction

## Calculate the gradient of the negative log likelihood of x with respect
## to the parameters of the generalized extreme value distribution for gevlike
function [G, kk_terms] = gevgrad (x, k, sigma, mu, k_terms)

  kk_terms = [];
  G = ones(3, 1);
  ## Use the expressions for first derivatives that are the limits as k --> 0
  if (k == 0)
    a = (x - mu) ./ sigma;
    f = exp(-a) - 1;
    ## k
    g = a .* (1 + a .* f / 2);
    G(1) = sum(g(:));
    ## sigma
    g = (a .* f + 1) ./ sigma;
    G(2) = sum(g(:));
    ## mu
    g = f ./ sigma;
    G(3) = sum(g(:));
    return
  endif

  a = (x - mu) ./ sigma;
  b = 1 + k .* a;
  ## Negative log likelihood is locally infinite
  if (any (b <= 0))
    G(:) = 0;
    return
  endif
  ## k
  c = log(b);
  d = 1 ./ k + 1;
  ## Use a series expansion to find the gradient more accurately when k is small
  if (nargin > 4 && ! isempty (k_terms))
    aa = k .* a;
    f = exp (-a .* k_terms);
    kk_terms = 0.5;
    sgn = 1;
    i = 0;
    while 1
      sgn = -sgn;
      i++;
      newterm = (sgn * (i + 1) / (i + 2)) * (aa .^ i);
      kk_terms = kk_terms + newterm;
      if (max (abs (newterm)) <= eps)
        break
      endif
    endwhile
    g = a .* ((a .* kk_terms) .* (f - 1 - k) + k_terms);
  else
    g = (c ./ k - a ./ b) ./ (k .* b .^ (1/k)) - c ./ (k .^ 2) + a .* d ./ b;
  endif
  G(1) = sum(g(:));

  ## sigma
  ## Use a series expansion to find the gradient more accurately when k is small
  if nargin > 4 && ~isempty(k_terms)
    g = (1 - a .* (a .* k .* kk_terms - k_terms) .* (f - k - 1)) ./ sigma;
  else
    g = (a .* b .^ (-d) - (k + 1) .* a ./ b + 1) ./ sigma;
  endif
  G(2) = sum(g(:));

  ## mu
  ## Use a series expansion to find the gradient more accurately when k is small
  if (nargin > 4 && ! isempty (k_terms))
    g = - (a .* k .* kk_terms - k_terms) .* (f - k - 1) ./ sigma;
  else
    g = (b .^ (-d) - (k + 1) ./ b) ./ sigma;
  end
  G(3) = sum(g(:));

endfunction

## Internal function to calculate the Fisher information matrix for gevlike
function ACOV = gevfim (x, k, sigma, mu, k_terms, kk_terms)

  ACOV = ones(3);
  ## Use the expressions for second derivatives that are the limits as k --> 0
  if (k == 0)
    ## k, k
    a = (x - mu) ./ sigma;
    f = exp(-a);
    der = (a .^ 2) .* (a .* (a/4 - 2/3) .* f + 2/3 * a - 1);
    ACOV(1, 1) = sum(der(:));

    ## sigma, sigma
    der = (sigma .^ -2) .* (a .* ((a - 2) .* f + 2) - 1);
    ACOV(2, 2) = sum(der(:));

    ## mu, mu
    der = (sigma .^ -2) .* f;
    ACOV(3, 3) = sum(der(:));

    ## k, sigma
    der = (-a ./ sigma) .* (a .* (1 - a/2) .* f - a + 1);
    ACOV(1, 2) = ACOV(2, 1) = sum(der(:));

    ## k, mu
    der = (-1 ./ sigma) .* (a .* (1 - a/2) .* f - a + 1);
    ACOV(1, 3) = ACOV(3, 1) = sum(der(:));

    ## sigma, mu
    der = (1 + (a - 1) .* f) ./ (sigma .^ 2);
    ACOV(2, 3) = ACOV(3, 2) = sum(der(:));

    return
  endif

  ## General case
  z = 1 + k .* (x - mu) ./ sigma;
  ## k, k
  a = (x - mu) ./ sigma;
  b = k .* a + 1;
  c = log(b);
  d = 1 ./ k + 1;
  ## Use a series expansion to find the derivatives more accurately
  ## when k is small
  if (nargin > 5 && ! isempty (kk_terms))
    aa = k .* a;
    f = exp (-a .* k_terms);
    kkk_terms = 2/3;
    sgn = 1;
    i = 0;
    while 1
      sgn = -sgn; i++;
      newterm = (sgn * (i + 1) * (i + 2) / (i + 3)) * (aa .^ i);
      kkk_terms = kkk_terms + newterm;
      if (max (abs (newterm)) <= eps)
        break
      endif
    endwhile
    der = (a .^ 2) .* (a .* (a .* kk_terms .^ 2 - kkk_terms) .* ...
                       f + a .* (1 + k) .* kkk_terms - 2 * kk_terms);
  else
    der = ((((c ./ k.^2) - (a ./ (k .* b))) .^ 2) ./ (b .^ (1 ./ k))) + ...
    ((-2*c ./ k.^3) + (2*a ./ (k.^2 .* b)) + ((a ./ b) .^ 2 ./ k)) ./ ...
    (b .^ (1 ./ k)) + 2*c ./ k.^3 - (2*a ./ (k.^2 .* b)) - (d .* (a ./ b) .^ 2);
  endif
  der(z <= 0) = 0; # no probability mass in this region
  ACOV(1, 1) = sum (der(:));

  ## sigma, sigma
  ## Use a series expansion to find the derivatives more accurately
  ## when k is small
  if (nargin > 5 && ! isempty (kk_terms))
    der = ((-2*a .* k_terms + 4 * a .^ 2 .* k .* kk_terms - a .^ 3 .* ...
          (k .^ 2) .* kkk_terms) .* (f - k - 1) + f .* ((a .* ...
          (k_terms - a .* k .* kk_terms)) .^ 2) - 1) ./ (sigma .^ 2);
  else
    der = (sigma .^ -2) .* (-2 * a .* b .^ (-d) + d .* k .* a .^ 2 .* ...
          (b .^ (-d-1)) + 2 .* d .* k .* a ./ b - d .* (k .* a ./ b) .^ 2 - 1);
  end
  der(z <= 0) = 0; # no probability mass in this region
  ACOV(2, 2) = sum (der(:));

  ## mu, mu
  ## Use a series expansion to find the derivatives more accurately
  ## when k is small
  if (nargin > 5 && ! isempty (kk_terms))
      der = (f .* (a .* k .* kk_terms - k_terms) .^ 2 - a .* k .^ 2 .* ...
      kkk_terms .* (f - k - 1)) ./ (sigma .^ 2);
  else
    der = (d .* (sigma .^ -2)) .*  (k .* (b .^ (-d-1)) - (k ./ b) .^ 2);
  endif
  der(z <= 0) = 0; # no probability mass in this region
  ACOV(3, 3) = sum (der(:));

  ## k, mu
  ## Use a series expansion to find the derivatives more accurately
  ## when k is small
  if (nargin > 5 && ! isempty (kk_terms))
    der = 2 * a .* kk_terms .* (f - 1 - k) - a .^ 2 .* k_terms .* ...
                   kk_terms .* f + k_terms;
    der = -der ./ sigma;
  else
    der = ((b .^ (-d)) .* (c ./ k  - a ./ b) ./ k - a .* (b .^ (-d-1)) + ...
          ((1 ./ k) - d) ./ b + a .* k .* d ./ (b .^ 2)) ./ sigma;
  endif
  der(z <= 0) = 0; # no probability mass in this region
  ACOV(1, 3) = ACOV(3, 1) = sum (der(:));

  ## k, sigma
  der = a .* der;
  der(z <= 0) = 0; # no probability mass in this region
  ACOV(1, 2) = ACOV(2, 1) = sum (der(:));

  ## sigma, mu
  ## Use a series expansion to find the derivatives more accurately
  ## when k is small
  if (nargin > 5 && ! isempty (kk_terms))
    der = ((-k_terms + 3 * a .* k .* kk_terms - (a .* k) .^ 2 .* ...
            kkk_terms) .* (f - k - 1) + a .* (k_terms - a .* k .* ...
            kk_terms) .^ 2 .* f) ./ (sigma .^ 2);
  else
    der = (-(b .^ (-d)) + a .* k .* d .* (b .^ (-d-1)) + ...
          (d .* k ./ b) - a .* (k./b).^2 .* d) ./ (sigma .^ 2);
  end
  der(z <= 0) = 0; # no probability mass in this region
  ACOV(2, 3) = ACOV(3, 2) = sum (der(:));

endfunction

## Test output
%!test
%! x = 1;
%! k = 0.2;
%! sigma = 0.3;
%! mu = 0.5;
%! [L, C] = gevlike ([k sigma mu], x);
%! expected_L = 0.75942;
%! expected_C = [-0.12547 1.77884 1.06731; 1.77884 16.40761 8.48877; 1.06731 8.48877 0.27979];
%! assert (L, expected_L, 0.001);
%! assert (C, inv (expected_C), 0.001);
%!test
%! x = 1;
%! k = 0;
%! sigma = 0.3;
%! mu = 0.5;
%! [L, C] = gevlike ([k sigma mu], x);
%! expected_L = 0.65157;
%! expected_C = [0.090036 3.41229 2.047337; 3.412229 24.760027 12.510190; 2.047337 12.510190 2.098618];
%! assert (L, expected_L, 0.001);
%! assert (C, inv (expected_C), 0.001);
%!test
%! x = -5:-1;
%! k = -0.2;
%! sigma = 0.3;
%! mu = 0.5;
%! [L, C] = gevlike ([k sigma mu], x);
%! expected_L = 3786.4;
%! expected_C = [1.6802e-07, 4.6110e-06, 8.7297e-05; ...
%!               4.6110e-06, 7.5693e-06, 1.2034e-05; ...
%!               8.7297e-05, 1.2034e-05, -0.0019125];
%! assert (L, expected_L, -0.001);
%! assert (C, expected_C, -0.001);
%!test
%! x = -5:0;
%! k = -0.2;
%! sigma = 0.3;
%! mu = 0.5;
%! [L, C] = gevlike ([k sigma mu], x, [1, 1, 1, 1, 1, 0]);
%! expected_L = 3786.4;
%! expected_C = [1.6802e-07, 4.6110e-06, 8.7297e-05; ...
%!               4.6110e-06, 7.5693e-06, 1.2034e-05; ...
%!               8.7297e-05, 1.2034e-05, -0.0019125];
%! assert (L, expected_L, -0.001);
%! assert (C, expected_C, -0.001);

## Test input validation
%!error<gevlike: function called with too few input arguments.> gevlike (3.25)
%!error<gevlike: X must be a vector.> gevlike ([1, 2, 3], ones (2))
%!error<gevlike: PARAMS must be a three-element vector.> ...
%! gevlike ([1, 2], [1, 3, 5, 7])
%!error<gevlike: PARAMS must be a three-element vector.> ...
%! gevlike ([1, 2, 3, 4], [1, 3, 5, 7])
%!error<gevlike: X and FREQ vectors mismatch.> ...
%! gevlike ([5, 0.2, 1], ones (10, 1), ones (8,1))
%!error<gevlike: FREQ must not contain negative values.> ...
%! gevlike ([5, 0.2, 1], ones (1, 8), [1 1 1 1 1 1 1 -1])
%!error<gevlike: FREQ must contain integer values.> ...
%! gevlike ([5, 0.2, 1], ones (1, 8), [1 1 1 1 1 1 1 1.5])
