## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{nlogL} =} stbllike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} stbllike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} stbllike (@var{params}, @var{x}, @var{freq})
##
## Negative log-likelihood for the stable distribution.
##
## @code{@var{nlogL} = stbllike (@var{params}, @var{x})} returns the negative
## log-likelihood of the data in @var{x} corresponding to the stable
## distribution, in the Nolan @qcode{S0} parameterization, with (1) tail index
## @var{alpha}, (2) skewness @var{beta}, (3) scale @var{gam}, and (4) location
## @var{delta} given in the four-element vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = stbllike (@var{params}, @var{x})} also
## returns the inverse of the observed Fisher information matrix, @var{acov}.  If
## the input parameter values in @var{params} are the maximum likelihood
## estimates, the diagonal elements of @var{acov} are their asymptotic
## variances.  @var{acov} is based on the numerically evaluated Hessian of the
## negative log-likelihood, since the stable density has no closed form.
##
## @code{[@dots{}] = stbllike (@var{params}, @var{x}, @var{freq})} accepts a
## frequency vector, @var{freq}, of the same size as @var{x}.  @var{freq} must
## contain non-negative integer frequencies for the corresponding elements in
## @var{x}.  By default, or if left empty, @qcode{@var{freq} = ones (size
## (@var{x}))}.
##
## Further information about the stable distribution can be found at
## @url{https://en.wikipedia.org/wiki/Stable_distribution}
##
## @seealso{stblfit, stblpdf, stblcdf, stblinv, stblrnd}
## @end deftypefn

function [nlogL, acov] = stbllike (params, x, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("stbllike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("stbllike: X must be a vector.");
  endif

  if (numel (params) != 4)
    error ("stbllike: PARAMS must be a four-element vector.");
  endif

  if (nargin < 3 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("stbllike: X and FREQ vectors mismatch.");
  elseif (any (freq < 0))
    error ("stbllike: FREQ must not contain negative values.");
  elseif (any (fix (freq) != freq))
    error ("stbllike: FREQ must contain integer values.");
  endif

  ## Force column vectors
  x = x(:);
  freq = freq(:);

  ## Negative log-likelihood at the given parameters
  nlogL = stbl_nll (params, x, freq);

  ## Optionally return the asymptotic covariance from the observed Fisher
  ## information, computed as the inverse of the numeric Hessian
  if (nargout > 1)
    theta = params(:).';
    hstep = (eps ^ (1/4)) .* max (abs (theta), 1);
    H = num_hessian (@(t) stbl_nll (t, x, freq), theta, hstep);
    H = (H + H') / 2;
    [~, notpd] = chol (H);
    if (notpd != 0)
      warning (strcat ("stbllike: Fisher information matrix not positive", ...
                       " definite; returning NaN covariance."));
      acov = NaN (4);
    else
      acov = inv (H);
      acov = (acov + acov') / 2;
    endif
  endif

endfunction

## Negative log-likelihood of the stable density at parameter vector THETA.
## Out-of-range parameters return Inf so that the enclosing optimizers (e.g.
## the profile-likelihood search in proflik) stay within the feasible region.
function nll = stbl_nll (theta, x, freq)
  alpha = theta(1);
  beta = theta(2);
  gam = theta(3);
  if (alpha <= 0 || alpha > 2 || abs (beta) > 1 || gam <= 0)
    nll = Inf;
    return
  endif
  y = __stable_pdf__ (x, theta(1), theta(2), theta(3), theta(4));
  y(y <= 0) = realmin;
  nll = -sum (freq .* log (y));
endfunction

## Central finite-difference Hessian of NLLFUN at THETA with per-parameter step
function H = num_hessian (nllfun, theta, hstep)
  p = numel (theta);
  H = zeros (p);
  f0 = nllfun (theta);
  for i = 1:p
    ei = zeros (1, p);
    ei(i) = hstep(i);
    fpi = nllfun (theta + ei);
    fmi = nllfun (theta - ei);
    H(i,i) = (fpi - 2 * f0 + fmi) / (hstep(i) ^ 2);
    for j = (i + 1):p
      ej = zeros (1, p);
      ej(j) = hstep(j);
      fpp = nllfun (theta + ei + ej);
      fpm = nllfun (theta + ei - ej);
      fmp = nllfun (theta - ei + ej);
      fmm = nllfun (theta - ei - ej);
      H(i,j) = (fpp - fpm - fmp + fmm) / (4 * hstep(i) * hstep(j));
      H(j,i) = H(i,j);
    endfor
  endfor
endfunction

%!demo
%! ## Negative log-likelihood of a stable fit to simulated data
%! rand ("seed", 42);
%! x = stblrnd (1.5, 0.5, 1, 0, 150, 1);
%! phat = stblfit (x);
%! nlogL = stbllike (phat, x)

## Test output
%!test  # matches -sum (log (pdf)); stbllike uses a fast CF-inversion density
%!      # that differs from stblpdf by a few parts in 1e-5
%! x = [-2.3, -0.9, 0.1, 0.4, 1.2, 2.8, 5.1];
%! nlogL = stbllike ([1.5, 0.5, 1, 0], x);
%! assert (nlogL, - sum (log (stblpdf (x, 1.5, 0.5, 1, 0))), 1e-3);
%!test  # frequency weights replicate observations
%! x = [-1, 0.5, 2];
%! f = [2, 1, 3];
%! xr = [-1, -1, 0.5, 2, 2, 2];
%! assert (stbllike ([1.2, 0, 1, 0], x, f), stbllike ([1.2, 0, 1, 0], xr), 1e-6);
%!test  # acov is symmetric positive (co)variance at a sensible parameter
%! rand ("seed", 1);
%! x = stblrnd (1.6, 0, 1, 0, 150, 1);
%! [~, acov] = stbllike ([1.6, 0, 1, 0], x);
%! assert (issymmetric (acov, 1e-10));
%! assert (all (diag (acov) > 0));

## Test input validation
%!error <stbllike: function called with too few input arguments.> ...
%! stbllike ([1.5, 0, 1, 0])
%!error <stbllike: X must be a vector.> stbllike ([1.5, 0, 1, 0], ones (2, 2))
%!error <stbllike: PARAMS must be a four-element vector.> ...
%! stbllike ([1.5, 0, 1], [1, 2, 3])
%!error <stbllike: X and FREQ vectors mismatch.> ...
%! stbllike ([1.5, 0, 1, 0], [1, 2, 3], [1, 2])
%!error <stbllike: FREQ must not contain negative values.> ...
%! stbllike ([1.5, 0, 1, 0], [1, 2, 3], [1, -1, 2])
%!error <stbllike: FREQ must contain integer values.> ...
%! stbllike ([1.5, 0, 1, 0], [1, 2, 3], [1, 1.5, 2])
