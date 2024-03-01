## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{p} =} gamcdf (@var{x}, @var{a})
## @deftypefnx {statistics} {@var{p} =} gamcdf (@var{x}, @var{a}, @var{b})
## @deftypefnx {statistics} {@var{p} =} gamcdf (@dots{}, @qcode{"upper"})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} gamcdf (@var{x}, @var{a}, @var{b}, @var{pcov})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} gamcdf (@var{x}, @var{a}, @var{b}, @var{pcov}, @var{alpha})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} gamcdf (@dots{}, @qcode{"upper"})
##
## Gamma cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the Gamma distribution with shape parameter @var{a} and scale
## parameter @var{b}.  When called with only one parameter, then @var{b}
## defaults to 1.  The size of @var{p} is the common size of @var{x}, @var{a},
## and @var{b}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## When called with three output arguments, i.e. @qcode{[@var{p}, @var{plo},
## @var{pup}]}, @code{gamcdf} computes the confidence bounds for @var{p} when
## the input parameters @var{a} and @var{b} are estimates.  In such case,
## @var{pcov}, a @math{2x2} matrix containing the covariance matrix of the
## estimated parameters, is necessary.  Optionally, @var{alpha}, which has a
## default value of 0.05, specifies the @qcode{100 * (1 - @var{alpha})} percent
## confidence bounds.  @var{plo} and @var{pup} are arrays of the same size as
## @var{p} containing the lower and upper confidence bounds.
##
## @code{[@dots{}] = gamcdf (@dots{}, "upper")} computes the upper tail
## probability of the Gamma distribution with parameters @var{a} and
## @var{b}, at the values in @var{x}.
##
## OCTAVE/MATLAB use the alternative parameterization given by the pair
## @math{α, β}, i.e. shape @var{a} and scale @var{b}.  In Wikipedia, the two
## common parameterizations use the pairs @math{k, θ}, as shape and scale, and
## @math{α, β}, as shape and rate, respectively.  The parameter names @var{a}
## and @var{b} used here (for MATLAB compatibility) correspond to the parameter
## notation @math{k, θ} instead of the @math{α, β} as reported in Wikipedia.
##
## Further information about the Gamma distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gamma_distribution}
##
## @seealso{gaminv, gampdf, gamrnd, gamfit, gamlike, gamstat}
## @end deftypefn

function [varargout] = gamcdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2 || nargin > 6)
    error ("gamcdf: invalid number of input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 2 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 2 && ischar (varargin{end}) && ...
          ! strcmpi (varargin{end}, "upper"))
    error ("gamcdf: invalid argument for upper tail.");
  elseif (nargin > 2 && isempty (varargin{end}))
    uflag = false;
    varargin(end) = [];
  else
    uflag = false;
  endif

  ## Get extra arguments (if they exist) or add defaults
  a = varargin{1};
  if (numel (varargin) > 1)
    b = varargin{2};
  else
    b = 1;
  endif
  if (numel (varargin) > 2)
    pcov = varargin{3};
    ## Check for valid covariance matrix 2x2
    if (! isequal (size (pcov), [2, 2]))
      error ("gamcdf: invalid size of covariance matrix.");
    endif
  else
    ## Check that cov matrix is provided if 3 output arguments are requested
    if (nargout > 1)
      error ("gamcdf: covariance matrix is required for confidence bounds.");
    endif
    pcov = [];
  endif
  if (numel (varargin) > 3)
    alpha = varargin{4};
    ## Check for valid alpha value
    if (! isnumeric (alpha) || numel (alpha) !=1 || alpha <= 0 || alpha >= 1)
      error ("gamcdf: invalid value for alpha.");
    endif
  else
    alpha = 0.05;
  endif

  ## Check for common size of X, A, and B
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b))
    [err, x, a, b] = common_size (x, a, b);
    if (err > 0)
      error ("gamcdf: X, A, and B must be of common size or scalars.");
    endif
  endif

  ## Check for X, A, and B being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("gamcdf: X, A, and B must not be complex.");
  endif

  ## Prepare parameters so that gammainc returns NaN for out of range parameters
  a(a < 0) = NaN;
  b(b < 0) = NaN;

  ## Prepare data so that gammainc returns 0 for negative X
  x(x < 0) = 0;

  ## Compute gammainc
  z = x ./ b;
  if (uflag)
    p = gammainc (z, a, "upper");
    ## Fix NaNs to gammainc output when a == NaN
    p(isnan (a)) = NaN;
  else
    p = gammainc (z, a);
    ## Fix NaNs to gammainc output when a == NaN
    p(isnan (a)) = NaN;
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (a, "single") || isa (b, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  varargout{1} = cast (p, is_class);
  if (nargout > 1)
    plo = NaN (size (z), is_class);
    pup = NaN (size (z), is_class);
  endif

  ## Compute confidence bounds (if requested)
  if (nargout >= 2)
    ## Approximate the variance of p on the logit scale
    logitp = log (p ./ (1 - p));
    dp = 1 ./ (p .* (1 - p));
    dk = dgammainc (z, a) .* dp;
    dt = -exp (a .* log (z) - z - gammaln (a) - log (b)) .* dp;
    varLogitp = pcov(1,1) .* dk .^ 2 + 2 .* pcov(1,2) .* dk .* dt + ...
                                            pcov(2,2) .* dt .^ 2;
    if (any (varLogitp(:) < 0))
        error ("gamcdf: bad covariance matrix.");
    endif
    ## Use a normal approximation on the logit scale, then transform back to
    ## the original CDF scale
    halfwidth = -norminv (alpha / 2) * sqrt (varLogitp);
    explogitplo = exp (logitp - halfwidth);
    explogitpup = exp (logitp + halfwidth);
    plo = explogitplo ./ (1 + explogitplo);
    pup = explogitpup ./ (1 + explogitpup);
    varargout{2} = plo;
    varargout{3} = pup;
  endif

endfunction

## Compute 1st derivative of the incomplete Gamma function
function dy = dgammainc (x, a)

  ## Initialize return variables
  dy = nan (size (x));

  ## Use approximation for A > 2^20
  ulim = 2^20;
  is_lim = find (a > ulim);
  if (! isempty (is_lim))
    x(is_lim) = max (ulim - 1/3 + sqrt (ulim ./ a(is_lim)) .* ...
                     (x(is_lim) - (a(is_lim) - 1/3)), 0);
    a(is_lim) = ulim;
  endif

  ## For x < a+1
  is_lo = find (x < a + 1 & x != 0);
  if (! isempty (is_lo))
    x_lo = x(is_lo);
    k_lo = a(is_lo);
    k_1 = k_lo;
    step = 1;
    d1st = 0;
    stsum = step;
    d1sum = d1st;
    while norm (step, "inf") >= 100 * eps (norm (stsum, "inf"))
      k_1 += 1;
      step = step .* x_lo ./ k_1;
      d1st = (d1st .* x_lo - step) ./ k_1;
      stsum = stsum + step;
      d1sum = d1sum + d1st;
    endwhile
    fklo = exp (-x_lo + k_lo .* log (x_lo) - gammaln (k_lo + 1));
    ## Compute 1st derivative
    dlogfklo = (log (x_lo) - psi (k_lo + 1));
    d1fklo = fklo .* dlogfklo;
    d1y_lo = d1fklo .* stsum + fklo .* d1sum;
    dy(is_lo) = d1y_lo;
  endif

  ## For x >= a+1
  is_hi = find (x >= a+1);
  if (! isempty (is_hi))
    x_hi = x(is_hi);
    k_hi = a(is_hi);
    zc = 0;
    k0 = 0;
    k1 = k_hi;
    x0 = 1;
    x1 = x_hi;
    d1k0 = 0;
    d1k1 = 1;
    d1x0 = 0;
    d1x1 = 0;
    kx = k_hi ./ x_hi;
    d1kx = 1 ./ x_hi;
    d2kx = 0;
    start = 1;
    while norm (d2kx - start, "Inf") > 100 * eps (norm (d2kx, "Inf"))
      rescale = 1 ./ x1;
      zc += 1;
      n_k = zc - k_hi;
      d1k0 = (d1k1 + d1k0 .* n_k - k0) .* rescale;
      d1x0 = (d1x1 + d1x0 .* n_k - x0) .* rescale;
      k0 = (k1 + k0 .* n_k) .* rescale;
      x0 = 1 + (x0 .* n_k) .* rescale;
      nrescale = zc .* rescale;
      d1k1 = d1k0 .* x_hi + d1k1 .* nrescale;
      d1x1 = d1x0 .* x_hi + d1x1 .* nrescale;
      k1 = k0 .* x_hi + k1 .* nrescale;
      x1 = x0 .* x_hi + zc;
      start = d2kx;
      kx = k1 ./ x1;
      d1kx = (d1k1 - kx .* d1x1) ./ x1;
    endwhile
    fkhi = exp (-x_hi + k_hi .* log (x_hi) - gammaln (k_hi+1));
    ## Compute 1st derivative
    dlogfkhi = (log (x_hi) - psi (k_hi + 1));
    d1fkhi = fkhi .* dlogfkhi;
    d1y_hi = d1fkhi .* kx + fkhi .* d1kx;
    dy(is_hi) = -d1y_hi;
  endif

  ## Handle x == 0
  is_x0 = find (x == 0);
  if (! isempty (is_x0))
    dy(is_x0) = 0;
  endif

  ## Handle a == 0
  is_k0 = find (a == 0);
  if (! isempty (is_k0))
    is_k0x0 = find (a == 0 & x == 0);
    dy(is_k0x0) = -Inf;
  endif
endfunction

%!demo
%! ## Plot various CDFs from the Gamma distribution
%! x = 0:0.01:20;
%! p1 = gamcdf (x, 1, 2);
%! p2 = gamcdf (x, 2, 2);
%! p3 = gamcdf (x, 3, 2);
%! p4 = gamcdf (x, 5, 1);
%! p5 = gamcdf (x, 9, 0.5);
%! p6 = gamcdf (x, 7.5, 1);
%! p7 = gamcdf (x, 0.5, 1);
%! plot (x, p1, "-r", x, p2, "-g", x, p3, "-y", x, p4, "-m", ...
%!       x, p5, "-k", x, p6, "-b", x, p7, "-c")
%! grid on
%! legend ({"α = 1, β = 2", "α = 2, β = 2", "α = 3, β = 2", ...
%!          "α = 5, β = 1", "α = 9, β = 0.5", "α = 7.5, β = 1", ...
%!          "α = 0.5, β = 1"}, "location", "southeast")
%! title ("Gamma CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y, u
%! x = [-1, 0, 0.5, 1, 2, Inf];
%! y = [0, gammainc(x(2:end), 1)];
%! u = [0, NaN, NaN, 1, 0.1353352832366127, 0];
%!assert (gamcdf (x, ones (1,6), ones (1,6)), y, eps)
%!assert (gamcdf (x, ones (1,6), ones (1,6), []), y, eps)
%!assert (gamcdf (x, 1, ones (1,6)), y, eps)
%!assert (gamcdf (x, ones (1,6), 1), y, eps)
%!assert (gamcdf (x, [0, -Inf, NaN, Inf, 1, 1], 1), [1, NaN, NaN, 0, y(5:6)], eps)
%!assert (gamcdf (x, [0, -Inf, NaN, Inf, 1, 1], 1, "upper"), u, eps)
%!assert (gamcdf (x, 1, [0, -Inf, NaN, Inf, 1, 1]), [NaN, NaN, NaN, 0, y(5:6)], eps)
%!assert (gamcdf ([x(1:2), NaN, x(4:6)], 1, 1), [y(1:2), NaN, y(4:6)], eps)

## Test class of input preserved
%!assert (gamcdf ([x, NaN], 1, 1), [y, NaN])
%!assert (gamcdf (single ([x, NaN]), 1, 1), single ([y, NaN]), eps ("single"))
%!assert (gamcdf ([x, NaN], single (1), 1), single ([y, NaN]), eps ("single"))
%!assert (gamcdf ([x, NaN], 1, single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<gamcdf: invalid number of input arguments.> gamcdf ()
%!error<gamcdf: invalid number of input arguments.> gamcdf (1)
%!error<gamcdf: invalid number of input arguments.> gamcdf (1, 2, 3, 4, 5, 6, 7)
%!error<gamcdf: invalid argument for upper tail.> gamcdf (1, 2, 3, "uper")
%!error<gamcdf: invalid argument for upper tail.> gamcdf (1, 2, 3, 4, 5, "uper")
%!error<gamcdf: invalid size of covariance matrix.> gamcdf (2, 3, 4, [1, 2])
%!error<gamcdf: covariance matrix is required for confidence bounds.> ...
%! [p, plo, pup] = gamcdf (1, 2, 3)
%!error<gamcdf: covariance matrix is required for confidence bounds.> ...
%! [p, plo, pup] = gamcdf (1, 2, 3, "upper")
%!error<gamcdf: invalid value for alpha.> [p, plo, pup] = ...
%! gamcdf (1, 2, 3, [1, 0; 0, 1], 0)
%!error<gamcdf: invalid value for alpha.> [p, plo, pup] = ...
%! gamcdf (1, 2, 3, [1, 0; 0, 1], 1.22)
%!error<gamcdf: invalid value for alpha.> [p, plo, pup] = ...
%! gamcdf (1, 2, 3, [1, 0; 0, 1], "alpha", "upper")
%!error<gamcdf: X, A, and B must be of common size or scalars.> ...
%! gamcdf (ones (3), ones (2), ones (2))
%!error<gamcdf: X, A, and B must be of common size or scalars.> ...
%! gamcdf (ones (2), ones (3), ones (2))
%!error<gamcdf: X, A, and B must be of common size or scalars.> ...
%! gamcdf (ones (2), ones (2), ones (3))
%!error<gamcdf: X, A, and B must not be complex.> gamcdf (i, 2, 2)
%!error<gamcdf: X, A, and B must not be complex.> gamcdf (2, i, 2)
%!error<gamcdf: X, A, and B must not be complex.> gamcdf (2, 2, i)
%!error<gamcdf: bad covariance matrix.> ...
%! [p, plo, pup] = gamcdf (1, 2, 3, [1, 0; 0, -inf], 0.04)
