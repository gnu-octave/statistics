## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{p} = gamcdf (@var{x}, @var{a})
## @deftypefnx {statistics} @var{p} = gamcdf (@var{x}, @var{a}, @var{b})
## @deftypefnx {statistics} @var{p} = gamcdf (@dots{}, "upper")
## @deftypefnx {statistics} [@var{p}, @var{plo}, @var{pup}] = evcdf (@var{x}, @var{a}, @var{b}, @var{pcov})
## @deftypefnx {statistics} [@var{p}, @var{plo}, @var{pup}] = evcdf (@var{x}, @var{a}, @var{b}, @var{pcov}, @var{alpha})
## @deftypefnx {statistics} [@var{p}, @var{plo}, @var{pup}] = evcdf (@dots{}, "upper")
##
## Gamma cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the Gamma distribution with shape parameter @var{a} and
## scale @var{b}.  When called with only one parameter, then @var{b} defaults to
## 1.  The size of @var{p} is the common size of @var{x}, @var{a}, and @var{b}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## When called with three output arguments, @code{[@var{p}, @var{plo},
## @var{pup}]} it computes the confidence bounds for @var{p} when the input
## parameters @var{a} and @var{b} are estimates.  In such case, @var{pcov},
## a 2-by-2 matrix containing the covariance matrix of the estimated parameters,
## is necessary.  Optionally, @var{alpha} has a default value of 0.05, and
## specifies 100 * (1 - @var{alpha})% confidence bounds. @var{plo} and @var{pup}
## are arrays of the same size as @var{p} containing the lower and upper
## confidence bounds.
##
## @code{[@dots{}] = gamcdf (@dots{}, "upper")} computes the upper tail
## probability of the gamma distribution.
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
    da = dgammainc (z, a) .* dp;
    db = -exp (a .* log (z) - z - gammaln (a) - log (b)) .* dp;
    varLogitp = pcov(1,1) .* da .^ 2 + 2 .* pcov(1,2) .* da .* db + ...
                                            pcov(2,2) .* db .^ 2;
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

## Incomplete gamma function for first derivative
function dy = dgammainc(x, a)
  dy = NaN (size (x));
  amax = 2 ^ 20;
  k = find (a > amax);
  if (! isempty (k))
    x(k) = max(amax-1/3 + sqrt(amax./a(k)).*(x(k)-(a(k)-1/3)),0);
    a(k) = amax;
  endif
  ## Series expansion for lower incomplete gamma when x < a+1
  k = find (x < a + 1 & x != 0);
  if (! isempty (k))
    xk = x(k);
    ak = a(k);
    aplusn = ak;
    del = 1;
    ddel = 0;
    d2del = 0;
    sum = del;
    dsum = ddel;
    d2sum = d2del;
    while (norm (del, "inf") >= 100 * eps (norm (sum, "inf")))
      aplusn = aplusn + 1;
      del = del .* xk ./ aplusn;
      ddel = (ddel .* xk - del) ./ aplusn;
      d2del = (d2del .* xk - 2 .* ddel) ./ aplusn;
      sum = sum + del;
      dsum = dsum + ddel;
      d2sum = d2sum + d2del;
    endwhile
    fac = exp (-xk + ak .* log (xk) - gammaln (ak + 1));
    yk = fac .* sum;
    yk(xk > 0 & yk > 1) = 1;
    dlogfac = (log (xk) - psi (ak + 1));
    dfac = fac .* dlogfac;
    dyk = dfac .* sum + fac .* dsum;
    dy(k) = dyk;
  endif
  ## Continued fraction for upper incomplete gamma when x >= a+1
  k = find (x >= a + 1);
  if (! isempty (k))
    xk = x(k);
    ak = a(k);
    n = 0;
    a0 = 0;
    a1 = ak;
    b0 = 1;
    b1 = xk;
    da0 = 0; db0 = 0; da1 = 1; db1 = 0;
    d2a0 = 0; d2b0 = 0; d2a1 = 0; d2b1 = 0;
    g = ak ./ xk;
    dg = 1 ./ xk;
    d2g = 0;
    d2gold = 1;
    while (norm (d2g - d2gold, "inf") > 100 * eps (norm (d2g, "inf")))
      rescale = 1 ./ b1;
      n = n + 1;
      nminusa = n - ak;
      d2a0 = (d2a1 + d2a0 .* nminusa - 2 .* da0) .* rescale;
      d2b0 = (d2b1 + d2b0 .* nminusa - 2 .* db0) .* rescale;
      da0 = (da1 + da0 .* nminusa - a0) .* rescale;
      db0 = (db1 + db0 .* nminusa - b0) .* rescale;
      a0 = (a1 + a0 .* nminusa) .* rescale;
      b0 = 1 + (b0 .* nminusa) .* rescale;
      nrescale = n .* rescale;
      d2a1 = d2a0 .* xk + d2a1 .* nrescale;
      d2b1 = d2b0 .* xk + d2b1 .* nrescale;
      da1 = da0 .* xk + da1 .* nrescale;
      db1 = db0 .* xk + db1 .* nrescale;
      a1 = a0 .* xk + a1 .* nrescale;
      b1 = b0 .* xk + n;
      d2gold = d2g;
      g = a1 ./ b1;
      dg = (da1 - g.*db1) ./ b1;
      d2g = (d2a1 - dg.*db1 - g.*d2b1 - dg.*db1) ./ b1;
    endwhile
    fac = exp(-xk + ak.*log(xk) - gammaln(ak+1));
    yk = fac.*g;
    dlogfac = (log(xk) - psi(ak+1));
    dfac = fac .* dlogfac;
    dyk = dfac.*g + fac.*dg;
    dy(k) = -dyk;
  endif
  kx0 = find(x == 0);
  if (! isempty (kx0))
    dy(kx0) = 0;
  endif
  ka0 = find (a == 0);
  if (! isempty (ka0))
    ka0x0 = find(a == 0 & x == 0);
    dy(ka0x0) = -Inf;
  endif
endfunction

%!shared x, y, u
%! x = [-1, 0, 0.5, 1, 2, Inf];
%! y = [0, gammainc(x(2:end), 1)];
%! u = [0, NaN, NaN, 1, 0.1353352832366127, 0];
%!assert (gamcdf (x, ones (1,6), ones (1,6)), y, eps)
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
