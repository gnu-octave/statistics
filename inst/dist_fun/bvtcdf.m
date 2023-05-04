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
## @deftypefn  {statistics} {@var{p} =} bvtcdf (@var{x}, @var{rho}, @var{df})
## @deftypefnx {statistics} {@var{p} =} bvtcdf (@var{x}, @var{rho}, @var{df}, @var{Tol})
##
## Bivariate Student's t cumulative distribution function (CDF).
##
## @code{@var{p} = bvtcdf (@var{x}, @var{rho}, @var{df})} will compute the
## bivariate student's t cumulative distribution function of @var{x}, which must
## be an @math{Nx2} matrix, given a correlation coefficient @var{rho}, which
## must be a scalar, and @var{df} degrees of freedom, which can be a scalar or a
## vector of positive numbers commensurate with @var{x}.
##
## @var{Tol} is the tolerance for numerical integration and by default
## @code{@var{Tol} = 1e-8}.
##
## @seealso{mvtcdf}
## @end deftypefn

function p = bvtcdf (x, rho, df, TolFun)

  narginchk (3,4);

  if (nargin < 4)
    TolFun = 1e-8;
  endif

  if (isa (x, "single") || isa (rho, "single") || isa (df, "single"))
    is_type = "single";
  else
    is_type = "double";
  endif

  if (abs (rho) < 1)
    largeNu = 1e4;
    general = ! (fix (df) == df & df < largeNu);

    if (isscalar (df))
      if general
        p = generalDF (x, rho, repmat (df, size (x, 1), 1), TolFun);
      else
        p = integerDF (x, rho, df);
      endif

    else
      p = zeros (size (x, 1), 1, is_type);
      ## For large of non-integer df
      if (any (general))
        p(general) = generalDF (x(general,:), rho, df(general), TolFun);
      endif

      ## For small integer df
      for i = find (! general(:)')
        p(i) = integerDF (x(i,:), rho, df(i));
      endfor
    endif

  elseif (rho == 1)
    p = tcdf (min (x, [], 2), df);
    p(any (isnan( x), 2)) = NaN;

  else
    p = tcdf (x(:,1), df) - tcdf (-x(:,2), df);
  endif

endfunction


## CDF for the bivariate t with integer degrees of freedom
function p = integerDF (x, rho, df)

  x1 = x(:,1);
  x2 = x(:,2);
  tau = 1 - rho .^ 2;
  x1rx2 = x1 - rho * x2;
  x2rx1 = x2 - rho * x1;
  sx1r2 = sign (x1rx2);
  sx2r1 = sign (x2rx1);
  dfx1s = df + x1 .^ 2;
  dfx2s = df + x2 .^ 2;
  tdfx1x2 = tau * dfx2s ./ x1rx2 .^ 2;
  tdfx2x1 = tau * dfx1s ./ x2rx1 .^ 2;
  x_tdf12 = 1 ./ (1 + tdfx1x2);
  x_tdf21 = 1 ./ (1 + tdfx2x1);
  y_tdf12 = 1 ./ (1 + 1 ./ tdfx1x2);
  y_tdf21 = 1 ./ (1 + 1 ./ tdfx2x1);
  sqrtDF = sqrt (df);
  halfDF = df/2;
  if (fix (halfDF) == halfDF)    # for even DF
    p1 = atan2 (sqrt (tau), -rho) ./ (2 * pi);
    c1 = x1 ./ (4 * sqrt (dfx1s));
    c2 = x2 ./ (4 * sqrt (dfx2s));
    beta12 = 2 *atan2 (sqrt (x_tdf12), sqrt (y_tdf12)) / pi;
    beta21 = 2 *atan2 (sqrt (x_tdf21), sqrt (y_tdf21)) / pi;
    p2 = (1 + sx1r2 .* beta12) .* c2 + (1 + sx2r1 .* beta21) .* c1;
    betaT12 = 2 * sqrt (x_tdf12 .* y_tdf12) / pi;
    betaT21 = 2 * sqrt (x_tdf21 .* y_tdf21) / pi;
    for j = 2:halfDF
      fact = df * (j - 1.5) / (j - 1);
      c2 = c2 .* fact ./ dfx2s;
      c1 = c1 .* fact ./ dfx1s;
      beta12 = beta12 + betaT12;
      beta21 = beta21 + betaT21;
      p2 = p2 + (1 + sx1r2 .* beta12) .* c2 + (1 + sx2r1 .* beta21) .* c1;
      fact = 2 * (j - 1) / (2 * (j - 1) + 1);
      betaT12 = fact * betaT12 .* y_tdf12;
      betaT21 = fact * betaT21 .* y_tdf21;
    endfor
  else                          # for odd DF
    x1x2p = x1.*x2;
    x1x2s = x1 + x2;
    t1 = sqrt (x1 .^ 2 - 2 * rho * x1x2p + x2 .^ 2 + tau * df);
    t2 = x1x2p + rho * df;
    t3 = x1x2p - df;
    p1 = atan2 (sqrtDF .* (-x1x2s .* t2 - t3 .* t1), ...
                t3 .* t2 - df .* x1x2s .* t1) ./ (2 * pi);
    p1 = p1 + (p1 < 0);
    p2 = 0;
    if (df > 1)
      c1 = sqrtDF .* x1 ./ (2 * pi .* dfx1s);
      c2 = sqrtDF .* x2 ./ (2 * pi .* dfx2s);
      betaT12 = sqrt (x_tdf12);
      betaT21 = sqrt (x_tdf21);
      beta12 = betaT12;
      beta21 = betaT21;
      p2 = (1 + sx1r2 .* beta12) .* c2 + (1 + sx2r1 .* beta21) .* c1;
      for j = 2:(halfDF - 0.5)
        fact = df * (j - 1) / (j - 0.5);
        c2 = fact * c2 ./ dfx2s;
        c1 = fact * c1 ./ dfx1s;
        fact = 1 - 0.5 / (j - 1);
        betaT12 = fact * betaT12 .* y_tdf12;
        betaT21 = fact * betaT21 .* y_tdf21;
        beta12 = beta12 + betaT12;
        beta21 = beta21 + betaT21;
        p2 = p2 + (1 + sx1r2 .* beta12) .* c2 + (1 + sx2r1 .* beta21) .* c1;
      endfor
    endif
  endif
  p = p1 + p2;

  ## Fix limit cases
  large = 1e10;
  p(x1 < -large | x2 < -large) = 0;
  p(x1 > large) = tcdf (x2(x1 > large), df);
  p(x2 > large) = tcdf (x1(x2 > large), df);

endfunction


## CDF for the bivariate t with arbitrary degrees of freedom.
function p = generalDF (x, rho, df, TolFun)

  n = size (x, 1);
    if (rho >= 0)
    p1 = tcdf (min (x, [], 2), df);
    p1(any (isnan (x), 2)) = NaN;
  else
    p1 = tcdf (x(:,1), df) - tcdf (-x(:,2), df);
    p1(p1 < 0) = 0;
  endif
  lo = asin (rho);
  hi = (sign (rho) + (rho == 0)) .* pi ./ 2;
  p2 = zeros (size (p1), class (rho));
  for i = 1:n
    b1 = x(i,1); b2 = x(i,2);
    v = df(i);
    if (isfinite (b1) && isfinite (b2))
      p2(i) = quadgk (@bvtIntegrand, lo, hi, "AbsTol", TolFun, "RelTol", 0);
    endif
  endfor
  p = p1 - p2 ./ (2 .* pi);

  function integrand = bvtIntegrand (theta)
    st = sin (theta);
    ct2 = cos (theta).^2;
    integrand = (1 ./ (1 + ((b1 * st - b2) .^ 2 ./ ct2 + b1 .^ 2) / v)) ...
                                           .^ (v / 2);
  endfunction

endfunction

## Test output
%!test
%! x = [1, 2];
%! rho = [1, 0.5; 0.5, 1];
%! df = 4;
%! assert (bvtcdf(x, rho(2), df), mvtcdf(x, rho, df), 1e-14);
%!test
%! x = [3, 2;2, 4;1, 5];
%! rho = [1, 0.5; 0.5, 1];
%! df = 4;
%! assert (bvtcdf(x, rho(2), df), mvtcdf(x, rho, df), 1e-14);

