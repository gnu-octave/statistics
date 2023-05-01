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
## @deftypefn  {statistics} {@var{nlogL} =} betalike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{avar}] =} betalike (@var{params}, @var{x})
##
## Negative log-likelihood for the Beta distribution.
##
## @code{@var{nlogL} = betalike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the Beta distribution
## with (1) shape parameter @math{α} and (2) shape parameter @math{β} given in
## the two-element vector @var{paramhat}.  Both parameters must be positive real
## numbers and the data in the range @math{[0,1]}.  Out of range parameters or
## data return @qcode{NaN}.
##
## @code{[@var{nlogL}, @var{avar}] = betalike (@var{params}, @var{x})} returns
## the inverse of Fisher's information matrix, @var{avar}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{params} are their asymptotic variances.
##
## The Beta distribution is defined on the open interval @math{(0,1)}.  However,
## @code{betafit} can also compute the unbounded beta likelihood function for
## data that include exact zeros or ones.  In such cases, zeros and ones are
## treated as if they were values that have been left-censored at
## @qcode{sqrt (realmin)} or right-censored at @qcode{1 - eps/2}, respectively.
##
## @seealso{betacdf, betainv, betapdf, betarnd, betafit, betastat}
## @end deftypefn

function [nlogL, avar] = betalike (params, x)

  ## Check input arguments and add defaults
  if (nargin < 2)
    error ("betalike: too few input arguments.");
  endif
  if (numel (params) != 2)
    error ("betalike: wrong parameters length.");
  endif

  ## Get α and β parameters
  a = params(1);
  b = params(2);

  ## Force X to column vector
  x = x(:);

  ## Return NaN for out of range parameters or data.
  a(a <= 0) = NaN;
  b(b <= 0) = NaN;
  xmin = min (x);
  xmax = max (x);
  x(! (0 <= x & x <= 1)) = NaN;

  ## Add tolerance for boundary conditions
  x_lo = sqrt (realmin (class (x)));
  x_hi = 1 - eps (class (x)) / 2;

  ## All values are strictly within the interval (0,1)
  if (all (x > x_lo) && all (x < x_hi))
    num0 = 0;
    num1 = 0;
    x_ct = x;
    numx = length (x_ct);

  ## Find boundary elements and process them separately
  else
    num0 = sum (x < x_lo);
    num1 = sum (x > x_hi);
    x_ct = x (x > x_lo & x < x_hi);
    numx = length (x_ct);
  endif

  ## Compute continuous log likelihood
  logx = log (x_ct);
  log1px = log1p (-x_ct);
  sumlogx = sum (logx);
  sumlog1px = sum (log1px);
  nlogL = numx * betaln (a, b) - (a - 1) * sumlogx - (b - 1) * sumlog1px;

  ## Include log likelihood for zeros
  if (num0 > 0)
    nlogL = nlogL - num0 * log (betainc (x_lo, a, b, "lower"));
  endif

  ## Include log likelihood for ones
  if (num1 > 0)
    nlogL = nlogL - num1 * log (betainc (x_hi, a, b, "upper"));
  endif

  ## Compute the asymptotic covariance
  if (nargout > 1)
    if (numel (x) < 2)
      error ("betalike: not enough data in X.");
    endif

    ## Compute the Jacobian of the likelihood for values (0,1)
    psiab = psi (a + b);
    psi_a = psi (a);
    psi_b = psi (b);
    J = [logx+psiab-psi_a, log1px+psiab-psi_b];

    ## Add terms into the Jacobian for the zero and one values.
    if (num0 > 0 || num1 > 0)
      dd = sqrt (eps (class (x)));
      aa = a + a * dd * [1, -1];
      bb = b + b * dd * [1, -1];
      ad = 2 * a *dd;
      bd = 2 * b *dd;
      if (num0 > 0)
        da = diff (log (betainc (x_lo, aa, b, "lower"))) / ad;
        db = diff (log (betainc (x_lo, a, bb, "lower"))) / bd;
        J = [J; repmat([da, db], num0, 1)];
      endif
      if num1 > 0
        da = diff (log (betainc (x_hi, aa, b, "upper"))) / ad;
        db = diff (log (betainc (x_hi, a, bb, "upper"))) / bd;
        J = [J; repmat([da, db], num1, 1)];
      endif
    endif

    ## Invert the inner product of the Jacobian to get the asymptotic covariance
    [~, R] = qr (J, 0);
    if (any (isnan (R(:))))
      avar = [NaN, NaN; NaN, NaN];
    else
      Rinv = R \ eye (2);
      avar = Rinv * Rinv';
    endif
  endif

endfunction

## Test output
%!test
%! x = 0.01:0.02:0.99;
%! [nlogL, avar] = betalike ([2.3, 1.2], x);
%! avar_out = [0.03691678, 0.02803056; 0.02803056, 0.03965629];
%! assert (nlogL, 17.873477715879040, 1e-14);
%! assert (avar, avar_out, 1e-7);
%!test
%! x = 0.01:0.02:0.99;
%! [nlogL, avar] = betalike ([1, 4], x);
%! avar_out = [0.02793282, 0.02717274; 0.02717274, 0.03993361];
%! assert (nlogL, 79.648061114839550, 1e-13);
%! assert (avar, avar_out, 1e-7);
%!test
%! x = 0.00:0.02:1;
%! [nlogL, avar] = betalike ([1, 4], x);
%! avar_out = [0.00000801564765, 0.00000131397245; ...
%!             0.00000131397245, 0.00070827639442];
%! assert (nlogL, 573.2008434477486, 1e-10);
%! assert (avar, avar_out, 1e-14);

## Test input validation
%!error<betalike: too few input arguments.> betalike ([12, 15]);
%!error betalike ([12, 15], 3, 5, 6, 8);
%!error<betalike: wrong parameters length.> betalike ([12, 15, 3], [1:50]);
