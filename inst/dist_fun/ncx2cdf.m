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
## @deftypefn  {statistics} {@var{p} =} ncx2cdf (@var{x}, @var{df}, @var{delta})
## @deftypefnx {statistics} {@var{p} =} ncx2cdf (@var{x}, @var{df}, @var{delta}, @qcode{"upper"})
##
## Noncentral Chi-Square cumulative distribution function (CDF).
##
## @code{@var{p} = ncx2cdf (@var{x}, @var{df}, @var{delta})} returns the
## noncentral chi-square cdf with @var{df} degrees of freedom and noncentrality
## parameter @var{delta} at the values of @var{x}.
##
## The size of @var{p} is the common size of the input arguments. Scalar input
## arguments @var{x}, @var{df}, @var{delta} are regarded as constant matrices of
## the same size as the other inputs.
##
## @code{@var{p} = ncx2cdf (@var{x}, @var{df}, @var{delta}, "upper"} returns the
## upper tail probability of the noncentral chi-square distribution with
## @var{df} degrees of freedom and noncentrality parameter @var{delta} at the
## values in @var{x}.
##
## @seealso{ncx2inv, ncx2pdf, ncx2rnd, ncx2stat}
## @end deftypefn

function p = ncx2cdf (x, df, delta, uflag)

  ## Check for valid input arguments
  if (nargin <  3)
    error ("ncx2cdf: too few input arguments.");
  endif

  ## Check and fix size of input arguments
  [err, x, df, delta] = common_size (x, df, delta);
  if (err > 0)
    error ("ncx2cdf: input size mismatch.");
  endif

  ## Check for upper tail option
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("ncx2cdf: improper definition of upper tail option.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Initialize p
  if (isa (x, "single") || isa (df, "single") || isa (delta, "single"))
    p = zeros (size (x), "single");
    c_eps = eps ("single");
    c_min = realmin ("single");
  else
    p = zeros (size (x));
    c_eps = eps;
    c_min = realmin;
  endif

  ## Find NaNs in input arguments (if any) and propagate them to p
  is_nan = isnan (x) | isnan (df) | isnan (delta);
  p(is_nan) = NaN;
  if (uflag)
    p(x == Inf & ! is_nan) = 0;
    p(x <= 0 & ! is_nan) = 1;
  else
    p(x == Inf & ! is_nan) = 1;
  endif

  ## Make P = NaN for negative values of noncentrality parameter and DF
  p(delta < 0) = NaN;
  p(df < 0) = NaN;

  ## For DF == 0 at x == 0
  k = df == 0 & x == 0 & delta >= 0 & ! is_nan;
  if (uflag)
    p(k) = -expm1 (-delta(k) / 2);
  else
    p(k) = exp (-delta(k) / 2);
  endif

  ## Central chi2cdf
  k = df >= 0 & x > 0 & delta == 0 & isfinite (x) & ! is_nan;
  if (uflag)
    p(k) = chi2cdf (x(k), df(k), "upper");
  else
    p(k) = chi2cdf (x(k), df(k));
  endif

  ## Keep only valid samples
  td = find (df >= 0 & x > 0 & delta > 0 & isfinite (x) & ! is_nan);
  delta = delta(td) / 2;
  df = df(td) / 2;
  x = x(td) / 2;

  ## Compute Chernoff bounds
  e0 = log(c_min);
  e1 = log(c_eps/4);
  t = 1 - (df + sqrt (df .^ 2 + 4 * delta .* x)) ./ (2 * x);
  q = delta .* t ./ (1 - t) - df .* log(1 - t) - t .* x;
  peq0 = x < delta + df & q < e0;
  peq1 = x > delta + df & q < e1;
  if (uflag)
    p(td(peq0)) = 1;
  else
    p(td(peq1)) = 1;
  endif
  td(peq0 | peq1) = [];
  x(peq0 | peq1) = [];
  df(peq0 | peq1) = [];
  delta(peq0 | peq1) = [];

  ## Find index K of the maximal term in the summation series.
  ## K1 and K2 are lower and upper bounds for K, respectively.
  ## Indexing of terms in the summation series starts at 0.
  K1 = ceil ((sqrt ((df + x) .^ 2 + 4 * x .* delta) - (df + x)) / 2);
  K = zeros (size (x));
  k1above1 = find (K1 > 1);
  K2 = floor (delta(k1above1) .* gammaincratio (x(k1above1), K1(k1above1)));
  fixK2 = isnan(K2) | isinf(K2);
  K2(fixK2) = K1(k1above1(fixK2));
  K(k1above1) = K2;

  ## Find Poisson and Poisson*chi2cdf parts for the maximal terms in the
  ## summation series.
  if (uflag)
    k0 = (K==0 & df==0);
    K(k0) = 1;
  endif
  pois = poisspdf (K, delta);
  if (uflag)
    full = pois .* gammainc (x, df + K, "upper");
  else
    full = pois .* gammainc (x, df + K);
  endif

  ## Sum the series. First go downward from K and then go upward.
  ## The term for K is added afterwards - it is not included in either sum.
  sumK = zeros (size (x));

  ## Downward. poisspdf(k-1,delta)/poisspdf(k,delta) = k/delta
  poisterm = pois;
  fullterm = full;
  keep = K > 0 & fullterm > 0;
  k = K;
  while any(keep)
    poisterm(keep) = poisterm(keep) .* k(keep) ./ delta(keep);
    k(keep) = k(keep) - 1;
    if (uflag)
      fullterm(keep) = poisterm(keep) .* ...
                       gammainc (x(keep), df(keep) + k(keep), "upper");
    else
      fullterm(keep) = poisterm(keep) .* ...
                       gammainc (x(keep), df(keep) + k(keep));
    endif
    sumK(keep) = sumK(keep) + fullterm(keep);
    keep = keep & k > 0 & fullterm > eps(sumK);
  endwhile

  ## Upward. poisspdf(k+1,delta)/poisspdf(k,delta) = delta/(k+1)
  poisterm = pois;
  fullterm = full;
  keep = fullterm > 0;
  k = K;
  while any(keep)
    k(keep) = k(keep)+1;
    poisterm(keep) = poisterm(keep) .* delta(keep) ./ k(keep);
    if (uflag)
      fullterm(keep) = poisterm(keep) .* ...
                       gammainc (x(keep), df(keep) + k(keep), "upper");
    else
      fullterm(keep) = poisterm(keep) .* ...
                       gammainc (x(keep), df(keep) + k(keep));
    end
    sumK(keep) = sumK(keep) + fullterm(keep);
    keep = keep & fullterm > eps(sumK);
  endwhile

  ## Get probabilities
  p(td) = full + sumK;
  p(p > 1) = 1;

endfunction

## Ratio of incomplete gamma function values at S and S-1.
function r = gammaincratio (x, s)
  ## Initialize
  r = zeros (size (s));
  ## Finf small
  small = s < 2 | s <= x;
  ## For small S, use the ratio computed directly
  if (any (small(:)))
    r(small) = gammainc (x(small), s(small)) ./ ...
               gammainc (x(small), s(small) - 1);
  endif
  ## For large S, estimate numerator and denominator using 'scaledlower' option
  if (any (! small(:)))
    idx = find (! small);
    x = x(idx);
    s = s(idx);
    r(idx) = gammainc (x, s, "scaledlower") ./ ...
             gammainc (x, s - 1, "scaledlower") .* x ./ s;
  endif
endfunction

%!demo
%! ## Compare the noncentral chi-square cdf with DELTA = 2 to the
%! ## chi-square cdf with the same number of degrees of freedom (4):
%!
%! x = (0:0.1:10)';
%! ncx2 = ncx2cdf (x, 4, 2);
%! chi2 = chi2cdf (x, 4);
%! plot(x, ncx2, "b-", "LineWidth", 2);
%! hold on
%! plot (x, chi2, "g--", "LineWidth", 2);
%! legend ("ncx2", "chi2", "Location", "NorthWest");

## Input validation tests
%!error<ncx2cdf: too few input arguments.> p = ncx2cdf (2);
%!error<ncx2cdf: too few input arguments.> p = ncx2cdf (2, 4);
%!error<ncx2cdf: input size mismatch.> p = ncx2cdf (2,  [4, 3], [3, 4, 5]);
%!error<ncx2cdf: improper definition of upper tail option.> ...
%! p = ncx2cdf (2, 4, 2, "lower");

## Output validation tests
%!test
%! x = (-2:0.1:2)';
%! p = ncx2cdf (x, 10, 1);
%! assert (p([1:21]), zeros (21, 1), 3e-84);
%! assert (p(22), 1.521400636466575e-09, 1e-14);
%! assert (p(30), 6.665480510026046e-05, 1e-14);
%! assert (p(41), 0.002406447308399836, 1e-14);
%!test
%! p = ncx2cdf (12, 10, 3);
%! assert (p, 0.4845555602398649, 1e-14);
%!test
%! p = ncx2cdf (2, 3, 2);
%! assert (p, 0.2207330870741212, 1e-14);
%!test
%! p = ncx2cdf (2, 3, 2, "upper");
%! assert (p, 0.7792669129258789, 1e-14);
%!test
%! p = ncx2cdf ([3, 6], 3, 2, "upper");
%! assert (p, [0.6423318186400054, 0.3152299878943012], 1e-14);
