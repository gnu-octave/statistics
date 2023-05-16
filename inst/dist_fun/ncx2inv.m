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
## @deftypefn  {statistics} {@var{x} =} ncx2inv (@var{p}, @var{df}, @var{lambda})
##
## Inverse of the noncentral chi-squared cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the noncentral chi-squared distribution with @var{df} degrees of freedom and
## noncentrality parameter @var{mu}.  The size of @var{x} is the common size of
## @var{p}, @var{df}, and @var{mu}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## @code{ncx2inv} uses Newton's method to converge to the solution.
##
## Further information about the noncentral chi-squared distribution can be
## found at @url{https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution}
##
## @seealso{ncx2cdf, ncx2pdf, ncx2rnd, ncx2stat, chi2inv}
## @end deftypefn

function x = ncx2inv (p, df, lambda)

  ## Check for valid number of input arguments
  if (nargin <  3)
    error ("ncx2inv: function called with too few input arguments.");
  endif

  ## Check for common size of P, DF, and LAMBDA
  [err, p, df, lambda] = common_size (p, df, lambda);
  if (err > 0)
    error ("ncx2inv: P, DF, and LAMBDA must be of common size or scalars.");
  endif

  ## Check for P, DF, and LAMBDA being reals
  if (iscomplex (p) || iscomplex (df) || iscomplex (lambda))
    error ("ncx2inv: P, DF, and LAMBDA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (df, "single") || isa (lambda, "single"))
    x = NaN (size (p), "single");
    crit = sqrt (eps ("single"));
  else
    x = NaN (size (p), "double");
    crit = sqrt (eps ("double"));
  endif

  ## For lambda == 0, call chi2inv
  d0 = lambda == 0;
  if (any (d0(:)))
    x(d0) = chi2inv (p(d0), df(d0));
    ## If lambda == 0 for all entries, then return
    if (all (d0(:)))
      return;
    endif
  endif

  ## CDF with 0 d.d0. has a step at x=0.
  ## Check if CDF at x=0 exceeds the requested p.
  df0 = df==0 & lambda > 0;
  if (any (df0(:)))
    p0 = zeros (size (p));
    p0(df0) = ncx2cdf (0, df(df0), lambda(df0));
    df0 = df0 & p0 >= p;
    x(df0) = 0;
  endif

  valid = ! df0 & df > 0 & lambda > 0;

  ## Force x = 0 for p == 0 and x = Inf for p == 1
  x(p == 0 & valid) = 0;
  x(p == 1 & valid) = Inf;
  ## Find valid samples within the range of 0 < p < 1
  k = find (p > 0 & p < 1 & valid);
  pk = p(k);

  ## Initialize counter
  count_limit = 100;
  count = 0;

  ## Supply a starting guess for the iteration.
  mn = df(k) + lambda(k);
  variance = 2 * (df(k) + 2 * lambda(k));
  temp = log (variance + mn .^ 2);
  mu = 2 * log (mn) - 0.5 * temp;
  sigma = -2 * log (mn) + temp;
  xk = exp (norminv (pk, mu, sigma));
  F = ncx2cdf (xk, df(k), lambda(k));
  h = ones(size(xk), class (xk));

  ## Start iteration with a break out loop
  while (count < count_limit)
    count = count + 1;
    h = (F - pk) ./ ncx2pdf (xk, df(k), lambda(k));
    xnew = max (xk / 50, min (5 * xk, xk - h));
    newF = ncx2cdf (xnew, df(k), lambda(k));
    while (true)
      worse = (abs (newF - pk) > abs (F - pk) * (1 + crit)) & ...
              (abs (xk - xnew) > crit * xk);
      if (! any (worse))
        break;
      endif
      xnew(worse) = 0.5 * (xnew(worse) + xk(worse));
      newF(worse) = ncx2cdf (xnew(worse), df(k(worse)), lambda(k(worse)));
    endwhile
    h = xk - xnew;
    x(k) = xnew;
    mask = (abs (h) > crit * abs (xk));
    if (! any (mask))
      break;
    endif
    k = k(mask);
    xk = xnew(mask);
    F = newF(mask);
    pk = pk(mask);
  endwhile

  if (count == count_limit)
    warning ("ncx2inv: did not converge.");
    fprintf ("ncx2inv: Last Step: %13.8f\n", h);
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the noncentral chi-squared distribution
%! p = 0.001:0.001:0.999;
%! x1 = ncx2inv (p, 2, 1);
%! x2 = ncx2inv (p, 2, 2);
%! x3 = ncx2inv (p, 2, 3);
%! x4 = ncx2inv (p, 4, 1);
%! x5 = ncx2inv (p, 4, 2);
%! x6 = ncx2inv (p, 4, 3);
%! plot (p, x1, "-r", p, x2, "-g", p, x3, "-k", ...
%!       p, x4, "-m", p, x5, "-c", p, x6, "-y")
%! grid on
%! ylim ([0, 10])
%! legend ({"df = 2, λ = 1", "df = 2, λ = 2", ...
%!          "df = 2, λ = 3", "df = 4, λ = 1", ...
%!          "df = 4, λ = 2", "df = 4, λ = 3"}, "location", "northwest")
%! title ("Noncentral chi-squared iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

%!demo
%! ## Compare the noncentral chi-squared CDF with LAMBDA = 2 to the
%! ## chi-squared CDF with the same number of degrees of freedom (4).
%!
%! p = 0.001:0.001:0.999;
%! x1 = ncx2inv (p, 4, 2);
%! x2 = chi2inv (p, 4);
%! plot (p, x1, "-", p, x2, "-");
%! grid on
%! ylim ([0, 10])
%! legend ({"Noncentral χ^2(4,2)", "χ^2(4)"}, "location", "northwest")
%! title ("Noncentral chi-squared vs chi-squared quantile functions")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!test
%! x = [0,0.3443,0.7226,1.1440,1.6220,2.1770,2.8436,3.6854,4.8447,6.7701,Inf];
%! assert (ncx2inv ([0:0.1:1], 2, 1), x, 1e-4);
%!test
%! x = [0,0.8295,1.6001,2.3708,3.1785,4.0598,5.0644,6.2765,7.8763,10.4199,Inf];
%! assert (ncx2inv ([0:0.1:1], 2, 3), x, 1e-4);
%!test
%! x = [0,0.5417,1.3483,2.1796,3.0516,4.0003,5.0777,6.3726,8.0748,10.7686,Inf];
%! assert (ncx2inv ([0:0.1:1], 1, 4), x, 1e-4);
%!test
%! x = [0.1808, 0.6456, 1.1842, 1.7650, 2.3760, 3.0105];
%! assert (ncx2inv (0.05, [1, 2, 3, 4, 5, 6], 4), x, 1e-4);
%!test
%! x = [0.4887, 0.6699, 0.9012, 1.1842, 1.5164, 1.8927];
%! assert (ncx2inv (0.05, 3, [1, 2, 3, 4, 5, 6]), x, 1e-4);
%!test
%! x = [1.3941, 1.6824, 2.0103, 2.3760, NaN, 3.2087];
%! assert (ncx2inv (0.05, 5, [1, 2, 3, 4, -1, 6]), x, 1e-4);
%!test
%! assert (ncx2inv (0.996, 5, 8), 35.51298862765576, 2e-13);

## Test input validation
%!error<ncx2inv: function called with too few input arguments.> ncx2inv ()
%!error<ncx2inv: function called with too few input arguments.> ncx2inv (1)
%!error<ncx2inv: function called with too few input arguments.> ncx2inv (1, 2)
%!error<ncx2inv: P, DF, and LAMBDA must be of common size or scalars.> ...
%! ncx2inv (ones (3), ones (2), ones (2))
%!error<ncx2inv: P, DF, and LAMBDA must be of common size or scalars.> ...
%! ncx2inv (ones (2), ones (3), ones (2))
%!error<ncx2inv: P, DF, and LAMBDA must be of common size or scalars.> ...
%! ncx2inv (ones (2), ones (2), ones (3))
%!error<ncx2inv: P, DF, and LAMBDA must not be complex.> ncx2inv (i, 2, 2)
%!error<ncx2inv: P, DF, and LAMBDA must not be complex.> ncx2inv (2, i, 2)
%!error<ncx2inv: P, DF, and LAMBDA must not be complex.> ncx2inv (2, 2, i)
