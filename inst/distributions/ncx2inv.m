## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{x} = ncx2inv (@var{p}, @var{df}, @var{delta})
##
## Inverse of the non-central chi-square cumulative distribution function (cdf).
##
## @code{@var{x} = ncx2inv (@var{p}, @var{df}, @var{delta})} returns the inverse
## of the noncentral chi-square distribution with @var{df} degrees of freedom
## and noncentrality parameter @var{delta}, at the probabilities of @var{p}.
##
## The size of @var{x} is the common size of @var{df} and @var{delta}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## @code{ncx2inv} uses Newton's method to converge to the solution.
##
## @seealso{ncx2cdf, ncx2pdf, ncx2rnd, ncx2stat}
## @end deftypefn

function x = ncx2inv (p, df, delta)

  ## Check for valid input arguments
  if (nargin <  3)
    error ("ncx2inv: too few input arguments.");
  endif

  ## Check and fix size of input arguments
  [err, p, df, delta] = common_size (p, df, delta);
  if (err > 0)
    error ("ncx2inv: input size mismatch.");
  endif

  ## Initialize x
  if (isa (p, "single") || isa (df, "single") || isa (delta, "single"))
    x = NaN (size (p), "single");
    crit = sqrt (eps ("single"));
  else
    x = NaN (size (p), "double");
    crit = sqrt (eps ("double"));
  endif

  ## For delta == 0, call chi2inv
  d0 = delta == 0;
  if (any (d0(:)))
    x(d0) = chi2inv (p(d0), df(d0));
    ## If delta == 0 for all entries, then return
    if (all (d0(:)))
      return;
    endif
  endif

  ## CDF with 0 d.d0. has a step at x=0.
  ## Check if CDF at x=0 exceeds the requested p.
  df0 = df==0 & delta > 0;
  if (any (df0(:)))
    p0 = zeros (size (p));
    p0(df0) = ncx2cdf (0, df(df0), delta(df0));
    df0 = df0 & p0 >= p;
    x(df0) = 0;
  endif

  valid = ! df0 & df > 0 & delta > 0;

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
  mn = df(k) + delta(k);
  variance = 2 * (df(k) + 2 * delta(k));
  temp = log (variance + mn .^ 2);
  mu = 2 * log (mn) - 0.5 * temp;
  sigma = -2 * log (mn) + temp;
  xk = exp (norminv (pk, mu, sigma));
  F = ncx2cdf (xk, df(k), delta(k));
  h = ones(size(xk),'like',xk);

  ## Start iteration with a break out loop
  while (count < count_limit)
    count = count + 1;
    h = (F - pk) ./ ncx2pdf (xk, df(k), delta(k));
    xnew = max (xk / 50, min (5 * xk, xk - h));
    newF = ncx2cdf (xnew, df(k), delta(k));
    while (true)
      worse = (abs (newF - pk) > abs (F - pk) * (1 + crit)) & ...
              (abs (xk - xnew) > crit * xk);
      if (! any (worse))
        break;
      endif
      xnew(worse) = 0.5 * (xnew(worse) + xk(worse));
      newF(worse) = ncx2cdf (xnew(worse), df(k(worse)), delta(k(worse)));
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

## Input validation tests
%!error<ncx2inv: too few input arguments.> p = ncx2inv ();
%!error<ncx2inv: too few input arguments.> p = ncx2inv (1);
%!error<ncx2inv: too few input arguments.> p = ncx2inv (1, 2);
%!error<ncx2inv: input size mismatch.> p = ncx2inv (1, [4, 3], [3, 4, 5]);

## Output validation tests
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
%! x = [0.1808   0.6456   1.1842   1.7650   2.3760   3.0105];
%! assert (ncx2inv ([0.05], [1, 2, 3, 4, 5, 6], 4), x, 1e-4);
%!test
%! x = [0.4887   0.6699   0.9012   1.1842   1.5164   1.8927];
%! assert (ncx2inv ([0.05], 3, [1, 2, 3, 4, 5, 6]), x, 1e-4);
%!test
%! x = [0.4887   0.6699   0.9012   1.1842   NaN   1.8927];
%! assert (ncx2inv ([0.05], 3, [1, 2, 3, 4, -1, 6]), x, 1e-4);
%!test
%! assert (ncx2inv (0.996, 5, 8), 35.51298862765576, 1e-14);
