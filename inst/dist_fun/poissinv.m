## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 2014 Mike Giles
## Copyright (C) 2016 Lachlan Andrew
## Copyright (C) 1995-2017 Kurt Hornik
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
## @deftypefn  {statistics} {@var{x} =} poissinv (@var{p}, @var{lambda})
##
## Inverse of the Poisson cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the Poisson distribution with rate parameter @var{lambda}.  The size of
## @var{x} is the common size of @var{p} and @var{lambda}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## Further information about the Poisson distribution can be found at
## @url{https://en.wikipedia.org/wiki/Poisson_distribution}
##
## @seealso{poisscdf, poisspdf, poissrnd, poissfit, poisslike, poisstat}
## @end deftypefn

function x = poissinv (p, lambda)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("poissinv: function called with too few input arguments.");
  endif

  ## Check for common size of P and LAMBDA
  if (! isscalar (p) || ! isscalar (lambda))
    [retval, p, lambda] = common_size (p, lambda);
    if (retval > 0)
      error ("poissinv: P and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Check for P and LAMBDA being reals
  if (iscomplex (p) || iscomplex (lambda))
    error ("poissinv: P and LAMBDA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (lambda, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  ## Force NaN for out of range parameters or p-values
  k = (p < 0) | (p > 1) | isnan (p) | ! (lambda > 0);
  x(k) = NaN;

  k = (p == 1) & (lambda > 0);
  x(k) = Inf;

  k = (p > 0) & (p < 1) & (lambda > 0);
  if (any (k(:)))
    limit = 20;                       # After 'limit' iterations, use approx
    if (isscalar (lambda))
      cdf = [(cumsum (poisspdf (0:limit-1,lambda))), 2];
      y = p(:);                       # force to column
      r = bsxfun (@le, y(k), cdf);
      [~, x(k)] = max (r, [], 2);     # find first instance of p <= cdf
      x(k) -= 1;
    else
      kk = find (k);
      cdf = exp (-lambda(kk));
      for i = 1:limit
        m = find (cdf < p(kk));
        if (isempty (m))
          break;
        else
          x(kk(m)) += 1;
          cdf(m) += poisspdf (i, lambda(kk(m)));
        endif
      endfor
    endif

    ## Use Mike Giles's magic when x isn't < limit
    k &= (x == limit);
    if (any (k(:)))
      if (isscalar (lambda))
        lam = repmat (lambda, size (p));
      else
        lam = lambda;
      endif
      x(k) = analytic_approx (p(k), lam(k));
    endif
  endif

endfunction


## The following is based on Mike Giles's CUDA implementation,
## [http://people.maths.ox.ac.uk/gilesm/codes/poissinv/poissinv_cuda.h]
## which is copyright by the University of Oxford
## and is provided under the terms of the GNU GPLv3 license:
## http://www.gnu.org/licenses/gpl.html

function x = analytic_approx (p, lambda)

  s = norminv (p, 0, 1) ./ sqrt (lambda);
  k = (s > -0.6833501) & (s < 1.777993);
  ## use polynomial approximations in central region
  if (any (k))
    lam = lambda(k);
    if (isscalar (s))
      sk = s;
    else
      sk = s(k);
    endif

    ## polynomial approximation to f^{-1}(s) - 1
    rm =  2.82298751e-07;
    rm = -2.58136133e-06 + rm.*sk;
    rm =  1.02118025e-05 + rm.*sk;
    rm = -2.37996199e-05 + rm.*sk;
    rm =  4.05347462e-05 + rm.*sk;
    rm = -6.63730967e-05 + rm.*sk;
    rm =  0.000124762566 + rm.*sk;
    rm = -0.000256970731 + rm.*sk;
    rm =  0.000558953132 + rm.*sk;
    rm =  -0.00133129194 + rm.*sk;
    rm =   0.00370367937 + rm.*sk;
    rm =   -0.0138888706 + rm.*sk;
    rm =     0.166666667 + rm.*sk;
    rm =         sk + sk.*(rm.*sk);

    ## polynomial approximation to correction c0(r)

    t  =   1.86386867e-05;
    t  =  -0.000207319499 + t.*rm;
    t  =     0.0009689451 + t.*rm;
    t  =   -0.00247340054 + t.*rm;
    t  =    0.00379952985 + t.*rm;
    t  =   -0.00386717047 + t.*rm;
    t  =    0.00346960934 + t.*rm;
    t  =   -0.00414125511 + t.*rm;
    t  =    0.00586752093 + t.*rm;
    t  =   -0.00838583787 + t.*rm;
    t  =     0.0132793933 + t.*rm;
    t  =     -0.027775536 + t.*rm;
    t  =      0.333333333 + t.*rm;

    ##  O(1/lam) correction

    y  =   -0.00014585224;
    y  =    0.00146121529 + y.*rm;
    y  =   -0.00610328845 + y.*rm;
    y  =     0.0138117964 + y.*rm;
    y  =    -0.0186988746 + y.*rm;
    y  =     0.0168155118 + y.*rm;
    y  =     -0.013394797 + y.*rm;
    y  =     0.0135698573 + y.*rm;
    y  =    -0.0155377333 + y.*rm;
    y  =     0.0174065334 + y.*rm;
    y  =    -0.0198011178 + y.*rm;
    y ./= lam;

    x(k) = floor (lam + (y+t)+lam.*rm);
  endif

  k = ! k & (s > -sqrt (2));
  if (any (k))
    ## Newton iteration
    r = 1 + s(k);
    r2 = r + 1;
    while (any (abs (r - r2) > 1e-5))
      t = log (r);
      r2 = r;
      s2 = sqrt (2 * ((1-r) + r.*t));
      s2(r<1) *= -1;
      r = r2 - (s2 - s(k)) .* s2 ./ t;
      if (r < 0.1 * r2)
        r = 0.1 * r2;
      endif
    endwhile
    t = log (r);
    y = lambda(k) .* r + log (sqrt (2*r.*((1-r) + r.*t)) ./ abs (r-1)) ./ t;
    x(k) = floor (y - 0.0218 ./ (y + 0.065 * lambda(k)));
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the Poisson distribution
%! p = 0.001:0.001:0.999;
%! x1 = poissinv (p, 13);
%! x2 = poissinv (p, 4);
%! x3 = poissinv (p, 10);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r")
%! grid on
%! ylim ([0, 20])
%! legend ({"λ = 1", "λ = 4", "λ = 10"}, "location", "northwest")
%! title ("Poisson iCDF")
%! xlabel ("probability")
%! ylabel ("values in x (number of occurences)")

## Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (poissinv (p, ones (1,5)), [NaN 0 1 Inf NaN])
%!assert (poissinv (p, 1), [NaN 0 1 Inf NaN])
%!assert (poissinv (p, [1 0 NaN 1 1]), [NaN NaN NaN Inf NaN])
%!assert (poissinv ([p(1:2) NaN p(4:5)], 1), [NaN 0 NaN Inf NaN])

## Test class of input preserved
%!assert (poissinv ([p, NaN], 1), [NaN 0 1 Inf NaN NaN])
%!assert (poissinv (single ([p, NaN]), 1), single ([NaN 0 1 Inf NaN NaN]))
%!assert (poissinv ([p, NaN], single (1)), single ([NaN 0 1 Inf NaN NaN]))

## Test input validation
%!error<poissinv: function called with too few input arguments.> poissinv ()
%!error<poissinv: function called with too few input arguments.> poissinv (1)
%!error<poissinv: P and LAMBDA must be of common size or scalars.> ...
%! poissinv (ones (3), ones (2))
%!error<poissinv: P and LAMBDA must be of common size or scalars.> ...
%! poissinv (ones (2), ones (3))
%!error<poissinv: P and LAMBDA must not be complex.> poissinv (i, 2)
%!error<poissinv: P and LAMBDA must not be complex.> poissinv (2, i)
