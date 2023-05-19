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
## @deftypefn  {statistics} {@var{x} =} ncx2inv (@var{p}, @var{mu}, @var{k})
##
## Inverse of the von Mises cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the von Mises distribution with location parameter @var{mu} and concentration
## parameter @var{k} on the interval @math{[-pi,pi]}.  The size of @var{x} is
## the common size of @var{p}, @var{mu}, and @var{k}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## Note: the quantile of the von Mises distribution is not analytic.  Hence, it
## is approximated by a custom searching algorithm using its CDF until it
## converges up to a tolerance of @qcode{1e-5} or 100 iterations.  As a result,
## balancing between performance and accuracy, the accuracy is about
## @qcode{5e-5} for @qcode{@var{k} = 1} and it drops to @qcode{5e-5} as @var{k}
## increases.
##
## Further information about the von Mises distribution can be found at
## @url{https://en.wikipedia.org/wiki/Von_Mises_distribution}
##
## @seealso{vmcdf, vmpdf, nmrnd}
## @end deftypefn

function x = vminv (p, mu, k)

  ## Check for valid number of input arguments
  if (nargin <  3)
    error ("vminv: function called with too few input arguments.");
  endif

  ## Check for common size of P, MU, and K
  [err, p, mu, k] = common_size (p, mu, k);
  if (err > 0)
    error ("vminv: P, MU, and K must be of common size or scalars.");
  endif

  ## Check for P, MU, and K being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (k))
    error ("vminv: P, MU, and K must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") || isa (k, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p), "double");
  endif

  ## Process edge cases p=0, p=0.5, p=1
  p_0 = p < eps (class (x)) & k > 0 & isfinite (mu);
  x(p_0) = -pi + mu(p_0);
  p_5 = abs (p - 0.5) < eps (class (x)) & k > 0 & isfinite (mu);
  x(p_5) = mu(p_5);
  p_1 = 1 - p < eps (class (x)) & k > 0 & isfinite (mu);
  x(p_1) = pi + mu(p_1);

  ## Get remaining valid cases
  valc = p > 0 & p < 1 & ! p_5 & k > 0 & isfinite (mu);
  if (! any (valc))
    return
  endif
  p = p(valc);
  mu = mu(valc);
  k = k(valc);

  ## Complement cases of p<0.5 to 1-p and keep track to invert them at the end
  comp = p < 0.5;
  p(comp) = 1 - p(comp);

  ## Initialize counter and threshold
  crit = 1e-5;
  count_limit = 100;
  count = 0;

  ## Supply a starting guess for the iteration by linear interpolation with k=0
  x0 = 2 * pi .* p - pi;
  xz = zeros (size (p));
  xa = xz;
  ## Compute p0 and compare to target p
  p0 =  vmcdf (x0, 0, k);

  ## Solution is always 0 < x < x0. Search for x until p == p0 within threshold
  while (any (abs (p - p0) > crit) && count < count_limit)
    count = count + 1;
    xnew = xz + (abs (x0) - abs (xz))  .* 0.5;
    p0 = vmcdf (xnew, 0, k);
    ## Prepare for next step
    xdec = (p0 - p) > crit;
    xinc = (p - p0) > crit;
    if (any (xdec))
      x0(xdec) = xnew(xdec);
    endif
    if (any (xinc))
      xz(xinc) = xnew(xinc);
    endif
  endwhile

  ## Return the converged value(s).
  xnew(comp) = -xnew(comp);
  x(valc) = xnew + mu;

  if (count == count_limit)
    warning ("vminv: did not converge.");
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the von Mises distribution
%! p1 = [0,0.005,0.01:0.01:0.1,0.15,0.2:0.1:0.8,0.85,0.9:0.01:0.99,0.995,1];
%! x1 = vminv (p1, 0, 0.5);
%! x2 = vminv (p1, 0, 1);
%! x3 = vminv (p1, 0, 2);
%! x4 = vminv (p1, 0, 4);
%! plot (p1, x1, "-r", p1, x2, "-g", p1, x3, "-b", p1, x4, "-c")
%! grid on
%! ylim ([-pi, pi])
%! legend ({"μ = 0, k = 0.5", "μ = 0, k = 1", ...
%!          "μ = 0, k = 2", "μ = 0, k = 4"}, "location", "northwest")
%! title ("Von Mises iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared x, p0, p1
%! x = [-pi:pi/2:pi];
%! p0 = [0, 0.10975, 0.5, 0.89025, 1];
%! p1 = [0, 0.03752, 0.5, 0.99622, 1];
%!assert (vminv (p0, 0, 1), x, 5e-5)
%!assert (vminv (p0, zeros (1,5), ones (1,5)), x, 5e-5)
%!assert (vminv (p1, 0, [1 2 3 4 5]), x, [5e-5, 5e-4, 5e-5, 5e-4, 5e-5])

## Test input validation
%!error<vminv: function called with too few input arguments.> vminv ()
%!error<vminv: function called with too few input arguments.> vminv (1)
%!error<vminv: function called with too few input arguments.> vminv (1, 2)
%!error<vminv: P, MU, and K must be of common size or scalars.> ...
%! vminv (ones (3), ones (2), ones (2))
%!error<vminv: P, MU, and K must be of common size or scalars.> ...
%! vminv (ones (2), ones (3), ones (2))
%!error<vminv: P, MU, and K must be of common size or scalars.> ...
%! vminv (ones (2), ones (2), ones (3))
%!error<vminv: P, MU, and K must not be complex.> vminv (i, 2, 2)
%!error<vminv: P, MU, and K must not be complex.> vminv (2, i, 2)
%!error<vminv: P, MU, and K must not be complex.> vminv (2, 2, i)
