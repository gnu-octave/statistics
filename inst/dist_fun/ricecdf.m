## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} ricecdf (@var{x}, @var{nu}, @var{sigma})
## @deftypefnx {statistics} {@var{p} =} ricecdf (@var{x}, @var{nu}, @var{sigma}, @qcode{"upper"})
##
## Rician cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the Rician distribution with non-centrality (distance) parameter
## @var{nu} and scale parameter @var{sigma}.  The size of @var{p} is the common
## size of @var{x}, @var{nu}, and @var{sigma}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## @code{@var{p} = ricecdf (@var{x}, @var{nu}, @var{sigma}, "upper")} computes
## the upper tail probability of the Rician distribution with parameters
## @var{nu} and @var{sigma}, at the values in @var{x}.
##
## Further information about the Rician distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rice_distribution}
##
## @seealso{riceinv, ricepdf, ricernd, ricefit, ricelike, ricestat}
## @end deftypefn

function p = ricecdf (x, nu, sigma, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("ricecdf: function called with too few input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin == 4 && strcmpi (uflag, "upper"))
    uflag = true;
  elseif (nargin == 4  && ! strcmpi (uflag, "upper"))
    error ("ricecdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Check for common size of X and SIGMA
  if (! isscalar (x) || ! isscalar (nu) || ! isscalar (sigma))
    [retval, x, nu, sigma] = common_size (x, nu, sigma);
    if (retval > 0)
      error ("ricecdf: X, NU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X and SIGMA being reals
  if (iscomplex (x) || iscomplex (nu) || iscomplex (sigma))
    error ("ricecdf: X, NU, and SIGMA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (nu, "single") || isa (sigma, "single"));
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force 1 for upper flag and X <= 0
  k0 = nu >= 0 & sigma >= 0 & x < 0;
  if (uflag && any (k0(:)))
    p(k0) = 1;
  end

  ## Calculate Rayleigh CDF for valid parameter and data range
  k = nu >= 0 & sigma >= 0 & x >= 0;
  if (any (k(:)))
    if (uflag)
        p(k) = marcumQ1 (nu(k) ./ sigma(k), x(k) ./ sigma(k));
    else
        p(k) = 1 - marcumQ1 (nu(k) ./ sigma(k), x(k) ./ sigma(k));
    endif
  endif

  ## Continue argument check
  p(! (k0 | k)) = NaN;

endfunction

## Marcum's "Q" function of order 1
function Q = marcumQ1 (a, b)

  ## Prepare output matrix
  if (isa (a, "single") || isa (b, "single"))
   Q = NaN (size (b), "single");
  else
   Q = NaN (size (b));
  endif

  ## Force marginal cases
  Q(a != Inf & b == 0) = 1;
  Q(a != Inf & b == Inf) = 0;
  Q(a == Inf & b != Inf) = 1;
  z = isnan (Q) & a == 0 & b != Inf;
  if (any(z))
    Q(z) = exp ((-b(z) .^ 2) ./ 2);
  end

  ## Compute the remaining cases
  z = isnan (Q) & ! isnan (a) & ! isnan (b);
  if (any(z(:)))
    aa = (a(z) .^ 2) ./ 2;
    bb = (b(z) .^ 2) ./ 2;
    eA = exp (-aa);
    eB = bb .* exp (-bb);
    h = eA;
    d = eB .* h;
    s = d;
    j = (d > s.*eps(class(d)));
    k = 1;
    while (any (j))
      eA = aa .* eA ./ k;
      h = h + eA;
      eB = bb .* eB ./ (k + 1);
      d = eB .* h;
      s(j) = s (j) + d(j);
      j = (d > s .* eps (class (d)));
      k = k + 1;
    endwhile
    Q(z) = 1 - s;
  endif
endfunction

%!demo
%! ## Plot various CDFs from the Rician distribution
%! x = 0:0.01:10;
%! p1 = ricecdf (x, 0, 1);
%! p2 = ricecdf (x, 0.5, 1);
%! p3 = ricecdf (x, 1, 1);
%! p4 = ricecdf (x, 2, 1);
%! p5 = ricecdf (x, 4, 1);
%! plot (x, p1, "-b", x, p2, "g", x, p3, "-r", x, p4, "-m", x, p5, "-k")
%! grid on
%! ylim ([0, 1])
%! xlim ([0, 8])
%! legend ({"ν = 0, σ = 1", "ν = 0.5, σ = 1", "ν = 1, σ = 1", ...
%!          "ν = 2, σ = 1", "ν = 4, σ = 1"}, "location", "southeast")
%! title ("Rician CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

%!demo
%! ## Plot various CDFs from the Rician distribution
%! x = 0:0.01:10;
%! p1 = ricecdf (x, 0, 0.5);
%! p2 = ricecdf (x, 0, 2);
%! p3 = ricecdf (x, 0, 3);
%! p4 = ricecdf (x, 2, 2);
%! p5 = ricecdf (x, 4, 2);
%! plot (x, p1, "-b", x, p2, "g", x, p3, "-r", x, p4, "-m", x, p5, "-k")
%! grid on
%! ylim ([0, 1])
%! xlim ([0, 8])
%! legend ({"ν = 0, σ = 0.5", "ν = 0, σ = 2", "ν = 0, σ = 3", ...
%!          "ν = 2, σ = 2", "ν = 4, σ = 2"}, "location", "southeast")
%! title ("Rician CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!test
%! x = 0:0.5:2.5;
%! nu = 1:6;
%! p = ricecdf (x, nu, 1);
%! expected_p = [0.0000, 0.0179, 0.0108, 0.0034, 0.0008, 0.0001];
%! assert (p, expected_p, 0.001);
%!test
%! x = 0:0.5:2.5;
%! sigma = 1:6;
%! p = ricecdf (x, 1, sigma);
%! expected_p = [0.0000, 0.0272, 0.0512, 0.0659, 0.0754, 0.0820];
%! assert (p, expected_p, 0.001);
%!test
%! x = 0:0.5:2.5;
%! p = ricecdf (x, 0, 1);
%! expected_p = [0.0000, 0.1175, 0.3935, 0.6753, 0.8647, 0.9561];
%! assert (p, expected_p, 0.001);
%!test
%! x = 0:0.5:2.5;
%! p = ricecdf (x, 1, 1);
%! expected_p = [0.0000, 0.0735, 0.2671, 0.5120, 0.7310, 0.8791];
%! assert (p, expected_p, 0.001);
%!shared x, p
%! x = [-1, 0, 1, 2, Inf];
%! p = [0, 0, 0.26712019620318, 0.73098793996409, 1];
%!assert (ricecdf (x, 1, 1), p, 1e-14)
%!assert (ricecdf (x, 1, 1, "upper"), 1 - p, 1e-14)

## Test input validation
%!error<ricecdf: function called with too few input arguments.> ricecdf ()
%!error<ricecdf: function called with too few input arguments.> ricecdf (1)
%!error<ricecdf: function called with too few input arguments.> ricecdf (1, 2)
%!error<ricecdf: invalid argument for upper tail.> ricecdf (1, 2, 3, "uper")
%!error<ricecdf: invalid argument for upper tail.> ricecdf (1, 2, 3, 4)
%!error<ricecdf: X, NU, and SIGMA must be of common size or scalars.> ...
%! ricecdf (ones (3), ones (2), ones (2))
%!error<ricecdf: X, NU, and SIGMA must be of common size or scalars.> ...
%! ricecdf (ones (2), ones (3), ones (2))
%!error<ricecdf: X, NU, and SIGMA must be of common size or scalars.> ...
%! ricecdf (ones (2), ones (2), ones (3))
%!error<ricecdf: X, NU, and SIGMA must not be complex.> ricecdf (i, 2, 3)
%!error<ricecdf: X, NU, and SIGMA must not be complex.> ricecdf (2, i, 3)
%!error<ricecdf: X, NU, and SIGMA must not be complex.> ricecdf (2, 2, i)
