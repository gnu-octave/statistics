## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} stblinv (@var{p}, @var{alpha}, @var{beta}, @var{gam}, @var{delta})
##
## Inverse of the stable cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the stable distribution with tail index (first shape parameter) @var{alpha},
## skewness (second shape parameter) @var{beta}, scale parameter @var{gam}, and
## location parameter @var{delta}, in the Nolan @qcode{S0} parameterization.
## The size of @var{x} is the size of @var{p}.
##
## @var{alpha} must be in the range @math{(0, 2]}, @var{beta} in @math{[-1, 1]},
## @var{gam} positive, and @var{delta} real.  The parameters must be scalars.
##
## The quantile has a closed form for @var{alpha} equal to @code{2} (normal) and
## for @code{1} with @var{beta} equal to @code{0} (Cauchy); otherwise it is
## found by numerical inversion of @code{stblcdf}.
##
## @seealso{stblcdf, stblpdf, stblrnd, makedist}
## @end deftypefn

function x = stblinv (p, alpha, beta, gam, delta)

  if (nargin != 5)
    print_usage ();
  endif

  msg = __stable_checkparams__ (alpha, beta, gam, delta);
  if (! isempty (msg))
    error ("stblinv: %s", msg);
  endif
  if (! (isreal (p) && all (p(:) >= 0 & p(:) <= 1 | isnan (p(:)))))
    error ("stblinv: P must contain real values in the range [0, 1].");
  endif

  x = nan (size (p));

  if (alpha == 2)
    x = delta + sqrt (2) .* gam .* (-sqrt (2) .* erfcinv (2 .* p));
  elseif (alpha == 1 && beta == 0)
    x = delta + gam .* tan (pi .* (p - 0.5));
  else
    for i = 1:numel (p)
      pp = p(i);
      if (isnan (pp))
        continue;
      elseif (pp == 0)
        x(i) = -Inf;
      elseif (pp == 1)
        x(i) = Inf;
      else
        ## Root-find on the standardized cdf, started from the (heavy-tailed)
        ## Cauchy quantile so the search stays near the solution.
        z0 = tan (pi .* (pp - 0.5));
        z = fzero (@(zz) stblcdf (zz, alpha, beta, 1, 0) - pp, z0);
        x(i) = delta + gam .* z;
      endif
    endfor
  endif

endfunction

%!demo
%! ## Quantiles of a skewed stable distribution
%! p = [0.1, 0.25, 0.5, 0.75, 0.9];
%! x = stblinv (p, 1.5, 0.5, 1, 0)

## Test output against MATLAB
%!test
%! p = [0.1, 0.25, 0.5, 0.75, 0.9];
%! x = stblinv (p, 1.5, 0.5, 1, 0);
%! exp_x = [-1.63127009138493, -0.783313648587273, 0.133853042315326, ...
%!          1.20341055131626, 2.58231785139714];
%! assert (x, exp_x, 1e-6);
%!test
%! p = [0.1, 0.25, 0.5, 0.75, 0.9];
%! x = stblinv (p, 0.8, 0.5, 1, 0);
%! exp_x = [-1.62200033048034, -0.553853652413272, 0.250487323305453, ...
%!          2.11601429745393, 8.03924696271835];
%! assert (x, exp_x, 1e-5);
%!test  # scaled and shifted (gam = 2, delta = 3)
%! p = [0.1, 0.25, 0.5, 0.75, 0.9];
%! x = stblinv (p, 1.5, 0.5, 2, 3);
%! exp_x = [-0.26254018276987, 1.43337270282545, 3.26770608463065, ...
%!          5.40682110263252, 8.16463570279428];
%! assert (x, exp_x, 1e-6);
%!test  # symmetric case (beta = 0): quantiles antisymmetric about delta
%! p = [0.1, 0.25, 0.5, 0.75, 0.9];
%! x = stblinv (p, 1.5, 0, 1, 0);
%! exp_x = [-2.06146263813919, -0.968933181710917, 0, ...
%!          0.968933181710917, 2.06146263813919];
%! assert (x, exp_x, 1e-6);
%!test  # normal and Cauchy special cases
%! p = [0.1, 0.3, 0.5, 0.7, 0.9];
%! assert (stblinv (p, 2, 0, 1, 0), norminv (p, 0, sqrt (2)), 1e-12);
%! assert (stblinv (p, 1, 0, 1, 0), tan (pi .* (p - 0.5)), 1e-12);
%!test  # inverts stblcdf
%! x0 = [-3, -0.5, 0.8, 4];
%! p = stblcdf (x0, 1.4, 0.3, 1.5, -1);
%! assert (stblinv (p, 1.4, 0.3, 1.5, -1), x0, 1e-6);
%!test  # boundaries
%! assert (stblinv ([0, 1], 1.5, 0.5, 1, 0), [-Inf, Inf]);

## Test input validation
%!error <Invalid call to stblinv> stblinv (0.5, 1.5, 0.5, 1)
%!error <stblinv: ALPHA must be a scalar in the range \(0, 2\].> ...
%! stblinv (0.5, 2.5, 0, 1, 0)
%!error <stblinv: P must contain real values in the range \[0, 1\].> ...
%! stblinv (1.2, 1.5, 0, 1, 0)
%!error <stblinv: P must contain real values in the range \[0, 1\].> ...
%! stblinv (0.5i, 1.5, 0, 1, 0)
