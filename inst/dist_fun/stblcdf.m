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
## @deftypefn  {statistics} {@var{p} =} stblcdf (@var{x}, @var{alpha}, @var{beta}, @var{gam}, @var{delta})
##
## Stable cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the stable distribution with tail index (first shape parameter)
## @var{alpha}, skewness (second shape parameter) @var{beta}, scale parameter
## @var{gam}, and location parameter @var{delta}, in the Nolan @qcode{S0}
## parameterization.  The size of @var{p} is the size of @var{x}.
##
## @var{alpha} must be in the range @math{(0, 2]}, @var{beta} in @math{[-1, 1]},
## @var{gam} positive, and @var{delta} real.  The parameters must be scalars.
##
## The cumulative probability has a closed form for @var{alpha} equal to
## @code{2} (normal) and for @code{1} with @var{beta} equal to @code{0}
## (Cauchy); otherwise it is computed by numerical inversion of the
## characteristic function (the Gil-Pelaez formula).
##
## @seealso{stblpdf, stblinv, stblrnd, makedist}
## @end deftypefn

function p = stblcdf (x, alpha, beta, gam, delta)

  if (nargin != 5)
    print_usage ();
  endif

  msg = __stable_checkparams__ (alpha, beta, gam, delta);
  if (! isempty (msg))
    error ("stblcdf: %s", msg);
  endif
  if (! isreal (x))
    error ("stblcdf: X must be real.");
  endif

  z = (x - delta) ./ gam;
  p = nan (size (z));
  ok = ! isnan (z);

  if (alpha == 2)
    p(ok) = 0.5 .* erfc (-z(ok) ./ 2);
  elseif (alpha == 1 && beta == 0)
    p(ok) = 0.5 + atan (z(ok)) ./ pi;
  else
    zk = z(ok);
    v = zeros (numel (zk), 1);
    for i = 1:numel (zk)
      v(i) = 0.5 - (1 ./ pi) .* quadgk (@(t) imag (exp (-1i .* t .* zk(i)) ...
                          .* __stable_cf__ (t, alpha, beta)) ./ t, 0, Inf, ...
                          "AbsTol", 1e-12, "RelTol", 1e-10);
    endfor
    p(ok) = min (max (v, 0), 1);
  endif

endfunction

%!demo
%! ## Stable cdf: Cauchy, a skewed stable, and the normal limit
%! x = linspace (-6, 6, 200);
%! plot (x, stblcdf (x, 1, 0, 1, 0), "-", ...
%!       x, stblcdf (x, 1.5, 0.5, 1, 0), "-", ...
%!       x, stblcdf (x, 2, 0, 1, 0), "-");
%! legend ("Cauchy", "alpha=1.5, beta=0.5", "normal", "location", "southeast");

## Test output against MATLAB
%!test
%! x = -5:5;
%! p = stblcdf (x, 1.5, 0.5, 1, 0);
%! exp_p = [0.00961772128347771, 0.0143422747723476, 0.0257902242195547, ...
%!          0.0657154294128386, 0.201576145758624, 0.462186560100778, ...
%!          0.712063555515659, 0.855535196378772, 0.921201224725992, ...
%!          0.951409668616683, 0.966845678836178];
%! assert (p, exp_p, 1e-8);
%!test
%! x = -5:5;
%! p = stblcdf (x, 0.8, 0.5, 1, 0);
%! exp_p = [0.0427283102762096, 0.0504255041089544, 0.0624716830177048, ...
%!          0.0849086757683013, 0.150275591315296, 0.431333711402679, ...
%!          0.641248581720908, 0.74188789948898, 0.797926083610673, ...
%!          0.833292740693924, 0.857610462691116];
%! assert (p, exp_p, 1e-8);
%!test  # scaled and shifted (gam = 2, delta = 3)
%! x = -5:5;
%! p = stblcdf (x, 1.5, 0.5, 2, 3);
%! exp_p = [0.0143422747723476, 0.0186030552365194, 0.0257902242195547, ...
%!          0.0392075905274278, 0.0657154294128386, 0.116299801968237, ...
%!          0.201576145758624, 0.321987153858349, 0.462186560100778, ...
%!          0.598389078433622, 0.712063555515659];
%! assert (p, exp_p, 1e-8);
%!test  # normal special case
%! x = -5:5;
%! assert (stblcdf (x, 2, 0, 1, 0), normcdf (x, 0, sqrt (2)), 1e-12);
%!test  # Cauchy special case
%! x = -5:5;
%! assert (stblcdf (x, 1, 0, 1, 0), 0.5 + atan (x) ./ pi, 1e-12);
%!test  # cdf is the integral of the pdf
%! assert (stblcdf (0.7, 1.3, -0.4, 1, 0) - stblcdf (-1.2, 1.3, -0.4, 1, 0), ...
%!         quadgk (@(x) stblpdf (x, 1.3, -0.4, 1, 0), -1.2, 0.7), 1e-8);

## Test input validation
%!error <Invalid call to stblcdf> stblcdf (1, 1.5, 0.5, 1)
%!error <stblcdf: ALPHA must be a scalar in the range \(0, 2\].> ...
%! stblcdf (1, 2.5, 0, 1, 0)
%!error <stblcdf: BETA must be a scalar in the range \[-1, 1\].> ...
%! stblcdf (1, 1.5, 2, 1, 0)
%!error <stblcdf: GAM must be a positive scalar.> stblcdf (1, 1.5, 0, 0, 0)
%!error <stblcdf: DELTA must be a real scalar.> stblcdf (1, 1.5, 0, 1, 1i)
%!error <stblcdf: X must be real.> stblcdf (1i, 1.5, 0, 1, 0)
