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
## @deftypefn  {statistics} {@var{y} =} stblpdf (@var{x}, @var{alpha}, @var{beta}, @var{gam}, @var{delta})
##
## Stable probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the stable distribution with tail index (first shape parameter)
## @var{alpha}, skewness (second shape parameter) @var{beta}, scale parameter
## @var{gam}, and location parameter @var{delta}, in the Nolan @qcode{S0}
## parameterization.  The size of @var{y} is the size of @var{x}.
##
## @var{alpha} must be in the range @math{(0, 2]}, @var{beta} in @math{[-1, 1]},
## @var{gam} positive, and @var{delta} real.  The parameters must be scalars.
##
## The density has a closed form for @var{alpha} equal to @code{2} (normal) and
## for @code{1} with @var{beta} equal to @code{0} (Cauchy); otherwise it is
## computed by numerical inversion of the characteristic function.
##
## @seealso{stblcdf, stblinv, stblrnd, makedist}
## @end deftypefn

function y = stblpdf (x, alpha, beta, gam, delta)

  if (nargin != 5)
    print_usage ();
  endif

  msg = __stable_checkparams__ (alpha, beta, gam, delta);
  if (! isempty (msg))
    error ("stblpdf: %s", msg);
  endif
  if (! isreal (x))
    error ("stblpdf: X must be real.");
  endif

  ## Standardize
  z = (x - delta) ./ gam;
  y = nan (size (z));
  ok = ! isnan (z);

  if (alpha == 2)
    ## Normal with variance 2
    y(ok) = exp (-z(ok) .^ 2 ./ 4) ./ (2 .* sqrt (pi));
  elseif (alpha == 1 && beta == 0)
    ## Cauchy
    y(ok) = 1 ./ (pi .* (1 + z(ok) .^ 2));
  else
    zk = z(ok);
    v = zeros (numel (zk), 1);
    for i = 1:numel (zk)
      v(i) = (1 ./ pi) .* quadgk (@(t) real (exp (-1i .* t .* zk(i)) ...
                          .* __stable_cf__ (t, alpha, beta)), 0, Inf, ...
                          "AbsTol", 1e-12, "RelTol", 1e-10);
    endfor
    y(ok) = max (v, 0);
  endif

  ## Undo the scale (density transforms by 1 / gam)
  y = y ./ gam;

endfunction

%!demo
%! ## Stable densities: Cauchy, a skewed stable, and the normal limit
%! x = linspace (-6, 6, 200);
%! plot (x, stblpdf (x, 1, 0, 1, 0), "-", ...
%!       x, stblpdf (x, 1.5, 0.5, 1, 0), "-", ...
%!       x, stblpdf (x, 2, 0, 1, 0), "-");
%! legend ("Cauchy", "alpha=1.5, beta=0.5", "normal");

## Test output against MATLAB
%!test
%! x = -5:5;
%! y = stblpdf (x, 1.5, 0.5, 1, 0);
%! exp_y = [0.00330549826030791, 0.00673588721821526, 0.0190320671951022, ...
%!          0.0729514702833168, 0.208194435543156, 0.284283800988578, ...
%!          0.198573023913399, 0.0958317325744725, 0.0428461930184788, ...
%!          0.0207819141087301, 0.0113306451818624];
%! assert (y, exp_y, 1e-9);
%!test
%! x = -5:5;
%! y = stblpdf (x, 0.8, 0.5, 1, 0);
%! exp_y = [0.00634335874934447, 0.0093623044752761, 0.0155686941305108, ...
%!          0.0326882516316453, 0.135673711418341, 0.298698147231422, ...
%!          0.139071606104264, 0.0722555300098094, 0.0433943212392174, ...
%!          0.0288186423689968, 0.020522989417733];
%! assert (y, exp_y, 1e-9);
%!test
%! x = -5:5;
%! y = stblpdf (x, 1.2, -0.5, 1, 0);
%! exp_y = [0.0166464288580291, 0.0264743439008481, 0.0459799636082411, ...
%!          0.0881296016218348, 0.177627321920986, 0.288106176914537, ...
%!          0.196803906514695, 0.0520585692918225, 0.0173156565634881, ...
%!          0.00836895075945555, 0.00490202402183536];
%! assert (y, exp_y, 1e-9);
%!test  # scaled and shifted (gam = 2, delta = 3)
%! x = -5:5;
%! y = stblpdf (x, 1.5, 0.5, 2, 3);
%! exp_y = [0.00336794360910763, 0.00537726811151198, 0.00951603359755111, ...
%!          0.0184406959152125, 0.0364757351416584, 0.0666533040480966, ...
%!          0.104097217771578, 0.134023248277231, 0.142141900494289, ...
%!          0.127056343301115, 0.0992865119566997];
%! assert (y, exp_y, 1e-9);
%!test  # normal special case (alpha = 2)
%! x = -5:5;
%! assert (stblpdf (x, 2, 0, 1, 0), normpdf (x, 0, sqrt (2)), 1e-12);
%!test  # Cauchy special case (alpha = 1, beta = 0)
%! x = -5:5;
%! assert (stblpdf (x, 1, 0, 1, 0), 1 ./ (pi .* (1 + x .^ 2)), 1e-12);
%!test  # Levy (alpha = 0.5, beta = 1): S0 support boundary at x = -1
%! x = -5:5;
%! y = stblpdf (x, 0.5, 1, 1, 0);
%! assert (y(x < 0), zeros (1, 5), 1e-6);
%! assert (y(x == 0), 0.241970724519143, 1e-9);

## Test input validation
%!error <Invalid call to stblpdf> stblpdf (1, 1.5, 0.5, 1)
%!error <stblpdf: ALPHA must be a scalar in the range \(0, 2\].> ...
%! stblpdf (1, 2.5, 0, 1, 0)
%!error <stblpdf: BETA must be a scalar in the range \[-1, 1\].> ...
%! stblpdf (1, 1.5, 2, 1, 0)
%!error <stblpdf: GAM must be a positive scalar.> stblpdf (1, 1.5, 0, 0, 0)
%!error <stblpdf: DELTA must be a real scalar.> stblpdf (1, 1.5, 0, 1, 1i)
%!error <stblpdf: X must be real.> stblpdf (1i, 1.5, 0, 1, 0)
