## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {statistics} {@var{p} =} tlspdf (@var{x}, @var{mu}, @var{sigma}, @var{df})
##
## Location-scale Student's T probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the location-scale Student's T distribution with location parameter
## @var{mu}, scale parameter @var{sigma}, and @var{df} degrees of freedom.  The
## size of @var{y} is the common size of @var{x}, @var{mu}, @var{sigma}, and
## @var{df}. A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## Further information about the location-scale Student's T distribution can be
## found at @url{https://en.wikipedia.org/wiki/Student%27s_t-distribution#Location-scale_t_distribution}
##
## @seealso{tlscdf, tlspdf, tlsrnd, tlsfit, tlslike, tlsstat}
## @end deftypefn

function y = tlspdf (x, mu, sigma, df)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("tlspdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, MU, SIGMA, and DF
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (sigma) || ! isscalar (df))
    [err, x, mu, sigma, df] = common_size (x, mu, sigma, df);
    if (err > 0)
      error ("tlspdf: X, MU, SIGMA, and DF must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, SIGMA, and DF being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma) || iscomplex (df))
    error ("tlspdf: X, MU, SIGMA, and DF must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") ||
      isa (sigma, "single") || isa (df, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Force invalid SIGMA parameter to NaN
  sigma(sigma <= 0) = NaN;

  ## Call tpdf to do the work
  y = tpdf ((x - mu) ./ sigma, df) ./ sigma;

  ## Force class type
  y = cast (y, cls);

endfunction

%!demo
%! ## Plot various PDFs from the Student's T distribution
%! x = -8:0.01:8;
%! y1 = tlspdf (x, 0, 1, 1);
%! y2 = tlspdf (x, 0, 2, 2);
%! y3 = tlspdf (x, 3, 2, 5);
%! y4 = tlspdf (x, -1, 3, Inf);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-m")
%! grid on
%! xlim ([-8, 8])
%! ylim ([0, 0.41])
%! legend ({"mu = 0, sigma = 1, df = 1", "mu = 0, sigma = 2, df = 2", ...
%!          "mu = 3, sigma = 2, df = 5", 'mu = -1, sigma = 3, df = \infty'}, ...
%!         "location", "northwest")
%! title ("Location-scale Student's T PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!test
%! x = rand (10,1);
%! y = 1./(pi * (1 + x.^2));
%! assert (tlspdf (x, 0, 1, 1), y, 5*eps);
%! assert (tlspdf (x+5, 5, 1, 1), y, 5*eps);
%! assert (tlspdf (x.*2, 0, 2, 1), y./2, 5*eps);
%!shared x, y
%! x = [-Inf 0 0.5 1 Inf];
%! y = 1./(pi * (1 + x.^2));
%!assert (tlspdf (x, 0, 1, ones (1,5)), y, eps)
%!assert (tlspdf (x, 0, 1, 1), y, eps)
%!assert (tlspdf (x, 0, 1, [0 NaN 1 1 1]), [NaN NaN y(3:5)], eps)
%!assert (tlspdf (x, 0, 1, Inf), normpdf (x))

## Test class of input preserved
%!assert (class (tlspdf ([x, NaN], 1, 1, 1)), "double")
%!assert (class (tlspdf (single ([x, NaN]), 1, 1, 1)), "single")
%!assert (class (tlspdf ([x, NaN], single (1), 1, 1)), "single")
%!assert (class (tlspdf ([x, NaN], 1, single (1), 1)), "single")
%!assert (class (tlspdf ([x, NaN], 1, 1, single (1))), "single")

## Test input validation
%!error<tlspdf: function called with too few input arguments.> tlspdf ()
%!error<tlspdf: function called with too few input arguments.> tlspdf (1)
%!error<tlspdf: function called with too few input arguments.> tlspdf (1, 2)
%!error<tlspdf: function called with too few input arguments.> tlspdf (1, 2, 3)
%!error<tlspdf: X, MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlspdf (ones (3), ones (2), 1, 1)
%!error<tlspdf: X, MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlspdf (ones (2), 1, ones (3), 1)
%!error<tlspdf: X, MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlspdf (ones (2), 1, 1, ones (3))
%!error<tlspdf: X, MU, SIGMA, and DF must not be complex.> tlspdf (i, 2, 1, 1)
%!error<tlspdf: X, MU, SIGMA, and DF must not be complex.> tlspdf (2, i, 1, 1)
%!error<tlspdf: X, MU, SIGMA, and DF must not be complex.> tlspdf (2, 1, i, 1)
%!error<tlspdf: X, MU, SIGMA, and DF must not be complex.> tlspdf (2, 1, 1, i)
