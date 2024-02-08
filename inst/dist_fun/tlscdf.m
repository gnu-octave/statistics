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
## @deftypefn  {statistics} {@var{p} =} tlscdf (@var{x}, @var{mu}, @var{sigma}, @var{df})
## @deftypefnx {statistics} {@var{p} =} tlscdf (@var{x}, @var{mu}, @var{sigma}, @var{df}, @qcode{"upper"})
##
## Location-scale Student's T cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the location-scale Student's T distribution with location parameter
## @var{mu}, scale parameter @var{sigma}, and @var{df} degrees of freedom.  The
## size of @var{p} is the common size of @var{x}, @var{mu}, @var{sigma}, and
## @var{df}. A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## @code{@var{p} = tlscdf (@var{x}, @var{mu}, @var{sigma}, @var{df}, "upper")}
## computes the upper tail probability of the location-scale Student's T
## distribution with parameters @var{mu}, @var{sigma}, and @var{df}, at the
## values in @var{x}.
##
## Further information about the location-scale Student's T distribution can be
## found at @url{https://en.wikipedia.org/wiki/Student%27s_t-distribution#Location-scale_t_distribution}
##
## @seealso{tlsinv, tlspdf, tlsrnd, tlsfit, tlslike, tlsstat}
## @end deftypefn

function p = tlscdf (x, mu, sigma, df, uflag)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("tlscdf: function called with too few input arguments.");
  endif

  ## Check for "upper" flag
  upper = false;
  if (nargin > 4 && strcmpi (uflag, "upper"))
    upper = true;
  elseif (nargin > 4  && ! strcmpi (uflag, "upper"))
    error ("tlscdf: invalid argument for upper tail.");
  endif

  ## Check for common size of X, MU, SIGMA, and DF
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (sigma) || ! isscalar (df))
    [err, x, mu, sigma, df] = common_size (x, mu, sigma, df);
    if (err > 0)
      error ("tlscdf: X, MU, SIGMA, and DF must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, SIGMA, and DF being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma) || iscomplex (df))
    error ("tlscdf: X, MU, SIGMA, and DF must not be complex.");
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

  ## Call tcdf to do the work
  if (upper)
    p = tcdf ((x - mu) ./ sigma, df, "upper");
  else
    p = tcdf ((x - mu) ./ sigma, df);
  endif

  ## Force class type
  p = cast (p, cls);

endfunction

%!demo
%! ## Plot various CDFs from the location-scale Student's T distribution
%! x = -8:0.01:8;
%! p1 = tlscdf (x, 0, 1, 1);
%! p2 = tlscdf (x, 0, 2, 2);
%! p3 = tlscdf (x, 3, 2, 5);
%! p4 = tlscdf (x, -1, 3, Inf);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-m")
%! grid on
%! xlim ([-8, 8])
%! ylim ([0, 1])
%! legend ({"mu = 0, sigma = 1, df = 1", "mu = 0, sigma = 2, df = 2", ...
%!          "mu = 3, sigma = 2, df = 5", 'mu = -1, sigma = 3, df = \infty'}, ...
%!         "location", "northwest")
%! title ("Location-scale Student's T CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x,y
%! x = [-Inf 0 1 Inf];
%! y = [0 1/2 3/4 1];
%!assert (tlscdf (x, 0, 1, ones (1,4)), y, eps)
%!assert (tlscdf (x, 0, 1, 1), y, eps)
%!assert (tlscdf (x, 0, 1, [0 1 NaN 1]), [NaN 1/2 NaN 1], eps)
%!assert (tlscdf ([x(1:2) NaN x(4)], 0, 1, 1), [y(1:2) NaN y(4)], eps)
%!assert (tlscdf (2, 0, 1, 3, "upper"), 0.0697, 1e-4)
%!assert (tlscdf (205, 0, 1, 5, "upper"), 2.6206e-11, 1e-14)

## Test class of input preserved
%!assert (tlscdf ([x, NaN], 0, 1, 1), [y, NaN], eps)
%!assert (tlscdf (single ([x, NaN]), 0, 1, 1), single ([y, NaN]), eps ("single"))
%!assert (tlscdf ([x, NaN], single (0), 1, 1), single ([y, NaN]), eps ("single"))
%!assert (tlscdf ([x, NaN], 0, single (1), 1), single ([y, NaN]), eps ("single"))
%!assert (tlscdf ([x, NaN], 0, 1, single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<tlscdf: function called with too few input arguments.> tlscdf ()
%!error<tlscdf: function called with too few input arguments.> tlscdf (1)
%!error<tlscdf: function called with too few input arguments.> tlscdf (1, 2)
%!error<tlscdf: function called with too few input arguments.> tlscdf (1, 2, 3)
%!error<tlscdf: invalid argument for upper tail.> tlscdf (1, 2, 3, 4, "uper")
%!error<tlscdf: invalid argument for upper tail.> tlscdf (1, 2, 3, 4, 5)
%!error<tlscdf: X, MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlscdf (ones (3), ones (2), 1, 1)
%!error<tlscdf: X, MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlscdf (ones (3), 1, ones (2), 1)
%!error<tlscdf: X, MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlscdf (ones (3), 1, 1, ones (2))
%!error<tlscdf: X, MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlscdf (ones (3), ones (2), 1, 1, "upper")
%!error<tlscdf: X, MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlscdf (ones (3), 1, ones (2), 1, "upper")
%!error<tlscdf: X, MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlscdf (ones (3), 1, 1, ones (2), "upper")
%!error<tlscdf: X, MU, SIGMA, and DF must not be complex.> tlscdf (i, 2, 1, 1)
%!error<tlscdf: X, MU, SIGMA, and DF must not be complex.> tlscdf (2, i, 1, 1)
%!error<tlscdf: X, MU, SIGMA, and DF must not be complex.> tlscdf (2, 1, i, 1)
%!error<tlscdf: X, MU, SIGMA, and DF must not be complex.> tlscdf (2, 1, 1, i)
