## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{y} =} laplacepdf (@var{x}, @var{mu}, @var{beta})
##
## Laplace probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Laplace distribution with location parameter @var{mu} and scale
## parameter (i.e. "diversity") @var{beta}.  The size of @var{y} is the common
## size of @var{x}, @var{mu}, and @var{beta}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## Both parameters must be reals and @qcode{@var{beta} > 0}.
## For @qcode{@var{beta} <= 0}, @qcode{NaN} is returned.
##
## Further information about the Laplace distribution can be found at
## @url{https://en.wikipedia.org/wiki/Laplace_distribution}
##
## @seealso{laplacecdf, laplacepdf, laplacernd}
## @end deftypefn

function y = laplacepdf (x, mu, beta)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("laplacepdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, MU, and BETA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(beta))
    [retval, x, mu, beta] = common_size (x, mu, beta);
    if (retval > 0)
      error (strcat (["laplacepdf: X, MU, and BETA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, MU, and BETA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (beta))
    error ("laplacepdf: X, MU, and BETA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (beta, "single"));
    y = NaN (size (x), "single");
  else
    y = NaN (size (x));
  endif

  ## Compute Laplace PDF
  k1 = ((x == -Inf) & (beta > 0)) | ((x == Inf) & (beta > 0));
  y(k1) = 0;

  k = ! k1 & (beta > 0);
  y(k) = exp (- abs (x(k) - mu(k)) ./ beta(k)) ./ (2 .* beta(k));

endfunction

%!demo
%! ## Plot various PDFs from the Laplace distribution
%! x = -10:0.01:10;
%! y1 = laplacepdf (x, 0, 1);
%! y2 = laplacepdf (x, 0, 2);
%! y3 = laplacepdf (x, 0, 4);
%! y4 = laplacepdf (x, -5, 4);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c")
%! grid on
%! xlim ([-10, 10])
%! ylim ([0, 0.6])
%! legend ({"μ = 0, β = 1", "μ = 0, β = 2", ...
%!          "μ = 0, β = 4", "μ = -5, β = 4"}, "location", "northeast")
%! title ("Laplace PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test results
%!shared x, y
%! x = [-Inf -log(2) 0 log(2) Inf];
%! y = [0, 1/4, 1/2, 1/4, 0];
%!assert (laplacepdf ([x, NaN], 0, 1), [y, NaN])
%!assert (laplacepdf (x, 0, [-2, -1, 0, 1, 2]), [nan(1, 3), 0.25, 0])

## Test class of input preserved
%!assert (laplacepdf (single ([x, NaN]), 0, 1), single ([y, NaN]))
%!assert (laplacepdf ([x, NaN], single (0), 1), single ([y, NaN]))
%!assert (laplacepdf ([x, NaN], 0, single (1)), single ([y, NaN]))

## Test input validation
%!error<laplacepdf: function called with too few input arguments.> laplacepdf ()
%!error<laplacepdf: function called with too few input arguments.> laplacepdf (1)
%!error<laplacepdf: function called with too few input arguments.> ...
%! laplacepdf (1, 2)
%!error<laplacepdf: function called with too many inputs> laplacepdf (1, 2, 3, 4)
%!error<laplacepdf: X, MU, and BETA must be of common size or scalars.> ...
%! laplacepdf (1, ones (2), ones (3))
%!error<laplacepdf: X, MU, and BETA must be of common size or scalars.> ...
%! laplacepdf (ones (2), 1, ones (3))
%!error<laplacepdf: X, MU, and BETA must be of common size or scalars.> ...
%! laplacepdf (ones (2), ones (3), 1)
%!error<laplacepdf: X, MU, and BETA must not be complex.> laplacepdf (i, 2, 3)
%!error<laplacepdf: X, MU, and BETA must not be complex.> laplacepdf (1, i, 3)
%!error<laplacepdf: X, MU, and BETA must not be complex.> laplacepdf (1, 2, i)
