## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} loglpdf (@var{x}, @var{a}, @var{b})
##
## Log-logistic probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the log-logistic distribution with with scale parameter @var{a} and shape
## parameter @var{b}.  The size of @var{y} is the common size of @var{x},
## @var{a}, and @var{b}.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Both parameters, @var{a} and @var{b}, must be positive reals, otherwise
## @qcode{NaN} is returned.   @var{x} is supported in the range @math{[0,Inf)},
## otherwise @qcode{0} is returned.
##
## Further information about the log-logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-logistic_distribution}
##
## MATLAB compatibility: MATLAB uses an alternative parameterization given by
## the pair @math{μ, s}, i.e. @var{mu} and @var{s}, in analogy with the logistic
## distribution.  Their relation to the @var{a} and @var{b} parameters is given
## below:
##
## @itemize
## @item @qcode{@var{a} = exp (@var{mu})}
## @item @qcode{@var{b} = 1 / @var{s}}
## @end itemize
##
## @seealso{loglcdf, loglinv, loglrnd, loglfit, logllike}
## @end deftypefn

function y = loglpdf (x, a, b)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("loglpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, A, and B
  if (! isscalar (x) || ! isscalar (a) || ! isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error (strcat (["loglpdf: X, A, and B must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, A, and B being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("loglpdf: X, A, and B must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single"));
    y = NaN (size (x), "single");
  else
    y = NaN (size (x));
  endif

  ## Compute log-logistic PDF
  k1 = ((x == Inf) | (x < 0)) & (a > 0) & (b > 0);
  y(k1) = 0;

  k = (! k1) & (a > 0) & (b > 0);
  y(k) = ((b(k) ./ a(k)) .* (x(k) ./ a(k)) .^ (b(k) -1)) ./ ...
         ((1 + (x(k) ./ a(k)) .^ b(k)) .^ 2);

endfunction

%!demo
%! ## Plot various PDFs from the log-logistic distribution
%! x = 0:0.001:2;
%! y1 = loglpdf (x, 1, 0.5);
%! y2 = loglpdf (x, 1, 1);
%! y3 = loglpdf (x, 1, 2);
%! y4 = loglpdf (x, 1, 4);
%! y5 = loglpdf (x, 1, 8);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
%! grid on
%! ylim ([0,3])
%! legend ({"β = 0.5", "β = 1", "β = 2", "β = 4", "β = 8"}, ...
%!         "location", "northeast")
%! title ("Log-logistic PDF")
%! xlabel ("values in x")
%! ylabel ("density")
%! text (0.5, 2.8, "α = 1, values of β as shown in legend")

## Test output
%!shared out1, out2
%! out1 = [0, 1, 0.2500, 0.1111, 0.0625, 0.0400, 0.0278, 0];
%! out2 = [0, Inf, 0.0811, 0.0416, 0.0278, 0.0207, 0.0165, 0];
%!assert (loglpdf ([-1:5,Inf], 1, 1), out1, 1e-4)
%!assert (loglpdf ([-1:5,Inf], exp (0), 1), out1, 1e-4)
%!assert (loglpdf ([-1:5,Inf], exp (1), 1 / 3), out2, 1e-4)

## Test class of input preserved
%!assert (class (loglpdf (single (1), 2, 3)), "single")
%!assert (class (loglpdf (1, single (2), 3)), "single")
%!assert (class (loglpdf (1, 2, single (3))), "single")

## Test input validation
%!error<loglpdf: function called with too few input arguments.> loglpdf (1)
%!error<loglpdf: function called with too few input arguments.> loglpdf (1, 2)
%!error<loglpdf: X, A, and B must be of common size or scalars.> ...
%! loglpdf (1, ones (2), ones (3))
%!error<loglpdf: X, A, and B must be of common size or scalars.> ...
%! loglpdf (ones (2), 1, ones (3))
%!error<loglpdf: X, A, and B must be of common size or scalars.> ...
%! loglpdf (ones (2), ones (3), 1)
%!error<loglpdf: X, A, and B must not be complex.> loglpdf (i, 2, 3)
%!error<loglpdf: X, A, and B must not be complex.> loglpdf (1, i, 3)
%!error<loglpdf: X, A, and B must not be complex.> loglpdf (1, 2, i)
