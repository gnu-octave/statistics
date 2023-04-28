## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{p} =} cauchycdf (@var{x}, @var{x0}, @var{gamma})
## @deftypefnx {statistics} {@var{p} =} cauchycdf (@var{x}, @var{x0}, @var{gamma}, @qcode{"upper"})
##
## Cauchy cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the Cauchy distribution with location parameter @var{x0}
## and scale parameter @var{gamma}.  The size of @var{p} is the common size of
## @var{x}, @var{x0}, and @var{gamma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## @code{@var{p} = cauchycdf (@var{x}, @var{x0}, @var{gamma}, "upper")} computes
## the upper tail probability of the Cauchy distribution with parameters
## @var{x0} and @var{gamma} at the values in @var{x}.
##
## Further information about the Cauchy distribution can be found at
## @url{https://en.wikipedia.org/wiki/Cauchy_distribution}
##
## @seealso{cauchyinv, cauchypdf, cauchyrnd}
## @end deftypefn

function p = cauchycdf (x, x0, gamma, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("cauchycdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("cauchycdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, X0, and GAMMA
  if (! isscalar (x) || ! isscalar (x0) || ! isscalar (gamma))
    [retval, x, x0, gamma] = common_size (x, x0, gamma);
    if (retval > 0)
      error (strcat (["cauchycdf: X, X0, and GAMMA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, X0, and GAMMA being reals
  if (iscomplex (x) || iscomplex (x0) || iscomplex (gamma))
    error ("cauchycdf: X, X0, and GAMMA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (x0, "single") || isa (gamma, "single"));
    p = NaN (size (x), "single");
  else
    p = NaN (size (x));
  endif

  ## Find valid values in parameters and data
  k = ! isinf (x0) & (gamma > 0) & (gamma < Inf);
  if (isscalar (x0) && isscalar (gamma))
    if (uflag)
      p = 0.5 + atan ((-x(k) + x0) / gamma) / pi;
    else
      p = 0.5 + atan ((x(k) - x0) / gamma) / pi;
    endif
  else
    if (uflag)
      p(k) = 0.5 + atan ((-x(k) + x0(k)) ./ gamma(k)) / pi;
    else
      p(k) = 0.5 + atan ((x(k) - x0(k)) ./ gamma(k)) / pi;
    endif
  endif

endfunction

%!demo
%! ## Plot various CDFs from the Cauchy distribution
%! x = -5:0.01:5;
%! p1 = cauchycdf (x, 0, 0.5);
%! p2 = cauchycdf (x, 0, 1);
%! p3 = cauchycdf (x, 0, 2);
%! p4 = cauchycdf (x, -2, 1);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c")
%! grid on
%! xlim ([-5, 5])
%! legend ({"x0 = 0, γ = 0.5", "x0 = 0, γ = 1", ...
%!          "x0 = 0, γ = 2", "x0 = -2, γ = 1"}, "location", "southeast")
%! title ("Cauchy CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test results
%!shared x, y
%! x = [-1 0 0.5 1 2];
%! y = 1/pi * atan ((x-1) / 2) + 1/2;
%!assert (cauchycdf (x, ones (1,5), 2*ones (1,5)), y)
%!assert (cauchycdf (x, 1, 2*ones (1,5)), y)
%!assert (cauchycdf (x, ones (1,5), 2), y)
%!assert (cauchycdf (x, [-Inf 1 NaN 1 Inf], 2), [NaN y(2) NaN y(4) NaN])
%!assert (cauchycdf (x, 1, 2*[0 1 NaN 1 Inf]), [NaN y(2) NaN y(4) NaN])
%!assert (cauchycdf ([x(1:2) NaN x(4:5)], 1, 2), [y(1:2) NaN y(4:5)])

## Test class of input preserved
%!assert (cauchycdf ([x, NaN], 1, 2), [y, NaN])
%!assert (cauchycdf (single ([x, NaN]), 1, 2), single ([y, NaN]), eps ("single"))
%!assert (cauchycdf ([x, NaN], single (1), 2), single ([y, NaN]), eps ("single"))
%!assert (cauchycdf ([x, NaN], 1, single (2)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<cauchycdf: function called with too few input arguments.> cauchycdf ()
%!error<cauchycdf: function called with too few input arguments.> cauchycdf (1)
%!error<cauchycdf: function called with too few input arguments.> ...
%! cauchycdf (1, 2)
%!error<cauchycdf: function called with too many inputs> ...
%! cauchycdf (1, 2, 3, 4, 5)
%!error<cauchycdf: invalid argument for upper tail.> cauchycdf (1, 2, 3, "tail")
%!error<cauchycdf: X, X0, and GAMMA must be of common size or scalars.> ...
%! cauchycdf (ones (3), ones (2), ones (2))
%!error<cauchycdf: X, X0, and GAMMA must be of common size or scalars.> ...
%! cauchycdf (ones (2), ones (3), ones (2))
%!error<cauchycdf: X, X0, and GAMMA must be of common size or scalars.> ...
%! cauchycdf (ones (2), ones (2), ones (3))
%!error<cauchycdf: X, X0, and GAMMA must not be complex.> cauchycdf (i, 2, 2)
%!error<cauchycdf: X, X0, and GAMMA must not be complex.> cauchycdf (2, i, 2)
%!error<cauchycdf: X, X0, and GAMMA must not be complex.> cauchycdf (2, 2, i)
