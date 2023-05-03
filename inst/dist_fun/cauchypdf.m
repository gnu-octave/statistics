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
## @deftypefn  {statistics} {@var{y} =} cauchypdf (@var{x}, @var{x0}, @var{gamma})
##
## Cauchy probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Cauchy distribution with location parameter @var{x0} and
## scale parameter @var{gamma}.  The size of @var{y} is the common size of
## @var{x}, @var{x0}, and @var{gamma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Further information about the Cauchy distribution can be found at
## @url{https://en.wikipedia.org/wiki/Cauchy_distribution}
##
## @seealso{cauchycdf, cauchypdf, cauchyrnd}
## @end deftypefn

function y = cauchypdf (x, x0, gamma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("cauchypdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, X0, and GAMMA
  if (! isscalar (x) || ! isscalar (x0) || ! isscalar (gamma))
    [retval, x, x0, gamma] = common_size (x, x0, gamma);
    if (retval > 0)
      error (strcat (["cauchypdf: X, X0, and GAMMA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, X0, and GAMMA being reals
  if (iscomplex (x) || iscomplex (x0) || iscomplex (gamma))
    error ("cauchypdf: X, X0, and GAMMA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (x0, "single") || isa (gamma, "single"))
    y = NaN (size (x), "single");
  else
    y = NaN (size (x));
  endif

  ## Find valid values in parameters
  k = ! isinf (x0) & (gamma > 0) & (gamma < Inf);

  if (isscalar (x0) && isscalar (gamma))
    y(k) = ((1 ./ (1 + ((x(k) - x0) / gamma) .^ 2)) / pi / gamma);
  else
    y(k) = ((1 ./ (1 + ((x(k) - x0(k)) ./ gamma(k)) .^ 2)) / pi ./ gamma(k));
  endif

endfunction

%!demo
%! ## Plot various PDFs from the Cauchy distribution
%! x = -5:0.01:5;
%! y1 = cauchypdf (x, 0, 0.5);
%! y2 = cauchypdf (x, 0, 1);
%! y3 = cauchypdf (x, 0, 2);
%! y4 = cauchypdf (x, -2, 1);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c")
%! grid on
%! xlim ([-5, 5])
%! ylim ([0, 0.7])
%! legend ({"x0 = 0, γ = 0.5", "x0 = 0, γ = 1", ...
%!          "x0 = 0, γ = 2", "x0 = -2, γ = 1"}, "location", "northeast")
%! title ("Cauchy PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 0.5 1 2];
%! y = 1/pi * ( 2 ./ ((x-1).^2 + 2^2) );
%!assert (cauchypdf (x, ones (1,5), 2*ones (1,5)), y)
%!assert (cauchypdf (x, 1, 2*ones (1,5)), y)
%!assert (cauchypdf (x, ones (1,5), 2), y)
%!assert (cauchypdf (x, [-Inf 1 NaN 1 Inf], 2), [NaN y(2) NaN y(4) NaN])
%!assert (cauchypdf (x, 1, 2*[0 1 NaN 1 Inf]), [NaN y(2) NaN y(4) NaN])
%!assert (cauchypdf ([x, NaN], 1, 2), [y, NaN])

## Test class of input preserved
%!assert (cauchypdf (single ([x, NaN]), 1, 2), single ([y, NaN]), eps ("single"))
%!assert (cauchypdf ([x, NaN], single (1), 2), single ([y, NaN]), eps ("single"))
%!assert (cauchypdf ([x, NaN], 1, single (2)), single ([y, NaN]), eps ("single"))

## Cauchy (0,1) == Student's T distribution with 1 DOF
%!test
%! x = rand (10, 1);
%! assert (cauchypdf (x, 0, 1), tpdf (x, 1), eps);

## Test input validation
%!error<cauchypdf: function called with too few input arguments.> cauchypdf ()
%!error<cauchypdf: function called with too few input arguments.> cauchypdf (1)
%!error<cauchypdf: function called with too few input arguments.> ...
%! cauchypdf (1, 2)
%!error<cauchypdf: function called with too many inputs> cauchypdf (1, 2, 3, 4)
%!error<cauchypdf: X, X0, and GAMMA must be of common size or scalars.> ...
%! cauchypdf (ones (3), ones (2), ones(2))
%!error<cauchypdf: X, X0, and GAMMA must be of common size or scalars.> ...
%! cauchypdf (ones (2), ones (3), ones(2))
%!error<cauchypdf: X, X0, and GAMMA must be of common size or scalars.> ...
%! cauchypdf (ones (2), ones (2), ones(3))
%!error<cauchypdf: X, X0, and GAMMA must not be complex.> cauchypdf (i, 4, 3)
%!error<cauchypdf: X, X0, and GAMMA must not be complex.> cauchypdf (1, i, 3)
%!error<cauchypdf: X, X0, and GAMMA must not be complex.> cauchypdf (1, 4, i)
