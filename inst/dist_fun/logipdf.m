## Copyright (C) 1995-2017 Kurt Hornik
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
## @deftypefn  {statistics} {@var{y} =} logipdf (@var{x})
## @deftypefnx {statistics} {@var{y} =} logipdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{y} =} logipdf (@var{x}, @var{mu}, @var{s})
##
## Logistic probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the logistic distribution with location parameter @var{mu} and scale
## parameter @var{s}.  The size of @var{p} is the common size of @var{x},
## @var{mu}, and @var{s}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Both parameters must be reals and @qcode{@var{s} > 0}.
## For @qcode{@var{s} <= 0}, @qcode{NaN} is returned.
##
## Further information about the logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Logistic_distribution}
##
## @seealso{logicdf, logiinv, logirnd, logifit, logilike, logistat}
## @end deftypefn

function y = logipdf (x, mu, s)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("logipdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, MU, and S
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(s))
    [retval, x, mu, s] = common_size (x, mu, s);
    if (retval > 0)
      error (strcat (["logipdf: X, MU, and S must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, MU, and S being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (s))
    error ("logipdf: X, MU, and S must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (mu, "single") || isa (s, "single"));
    y = NaN (size (x), "single");
  else
    y = NaN (size (x));
  endif

  ## Compute logistic PDF
  k1 = ((x == -Inf) & (s > 0)) | ((x == Inf) & (s > 0));
  y(k1) = 0;

  k = ! k1 & (s > 0);
  y(k) = (1 ./ (4 .* s(k))) .* ...
         (sech ((x(k) - mu(k)) ./ (2 .* s(k))) .^ 2);

endfunction

%!demo
%! ## Plot various PDFs from the logistic distribution
%! x = -5:0.01:20;
%! y1 = logipdf (x, 5, 2);
%! y2 = logipdf (x, 9, 3);
%! y3 = logipdf (x, 9, 4);
%! y4 = logipdf (x, 6, 2);
%! y5 = logipdf (x, 2, 1);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
%! grid on
%! ylim ([0, 0.3])
%! legend ({"μ = 5, s = 2", "μ = 9, s = 3", "μ = 2, s = 4", ...
%!          "μ = 6, s = 2", "μ = 2, s = 1"}, "location", "northeast")
%! title ("Logistic PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-Inf -log(4) 0 log(4) Inf];
%! y = [0, 0.16, 1/4, 0.16, 0];
%!assert (logipdf ([x, NaN], 0, 1), [y, NaN], eps)
%!assert (logipdf (x, 0, [-2, -1, 0, 1, 2]), [nan(1, 3), y([4:5])], eps)

## Test class of input preserved
%!assert (logipdf (single ([x, NaN]), 0, 1), single ([y, NaN]), eps ("single"))
%!assert (logipdf ([x, NaN], single (0), 1), single ([y, NaN]), eps ("single"))
%!assert (logipdf ([x, NaN], 0, single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<logipdf: function called with too few input arguments.> logipdf ()
%!error<logipdf: function called with too few input arguments.> logipdf (1)
%!error<logipdf: function called with too few input arguments.> ...
%! logipdf (1, 2)
%!error<logipdf: X, MU, and S must be of common size or scalars.> ...
%! logipdf (1, ones (2), ones (3))
%!error<logipdf: X, MU, and S must be of common size or scalars.> ...
%! logipdf (ones (2), 1, ones (3))
%!error<logipdf: X, MU, and S must be of common size or scalars.> ...
%! logipdf (ones (2), ones (3), 1)
%!error<logipdf: X, MU, and S must not be complex.> logipdf (i, 2, 3)
%!error<logipdf: X, MU, and S must not be complex.> logipdf (1, i, 3)
%!error<logipdf: X, MU, and S must not be complex.> logipdf (1, 2, i)
