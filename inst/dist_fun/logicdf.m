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
## @deftypefn  {statistics} {@var{p} =} logicdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} logicdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} logicdf (@var{x}, @var{mu}, @var{s})
##
## Logistic cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the logistic distribution with location parameter
## @var{mu} and scale parameter @var{s}.  The size of @var{p} is the common size
## of @var{x}, @var{mu}, and @var{s}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Default values are @qcode{@var{mu} = 0} and @qcode{@var{s} = 1}.
## Both parameters must be reals and @qcode{@var{beta} > 0}.
## For @qcode{@var{beta} <= 0}, @qcode{NaN} is returned.
##
## Further information about the log-logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Logistic_distribution}
##
## @seealso{logiinv, logipdf, logirnd, logifit, logilike, logistat}
## @end deftypefn

function p = logicdf (x, mu = 0, s = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of X, MU, and S
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(s))
    [retval, x, mu, s] = common_size (x, mu, s);
    if (retval > 0)
      error (strcat (["logicdf: X, MU, and S must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, MU, and S being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (s))
    error ("logicdf: X, MU, and S must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (mu, "single") || isa (s, "single"));
    p = NaN (size (x), "single");
  else
    p = NaN (size (x));
  endif

  ## Compute logistic CDF
  k1 = (x == -Inf) & (s > 0);
  p(k1) = 0;

  k2 = (x == Inf) & (s > 0);
  p(k2) = 1;

  k = ! k1 & ! k2 & (s > 0);
  p(k) = 1 ./ (1 + exp (- (x(k) - mu(k)) ./ s(k)));

endfunction

%!demo
%! ## Plot various CDFs from the logistic distribution
%! x = -5:0.01:20;
%! p1 = logicdf (x, 5, 2);
%! p2 = logicdf (x, 9, 3);
%! p3 = logicdf (x, 9, 4);
%! p4 = logicdf (x, 6, 2);
%! p5 = logicdf (x, 2, 1);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-m")
%! grid on
%! legend ({"μ = 5, s = 2", "μ = 9, s = 3", "μ = 2, s = 4", ...
%!          "μ = 6, s = 2", "μ = 2, s = 1"}, "location", "southeast")
%! title ("Logistic CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x,y
%! x = [-Inf -log(3) 0 log(3) Inf];
%! y = [0, 1/4, 1/2, 3/4, 1];
%!assert (logicdf ([x, NaN]), [y, NaN], eps)
%!assert (logicdf (x, 0, [-2, -1, 0, 1, 2]), [nan(1, 3), 0.75, 1])

## Test class of input preserved
%!assert (logicdf (single ([x, NaN])), single ([y, NaN]), eps ("single"))

## Test input validation
%!error logicdf ()
%!error logicdf (1, 2, 3, 4)
%!error<logicdf: X, MU, and S must be of common size or scalars.> ...
%! logicdf (1, ones (2), ones (3))
%!error<logicdf: X, MU, and S must be of common size or scalars.> ...
%! logicdf (ones (2), 1, ones (3))
%!error<logicdf: X, MU, and S must be of common size or scalars.> ...
%! logicdf (ones (2), ones (3), 1)
%!error<logicdf: X, MU, and S must not be complex.> ...
%! logicdf (i, 2, 3)
%!error<logicdf: X, MU, and S must not be complex.> ...
%! logicdf (1, i, 3)
%!error<logicdf: X, MU, and S must not be complex.> ...
%! logicdf (1, 2, i)
