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
## @deftypefn  {statistics} {@var{y} =} loglpdf (@var{x})
## @deftypefnx {statistics} {@var{y} =} loglpdf (@var{x}, @var{alpha})
## @deftypefnx {statistics} {@var{y} =} loglpdf (@var{x}, @var{alpha}, @var{beta})
##
## Log-logistic probability density function (PDF).
##
## For each element of @var{x}, compute the PDF at @var{x} of the log-logistic
## distribution with with scale parameter @var{alpha} and shape parameter
## @var{beta}.  The size of @var{y} is the common size of @var{x}, @var{alpha},
## and @var{beta}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## Both parameters, @math{α} and @math{β}, must be positive reals, otherwise
## @qcode{NaN} is returned.   @var{x} is supported in the range @math{[0,Inf)},
## otherwise @qcode{0} is returned.  By default, @qcode{@var{alpha} = 1} and
## @qcode{@var{beta} = 1}.
##
## Further information about the log-logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-logistic_distribution}
##
## MATLAB compatibility: MATLAB uses an alternative parameterization given by
## the pair @math{μ, s}, i.e. @var{mu} and @var{scale}, in analogy with the
## logistic distribution.  Their relation to the @var{alpha} and @var{beta}
## parameters is given below:
##
## @itemize
## @item @qcode{@var{alpha} = exp (@var{mu})}
## @item @qcode{@var{beta} = 1 / @var{scale}}
## @end itemize
##
## @seealso{loglcdf, loglinv, loglrnd, loglfit, logllike, loglstat}
## @end deftypefn

function y = loglpdf (x, alpha = 1, beta = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of X, ALPHA, and BETA
  if (! isscalar (x) || ! isscalar (alpha) || ! isscalar(beta))
    [retval, x, alpha, beta] = common_size (x, alpha, beta);
    if (retval > 0)
      error (strcat (["loglpdf: X, ALPHA, and BETA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, ALPHA, and BETA being reals
  if (iscomplex (x) || iscomplex (alpha) || iscomplex (beta))
    error ("loglpdf: X, ALPHA, and BETA must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (alpha, "single") || isa (beta, "single"));
    y = NaN (size (x), "single");
  else
    y = NaN (size (x));
  endif

  ## Compute log-logistic PDF
  k1 = ((x == Inf) | (x <= 0)) & (alpha > 0) & (beta > 0);
  y(k1) = 0;

  k = (! k1) & (alpha > 0) & (beta > 0);
  y(k) = ((beta(k) ./ alpha(k)) .* (x(k) ./ alpha(k)) .^ (beta(k) -1)) ./ ...
         ((1 + (x(k) ./ alpha(k)) .^ beta(k)) .^ 2);

endfunction

%!demo
%! ## Plot various PDFs from the log-logistic distribution
%! x = 0:0.001:2;
%! y1 = loglpdf (x, 1, 0.5);
%! y2 = loglpdf (x, 1, 1);
%! y3 = loglpdf (x, 1, 2);
%! y4 = loglpdf (x, 1, 4);
%! y5 = loglpdf (x, 1, 8);
%! plot (x, y1, "-b", x, y2,"-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
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
%! out1 = [0, 0, 0.2500, 0.1111, 0.0625, 0.0400, 0.0278, 0];
%! out2 = [0, 0, 0.0811, 0.0416, 0.0278, 0.0207, 0.0165, 0];
%!assert (loglpdf ([-1:5,Inf]), out1, 1e-4)
%!assert (loglpdf ([-1:5,Inf], exp (0), 1), out1, 1e-4)
%!assert (loglpdf ([-1:5,Inf], exp (1), 1 / 3), out2, 1e-4)

## Test class of input preserved
%!assert (class (loglpdf (single (1), 2, 3)), "single")
%!assert (class (loglpdf (1, single (2), 3)), "single")
%!assert (class (loglpdf (1, 2, single (3))), "single")

## Test input validation
%!error loglpdf (1, 2, 3, 4)
%!error<loglpdf: X, ALPHA, and BETA must be of common size or scalars.> ...
%! loglpdf (1, ones (2), ones (3))
%!error<loglpdf: X, ALPHA, and BETA must be of common size or scalars.> ...
%! loglpdf (ones (2), 1, ones (3))
%!error<loglpdf: X, ALPHA, and BETA must be of common size or scalars.> ...
%! loglpdf (ones (2), ones (3), 1)
%!error<loglpdf: X, ALPHA, and BETA must not be complex.> loglpdf (i, 2, 3)
%!error<loglpdf: X, ALPHA, and BETA must not be complex.> loglpdf (1, i, 3)
%!error<loglpdf: X, ALPHA, and BETA must not be complex.> loglpdf (1, 2, i)
