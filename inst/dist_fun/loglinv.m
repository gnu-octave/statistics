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
## @deftypefn  {statistics} {@var{x} =} loglinv (@var{p})
## @deftypefnx {statistics} {@var{x} =} loglinv (@var{p}, @var{alpha})
## @deftypefnx {statistics} {@var{x} =} loglinv (@var{p}, @var{alpha}, @var{beta})
##
## Inverse of the log-logistic cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the log-logistic distribution with scale parameter @var{alpha}
## and shape parameter @var{beta}.  The size of @var{x} is the common size of
## @var{p}, @var{alpha}, and @var{beta}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Both parameters, @math{α} and @math{β}, must be positive reals and @var{p} is
## supported in the range @math{[0,1]}, otherwise @qcode{NaN} is returned.  By
## default, @qcode{@var{alpha} = 1} and @qcode{@var{beta} = 1}.
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
## @seealso{loglcdf, loglpdf, loglrnd, loglfit, logllike, loglstat}
## @end deftypefn

function x = loglinv (p, alpha = 1, beta = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of P, ALPHA, and BETA
  if (! isscalar (p) || ! isscalar (alpha) || ! isscalar(beta))
    [retval, p, alpha, beta] = common_size (p, alpha, beta);
    if (retval > 0)
      error (strcat (["loglinv: P, ALPHA, and BETA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, ALPHA, and BETA being reals
  if (iscomplex (p) || iscomplex (alpha) || iscomplex (beta))
    error ("loglinv: P, ALPHA, and BETA must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (p, "single") || isa (alpha, "single") || isa (beta, "single"));
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Check for valid points
  k = (p >= 0) & (p <= 1) & (alpha > 0) & (beta > 0);

  ## Compute the log-logistic iCDF
  x(k) = alpha(k) .* (p(k) ./ (1 - p(k))) .^ (1 ./ beta(k));

endfunction

%!demo
%! ## Plot various iCDFs from the log-logistic distribution
%! p = 0.001:0.001:0.999;
%! x1 = loglinv (p, 1, 0.5);
%! x2 = loglinv (p, 1, 1);
%! x3 = loglinv (p, 1, 2);
%! x4 = loglinv (p, 1, 4);
%! x5 = loglinv (p, 1, 8);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-m")
%! ylim ([0, 20])
%! grid on
%! legend ({"β = 0.5", "β = 1", "β = 2", "β = 4", "β = 8"}, ...
%!         "location", "northwest")
%! title ("Log-logistic iCDF")
%! xlabel ("probability")
%! ylabel ("x")
%! text (0.03, 12.5, "α = 1, values of β as shown in legend")

## Test output
%!shared p, out1, out2
%! p = [-1, 0, 0.2, 0.5, 0.8, 0.95, 1, 2];
%! out1 = [NaN, 0, 0.25, 1, 4, 19, Inf, NaN];
%! out2 = [NaN, 0, 0.0424732, 2.718282, 173.970037, 18644.695061, Inf, NaN];
%!assert (loglinv (p), out1, 1e-8)
%!assert (loglinv (p, exp (0), 1), out1, 1e-8)
%!assert (loglinv (p, exp (1), 1 / 3), out2, 1e-6)

## Test class of input preserved
%!assert (class (loglinv (single (1), 2, 3)), "single")
%!assert (class (loglinv (1, single (2), 3)), "single")
%!assert (class (loglinv (1, 2, single (3))), "single")

## Test input validation
%!error loglinv (1, 2, 3, 4)
%!error<loglinv: P, ALPHA, and BETA must be of common size or scalars.> ...
%! loglinv (1, ones (2), ones (3))
%!error<loglinv: P, ALPHA, and BETA must be of common size or scalars.> ...
%! loglinv (ones (2), 1, ones (3))
%!error<loglinv: P, ALPHA, and BETA must be of common size or scalars.> ...
%! loglinv (ones (2), ones (3), 1)
%!error<loglinv: P, ALPHA, and BETA must not be complex.> loglinv (i, 2, 3)
%!error<loglinv: P, ALPHA, and BETA must not be complex.> loglinv (1, i, 3)
%!error<loglinv: P, ALPHA, and BETA must not be complex.> loglinv (1, 2, i)
