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
## @deftypefn  {statistics} {@var{p} =} loglcdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} loglcdf (@var{x}, @var{alpha})
## @deftypefnx {statistics} {@var{p} =} loglcdf (@var{x}, @var{alpha}, @var{beta})
## @deftypefnx {statistics} {@var{p} =} loglcdf (@dots{}, @var{uflag})
##
## Log-logistic cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the log-logistic distribution with scale parameter
## @var{alpha} and shape parameter @var{beta}.  The size of @var{p} is the
## common size of @var{x}, @var{alpha}, and @var{beta}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## Both parameters, @math{α} and @math{β}, must be positive reals and @var{x} is
## supported in the range @math{[0,inf)}, otherwise @qcode{NaN} is returned.  By
## default, @qcode{@var{alpha} = 1} and @qcode{@var{beta} = 1}.
##
## @code{@var{p} = loglcdf (@var{x}, @var{alpha}, @var{beta}, "upper")} returns
## the upper tail probability of the log-logistic distribution with parameters
## @var{alpha}, and @var{beta} at the values in @var{x}.
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
## @seealso{loglinv, loglpdf, loglrnd, loglfit, logllike, loglstat}
## @end deftypefn

function p = loglcdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin > 4)
    error ("loglcdf: too many input arguments.");
  endif

  ## Check for 'upper' flag
  if (nargin > 1 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 1 && ischar (varargin{end}) && ...
          ! strcmpi (varargin{end}, "upper"))
    error ("loglcdf: invalid argument for upper tail.");
  elseif (nargin > 3)
    error ("loglcdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Get extra arguments (if they exist) or add defaults
  if (numel (varargin) > 0)
    alpha = varargin{1};
  else
    alpha = 1;
  endif
  if (numel (varargin) > 1)
    beta = varargin{2};
  else
    beta = 1;
  endif

  ## Check for invalid points
  alpha(alpha <= 0) = NaN;
  beta(beta <= 0) = NaN;
  x(x < 0) = NaN;

  ## Check for common size of X, ALPHA, and BETA
  if (! isscalar (x) || ! isscalar (alpha) || ! isscalar(beta))
    [retval, x, alpha, beta] = common_size (x, alpha, beta);
    if (retval > 0)
      error (strcat (["loglcdf: X, ALPHA, and BETA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, ALPHA, and BETA being reals
  if (iscomplex (x) || iscomplex (alpha) || iscomplex (beta))
    error ("loglcdf: X, ALPHA, and BETA must not be complex.");
  endif

  ## Compute log-logistic CDF
  z = (x ./ alpha) .^ -beta;
  if (uflag)
    p = 1 - (1 ./ (1 + z));
  else
    p = 1 ./ (1 + z);
  endif

  ## Check for "single" class
  if (isa (x, "single") || isa (alpha, "single") || isa (beta, "single"));
    p = cast (p, "single");
  endif

endfunction

%!demo
%! ## Plot various CDFs from the log-logistic distribution
%! x = 0:0.001:2;
%! p1 = loglcdf (x, 1, 0.5);
%! p2 = loglcdf (x, 1, 1);
%! p3 = loglcdf (x, 1, 2);
%! p4 = loglcdf (x, 1, 4);
%! p5 = loglcdf (x, 1, 8);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-m")
%! legend ({"β = 0.5", "β = 1", "β = 2", "β = 4", "β = 8"}, ...
%!         "location", "northwest")
%! grid on
%! title ("Log-logistic CDF")
%! xlabel ("values in x")
%! ylabel ("probability")
%! text (0.05, 0.64, "α = 1, values of β as shown in legend")

## Test output
%!shared out1, out2
%! out1 = [0, 0.5, 0.66666667, 0.75, 0.8, 0.83333333];
%! out2 = [0, 0.4174, 0.4745, 0.5082, 0.5321, 0.5506];
%!assert (loglcdf ([0:5]), out1, 1e-8)
%!assert (loglcdf ([0:5], "upper"), 1 - out1, 1e-8)
%!assert (loglcdf ([0:5], exp (0), 1), out1, 1e-8)
%!assert (loglcdf ([0:5], exp (0), 1, "upper"), 1 - out1, 1e-8)
%!assert (loglcdf ([0:5], exp (1), 1 / 3), out2, 1e-4)
%!assert (loglcdf ([0:5], exp (1), 1 / 3, "upper"), 1 - out2, 1e-4)

## Test class of input preserved
%!assert (class (loglcdf (single (1), 2, 3)), "single")
%!assert (class (loglcdf (1, single (2), 3)), "single")
%!assert (class (loglcdf (1, 2, single (3))), "single")

## Test input validation
%!error<loglcdf: too many input arguments.> loglcdf (1, 2, 3, 4, 5)
%!error<loglcdf: invalid argument for upper tail.> ...
%! loglcdf (1, 2, 3, 4)
%!error<loglcdf: invalid argument for upper tail.> ...
%! loglcdf (1, 2, 3, "uper")
%!error<loglcdf: X, ALPHA, and BETA must be of common size or scalars.> ...
%! loglcdf (1, ones (2), ones (3))
%!error<loglcdf: X, ALPHA, and BETA must be of common size or scalars.> ...
%! loglcdf (1, ones (2), ones (3), "upper")
%!error<loglcdf: X, ALPHA, and BETA must be of common size or scalars.> ...
%! loglcdf (ones (2), 1, ones (3))
%!error<loglcdf: X, ALPHA, and BETA must be of common size or scalars.> ...
%! loglcdf (ones (2), 1, ones (3), "upper")
%!error<loglcdf: X, ALPHA, and BETA must be of common size or scalars.> ...
%! loglcdf (ones (2), ones (3), 1)
%!error<loglcdf: X, ALPHA, and BETA must be of common size or scalars.> ...
%! loglcdf (ones (2), ones (3), 1, "upper")
%!error<loglcdf: X, ALPHA, and BETA must not be complex.> ...
%! loglcdf (i, 2, 3)
%!error<loglcdf: X, ALPHA, and BETA must not be complex.> ...
%! loglcdf (i, 2, 3, "upper")
%!error<loglcdf: X, ALPHA, and BETA must not be complex.> ...
%! loglcdf (1, i, 3)
%!error<loglcdf: X, ALPHA, and BETA must not be complex.> ...
%! loglcdf (1, i, 3, "upper")
%!error<loglcdf: X, ALPHA, and BETA must not be complex.> ...
%! loglcdf (1, 2, i)
%!error<loglcdf: X, ALPHA, and BETA must not be complex.> ...
%! loglcdf (1, 2, i, "upper")
