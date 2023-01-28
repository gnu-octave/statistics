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
## @deftypefn  {statistics} @var{p} = laplace_cdf (@var{x})
## @deftypefnx {statistics} @var{p} = laplace_cdf (@var{x}, @var{mu})
## @deftypefnx {statistics} @var{p} = laplace_cdf (@var{x}, @var{mu}, @var{beta})
##
## Laplace cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the Laplace distribution with a location parameter
## @var{mu} and a scale parameter (i.e. "diversity") @var{beta}.  The size of
## @var{p} is the common size of @var{x}, @var{mu}, and @var{beta}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{beta} = 1.
##
## @seealso{laplace_inv, laplace_pdf, laplace_rnd}
## @end deftypefn

function p = laplace_cdf (x, mu = 0, beta = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of X, MU, and BETA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(beta))
    [retval, x, mu, beta] = ...
        common_size (x, mu, beta);
    if (retval > 0)
      error (strcat (["laplace_cdf: X, MU, and BETA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, MU, and BETA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (beta))
    error ("laplace_cdf: X, MU, and BETA must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (mu, "single") || isa (beta, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Compute Laplace CDF
  p = (1 + sign (x - mu) .* (1 - exp (- abs (x - mu) ./ beta))) ./ 2;

  ## Cast to appropriate class
  p = cast (p, is_class);

endfunction


%!shared x,y
%! x = [-Inf -log(2) 0 log(2) Inf];
%! y = [0, 1/4, 1/2, 3/4, 1];
%!assert (laplace_cdf ([x, NaN]), [y, NaN])

## Test class of input preserved
%!assert (laplace_cdf (single ([x, NaN])), single ([y, NaN]))

## Test input validation
%!error laplace_cdf ()
%!error laplace_cdf (1, 2, 3, 4)
%!error<laplace_cdf: X, MU, and BETA must be of common size or scalars.> ...
%! laplace_cdf (1, ones (2), ones (3))
%!error<laplace_cdf: X, MU, and BETA must be of common size or scalars.> ...
%! laplace_cdf (ones (2), 1, ones (3))
%!error<laplace_cdf: X, MU, and BETA must be of common size or scalars.> ...
%! laplace_cdf (ones (2), ones (3), 1)
%!error<laplace_cdf: X, MU, and BETA must not be complex.> laplace_cdf (i, 2, 3)
%!error<laplace_cdf: X, MU, and BETA must not be complex.> laplace_cdf (1, i, 3)
%!error<laplace_cdf: X, MU, and BETA must not be complex.> laplace_cdf (1, 2, i)
