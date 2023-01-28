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
## @deftypefn  {statistics} @var{x} = laplace_inv (@var{p})
## @deftypefnx {statistics} @var{x} = laplace_inv (@var{p}, @var{mu})
## @deftypefnx {statistics} @var{x} = laplace_inv (@var{p}, @var{mu}, @var{beta})
##
## Inverse of the Laplace cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Laplace distribution with a location parameter @var{mu} and
## a scale parameter (i.e. "diversity") @var{beta}.  The size of @var{x} is the
## common size of @var{p}, @var{mu}, and @var{beta}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{beta} = 1.
##
## @seealso{laplace_inv, laplace_pdf, laplace_rnd}
## @end deftypefn

function x = laplace_inv (p, mu = 0, beta = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of P, MU, and BETA
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar(beta))
    [retval, p, mu, beta] = common_size (p, mu, beta);
    if (retval > 0)
      error (strcat (["laplace_inv: P, MU, and BETA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, MU, and BETA being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (beta))
    error ("laplace_inv: P, MU, and BETA must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (p, "single") || isa (mu, "single") || isa (beta, "single"));
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Compute Laplace iCDF
  k = (p >= 0) & (p <= 1);
  x(k) = mu(k) + beta(k) .* ((p(k) < 1/2) .* log (2 .* p(k)) - ...
                             (p(k) > 1/2) .* log (2 .* (1 - p(k))));

endfunction


%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (laplace_inv (p), [NaN -Inf 0 Inf NaN])

## Test class of input preserved
%!assert (laplace_inv ([p, NaN]), [NaN -Inf 0 Inf NaN NaN])
%!assert (laplace_inv (single ([p, NaN])), single ([NaN -Inf 0 Inf NaN NaN]))

## Test input validation
%!error laplace_inv ()
%!error laplace_inv (1, 2, 3, 4)
%!error<laplace_inv: P, MU, and BETA must be of common size or scalars.> ...
%! laplace_inv (1, ones (2), ones (3))
%!error<laplace_inv: P, MU, and BETA must be of common size or scalars.> ...
%! laplace_inv (ones (2), 1, ones (3))
%!error<laplace_inv: P, MU, and BETA must be of common size or scalars.> ...
%! laplace_inv (ones (2), ones (3), 1)
%!error<laplace_inv: P, MU, and BETA must not be complex.> laplace_inv (i, 2, 3)
%!error<laplace_inv: P, MU, and BETA must not be complex.> laplace_inv (1, i, 3)
%!error<laplace_inv: P, MU, and BETA must not be complex.> laplace_inv (1, 2, i)
