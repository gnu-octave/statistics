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
## @deftypefn  {statistics} {@var{x} =} logistic_inv (@var{p})
## @deftypefnx {statistics} {@var{x} =} logistic_inv (@var{p}, @var{mu})
## @deftypefnx {statistics} {@var{x} =} logistic_inv (@var{p}, @var{mu}, @var{scale})
##
## Inverse of the logistic cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the logistic distribution with mean @var{mu} and scale
## parameter @var{scale}.  The size of @var{p} is the common size of @var{p},
## @var{mu}, and @var{scale}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{scale} = 1.  Both parameters must be
## reals and @var{beta} > 0.  For @var{beta} <= 0, NaN is returned.
##
## @seealso{logistic_cdf, logistic_pdf, logistic_rnd}
## @end deftypefn

function x = logistic_inv (p, mu = 0, scale = 1)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  ## Check for common size of P, MU, and SCALE
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar(scale))
    [retval, p, mu, scale] = common_size (p, mu, scale);
    if (retval > 0)
      error (strcat (["logistic_inv: P, MU, and SCALE must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, MU, and SCALE being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (scale))
    error ("logistic_inv: P, MU, and SCALE must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (p, "single") || isa (mu, "single") || isa (scale, "single"));
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  k = (p == 0) & (scale > 0);
  x(k) = -Inf;

  k = (p == 1) & (scale > 0);
  x(k) = Inf;

  k = (p > 0) & (p < 1) & (scale > 0);
  x(k) = mu(k) + scale(k) .* log (p(k) ./ (1 - p(k)));

endfunction


%!test
%! p = [0.01:0.01:0.99];
%! assert (logistic_inv (p), log (p ./ (1-p)), 25*eps);
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (logistic_inv (p), [NaN -Inf 0 Inf NaN])
%!assert (logistic_inv (p, 0, [-1, 0, 1, 2, 3]), [NaN NaN 0 Inf NaN])

## Test class of input preserved
%!assert (logistic_inv ([p, NaN]), [NaN -Inf 0 Inf NaN NaN])
%!assert (logistic_inv (single ([p, NaN])), single ([NaN -Inf 0 Inf NaN NaN]))

## Test input validation
%!error logistic_inv ()
%!error logistic_inv (1, 2, 3, 4)
%!error<logistic_inv: P, MU, and SCALE must be of common size or scalars.> ...
%! logistic_inv (1, ones (2), ones (3))
%!error<logistic_inv: P, MU, and SCALE must be of common size or scalars.> ...
%! logistic_inv (ones (2), 1, ones (3))
%!error<logistic_inv: P, MU, and SCALE must be of common size or scalars.> ...
%! logistic_inv (ones (2), ones (3), 1)
%!error<logistic_inv: P, MU, and SCALE must not be complex.> ...
%! logistic_inv (i, 2, 3)
%!error<logistic_inv: P, MU, and SCALE must not be complex.> ...
%! logistic_inv (1, i, 3)
%!error<logistic_inv: P, MU, and SCALE must not be complex.> ...
%! logistic_inv (1, 2, i)
