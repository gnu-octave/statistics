## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{x} =} raylinv (@var{p}, @var{sigma})
##
## Inverse of the Rayleigh cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Rayleigh distribution with parameter @var{sigma}.  The size
## of @var{p} is the common size of @var{x} and @var{sigma}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## Default value is @var{sigma} = 1.
##
## @subheading References
##
## @enumerate
## @item
## Wendy L. Martinez and Angel R. Martinez. @cite{Computational Statistics
## Handbook with MATLAB}. Appendix E, pages 547-557, Chapman & Hall/CRC,
## 2001.
##
## @item
## Athanasios Papoulis. @cite{Probability, Random Variables, and Stochastic
## Processes}. pages 104 and 148, McGraw-Hill, New York, second edition,
## 1984.
## @end enumerate
##
## @seealso{raylcdf, raylpdf, raylrnd, raylstat}
## @end deftypefn

function x = raylinv (p, sigma)

  ## Check arguments
  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif

  if (nargin < 1)
    sigma = 1;
  endif

  ## Check for common size of P and SIGMA
  if (! isscalar (p) || ! isscalar (sigma))
    [retval, p, sigma] = common_size (p, sigma);
    if (retval > 0)
      error ("raylinv: P and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X and SIGMA being reals
  if (iscomplex (p) || iscomplex (sigma))
    error ("raylcdf: P and SIGMA must not be complex.");
  endif

  ## Calculate Rayleigh iCDF
  x = sqrt (-2 .* log (1 - p) .* sigma .^ 2);

  ## Check for valid parameter and support
  k = find (p == 1);
  if (any (k))
    x(k) = Inf;
  endif

  k = find (! (p >= 0) | ! (p <= 1) | ! (sigma > 0));
  if (any (k))
    x(k) = NaN;
  endif

endfunction

%!test
%! p = 0:0.1:0.5;
%! sigma = 1:6;
%! x = raylinv (p, sigma);
%! expected_x = [0.0000, 0.9181, 2.0041, 3.3784, 5.0538, 7.0645];
%! assert (x, expected_x, 0.001);

%!test
%! p = 0:0.1:0.5;
%! x = raylinv (p, 0.5);
%! expected_x = [0.0000, 0.2295, 0.3340, 0.4223, 0.5054, 0.5887];
%! assert (x, expected_x, 0.001);

## Test input validation
%!error poissinv ()
%!error poissinv (1)
%!error poissinv (1,2,3)
%!error poissinv (ones (3), ones (2))
%!error poissinv (ones (2), ones (3))
%!error poissinv (i, 2)
%!error poissinv (2, i)
