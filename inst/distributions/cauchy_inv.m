## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} cauchy_inv (@var{p})
## @deftypefnx {statistics} {@var{x} =} cauchy_inv (@var{p}, @var{location}, @var{scale})
##
## Inverse of the Cauchy cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Cauchy distribution with location parameter @var{location}
## and scale parameter @var{scale}.  The size of @var{x} is the common size of
## @var{p}, @var{location}, and @var{scale}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## Default values are @var{location} = 0, @var{scale} = 1.
##
## @seealso{cauchy_cdf, cauchy_pdf, cauchy_rnd}
## @end deftypefn

function x = cauchy_inv (p, location = 0, scale = 1)

  if (nargin != 1 && nargin != 3)
    print_usage ();
  endif

  if (! isscalar (p) || ! isscalar (location) || ! isscalar (scale))
    [retval, p, location, scale] = common_size (p, location, scale);
    if (retval > 0)
      error (strcat (["cauchy_inv: P, LOCATION, and SCALE must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  if (iscomplex (p) || iscomplex (location) || iscomplex (scale))
    error ("cauchy_inv: P, LOCATION, and SCALE must not be complex.");
  endif

  if (isa (p, "single") || isa (location, "single") || isa (scale, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ok = ! isinf (location) & (scale > 0) & (scale < Inf);

  k = (p == 0) & ok;
  x(k) = -Inf;

  k = (p == 1) & ok;
  x(k) = Inf;

  k = (p > 0) & (p < 1) & ok;
  if (isscalar (location) && isscalar (scale))
    x(k) = location - scale * cot (pi * p(k));
  else
    x(k) = location(k) - scale(k) .* cot (pi * p(k));
  endif

endfunction


%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (cauchy_inv (p, ones (1,5), 2*ones (1,5)), [NaN -Inf 1 Inf NaN], eps)
%!assert (cauchy_inv (p, 1, 2*ones (1,5)), [NaN -Inf 1 Inf NaN], eps)
%!assert (cauchy_inv (p, ones (1,5), 2), [NaN -Inf 1 Inf NaN], eps)
%!assert (cauchy_inv (p, [1 -Inf NaN Inf 1], 2), [NaN NaN NaN NaN NaN])
%!assert (cauchy_inv (p, 1, 2*[1 0 NaN Inf 1]), [NaN NaN NaN NaN NaN])
%!assert (cauchy_inv ([p(1:2) NaN p(4:5)], 1, 2), [NaN -Inf NaN Inf NaN])

## Test class of input preserved
%!assert (cauchy_inv ([p, NaN], 1, 2), [NaN -Inf 1 Inf NaN NaN], eps)
%!assert (cauchy_inv (single ([p, NaN]), 1, 2), ...
%! single ([NaN -Inf 1 Inf NaN NaN]), eps ("single"))
%!assert (cauchy_inv ([p, NaN], single (1), 2), ...
%! single ([NaN -Inf 1 Inf NaN NaN]), eps ("single"))
%!assert (cauchy_inv ([p, NaN], 1, single (2)), ...
%! single ([NaN -Inf 1 Inf NaN NaN]), eps ("single"))

## Test input validation
%!error cauchy_inv ()
%!error cauchy_inv (1,2)
%!error cauchy_inv (1,2,3,4)
%!error cauchy_inv (ones (3), ones (2), ones (2))
%!error cauchy_inv (ones (2), ones (3), ones (2))
%!error cauchy_inv (ones (2), ones (2), ones (3))
%!error cauchy_inv (i, 2, 2)
%!error cauchy_inv (2, i, 2)
%!error cauchy_inv (2, 2, i)
