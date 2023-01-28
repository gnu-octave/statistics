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
## @deftypefn  {statistics} @var{y} = laplace_pdf (@var{x})
## @deftypefnx {statistics} @var{y} = laplace_pdf (@var{x}, @var{mu})
## @deftypefnx {statistics} @var{y} = laplace_pdf (@var{x}, @var{mu}, @var{beta})
##
## Laplace probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Laplace distribution with a location parameter @var{mu} and
## a scale parameter (i.e. "diversity") @var{beta}.  The size of @var{y} is the
## common size of @var{x}, @var{mu}, and @var{beta}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{beta} = 1.
##
## @seealso{laplace_cdf, laplace_inv, laplace_rnd}
## @end deftypefn

function pdf = laplace_pdf (x, mu = 0, beta = 1)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  if (iscomplex (x))
    error ("laplace_pdf: X must not be complex");
  endif

  pdf = exp (- abs (x - mu) / beta) / (2 * beta);

endfunction


%!shared x,y
%! x = [-Inf -log(2) 0 log(2) Inf];
%! y = [0, 1/4, 1/2, 1/4, 0];
%!assert (laplace_pdf ([x, NaN]), [y, NaN])

## Test class of input preserved
%!assert (laplace_pdf (single ([x, NaN])), single ([y, NaN]))

## Test input validation
%!error laplace_pdf ()
%!error laplace_pdf (1,2,3,4)
%!error laplace_pdf (i)
