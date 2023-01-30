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
## @deftypefn  {statistics} @var{x} = stdnormal_inv (@var{p})
##
## Inverse of the standard normal cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) at
## @var{p} of the standard normal distribution (mean = 0, standard
## deviation = 1).
##
## @seealso{norminv, stdnormal_cdf, stdnormal_pdf, stdnormal_rnd}
## @end deftypefn

function x = stdnormal_inv (p)

  if (nargin != 1)
    print_usage ();
  endif

  if (iscomplex (p))
    error ("stdnormal_inv: P must not be complex.");
  endif

  x = - sqrt (2) * erfcinv (2 * p);

endfunction


%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (stdnormal_inv (p), [NaN -Inf 0 Inf NaN])

## Test class of input preserved
%!assert (stdnormal_inv ([p, NaN]), [NaN -Inf 0 Inf NaN NaN])
%!assert (stdnormal_inv (single ([p, NaN])), single ([NaN -Inf 0 Inf NaN NaN]))

## Test input validation
%!error stdnormal_inv ()
%!error stdnormal_inv (1,2)
%!error stdnormal_inv (i)
