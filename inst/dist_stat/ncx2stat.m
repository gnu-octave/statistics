## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} ncx2stat (@var{df}, @var{delta})
##
## Compute statistics for the noncentral @math{χ^2} distribution.
##
## @code{[@var{m}, @var{v}] = ncx2stat (@var{df}, @var{delta})} returns the mean
## and variance of the noncentral chi-square distribution with @var{df} degrees
## of freedom and noncentrality parameter @var{delta}.
##
## The size of @var{m} and @var{v} is the common size of the input arguments.
## Scalar input arguments @var{df} and @var{delta} are regarded as constant
## matrices of the same size as the other input.
##
## @seealso{ncx2cdf, ncx2inv, ncx2pdf, ncx2rnd}
## @end deftypefn

function [m, v] = ncx2stat (df, delta)

  ## Check for valid input arguments
  if (nargin <  2)
    error ("ncx2stat: too few input arguments.");
  endif

  ## Check and fix size of input arguments
  [err, df, delta] = common_size (df, delta);
  if (err > 0)
    error ("ncx2stat: input size mismatch.");
  endif

  ## Initialize mean and variance
  if (isa (df, "single") || isa (delta, "single"))
    m = NaN (size (df), "single");
    v = m;
  else
    m = NaN (size (df));
    v = m;
  endif

  ## Compute mean and variance for valid parameter values.
  k = (df > 0 & delta >= 0);
  if (any (k(:)))
    m(k) = delta(k) + df(k);
    v(k) = 2 * (df(k) + 2 * (delta(k)));
  endif

endfunction

## Input validation tests
%!error<ncx2stat: too few input arguments.> p = ncx2stat ();
%!error<ncx2stat: too few input arguments.> p = ncx2stat (1);
%!error<ncx2stat: input size mismatch.> p = ncx2stat ([4, 3], [3, 4, 5]);

## Output validation tests
%!shared df, d1
%! df = [2, 0, -1, 1, 4];
%! d1 = [1, NaN, 3, -1, 2];
%!assert (ncx2stat (df, d1), [3, NaN, NaN, NaN, 6]);
%!assert (ncx2stat ([df(1:2), df(4:5)], 1), [3, NaN, 2, 5]);
%!assert (ncx2stat ([df(1:2), df(4:5)], 3), [5, NaN, 4, 7]);
%!assert (ncx2stat ([df(1:2), df(4:5)], 2), [4, NaN, 3, 6]);
%!assert (ncx2stat (2, [d1(1), d1(3:5)]), [3, 5, NaN, 4]);
%!assert (ncx2stat (0, [d1(1), d1(3:5)]), [NaN, NaN, NaN, NaN]);
%!assert (ncx2stat (1, [d1(1), d1(3:5)]), [2, 4, NaN, 3]);
%!assert (ncx2stat (4, [d1(1), d1(3:5)]), [5, 7, NaN, 6]);
