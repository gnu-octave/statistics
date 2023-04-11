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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} nctstat (@var{df}, @var{delta})
##
## Compute statistics for the noncentral @math{t} distribution.
##
## @code{[@var{m}, @var{v}] = nctstat (@var{df}, @var{delta})} returns the mean
## and variance of the noncentral @math{t} distribution with @var{df} degrees
## of freedom and noncentrality parameter @var{delta}.
##
## The size of @var{m} and @var{v} is the common size of the input arguments.
## Scalar input arguments @var{df} and @var{delta} are regarded as constant
## matrices of the same size as the other input.
##
## @seealso{nctcdf, nctinv, nctpdf, nctrnd}
## @end deftypefn

function [m, v] = nctstat (df, delta)

  ## Check for valid input arguments
  if (nargin <  2)
    error ("nctstat: too few input arguments.");
  endif

  ## Check and fix size of input arguments
  [err, df, delta] = common_size (df, delta);
  if (err > 0)
    error ("nctstat: input size mismatch.");
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
  mk = df > 1;
  if (any (mk(:)))
    m(mk) = delta(mk) .* sqrt ((df(mk) / 2)) .* ...
            gamma ((df(mk) - 1) / 2) ./ gamma (df(mk) / 2);
  endif
  vk = df > 2;
  if (any (vk(:)))
    v(vk) = (df(vk) ./ (df(vk) - 2)) .* ...
            (1 + delta(vk) .^2) - 0.5 * (df(vk) .* delta(vk) .^ 2) .* ...
            exp (2 * (gammaln ((df(vk) - 1) / 2) - gammaln (df(vk) / 2)));
  endif

endfunction

## Input validation tests
%!error<nctstat: too few input arguments.> p = nctstat ();
%!error<nctstat: too few input arguments.> p = nctstat (1);
%!error<nctstat: input size mismatch.> p = nctstat ([4, 3], [3, 4, 5]);

## Output validation tests
%!shared df, d1
%! df = [2, 0, -1, 1, 4];
%! d1 = [1, NaN, 3, -1, 2];
%!assert (nctstat (df, d1), [1.7725, NaN, NaN, NaN, 2.5066], 1e-4);
%!assert (nctstat ([df(1:2), df(4:5)], 1), [1.7725, NaN, NaN, 1.2533], 1e-4);
%!assert (nctstat ([df(1:2), df(4:5)], 3), [5.3174, NaN, NaN, 3.7599], 1e-4);
%!assert (nctstat ([df(1:2), df(4:5)], 2), [3.5449, NaN, NaN, 2.5066], 1e-4);
%!assert (nctstat (2, [d1(1), d1(3:5)]), [1.7725,5.3174,-1.7725,3.5449], 1e-4);
%!assert (nctstat (0, [d1(1), d1(3:5)]), [NaN, NaN, NaN, NaN]);
%!assert (nctstat (1, [d1(1), d1(3:5)]), [NaN, NaN, NaN, NaN]);
%!assert (nctstat (4, [d1(1), d1(3:5)]), [1.2533,3.7599,-1.2533,2.5066], 1e-4);
