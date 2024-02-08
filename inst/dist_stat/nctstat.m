## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} nctstat (@var{df}, @var{mu})
##
## Compute statistics for the noncentral @math{t}-distribution.
##
## @code{[@var{m}, @var{v}] = nctstat (@var{df}, @var{mu})} returns the mean
## and variance of the noncentral @math{t}-distribution with @var{df} degrees
## of freedom and noncentrality parameter @var{mu}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## Further information about the noncentral @math{t}-distribution can be found
## at @url{https://en.wikipedia.org/wiki/Noncentral_t-distribution}
##
## @seealso{nctcdf, nctinv, nctpdf, nctrnd, tstat}
## @end deftypefn

function [m, v] = nctstat (df, mu)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("nctstat: function called with too few input arguments.");
  endif

  ## Check for DF and MU being numeric
  if (! (isnumeric (df) && isnumeric (mu)))
    error ("nctstat: DF and MU must be numeric.");
  endif

  ## Check for DF and MU being real
  if (iscomplex (df) || iscomplex (mu))
    error ("nctstat: DF and MU must not be complex.");
  endif

  ## Check for common size of DF and MU
  if (! isscalar (df) || ! isscalar (mu))
    [retval, df, mu] = common_size (df, mu);
    if (retval > 0)
      error ("nctstat: DF and MU must be of common size or scalars.");
    endif
  endif

  ## Initialize mean and variance
  if (isa (df, "single") || isa (mu, "single"))
    m = NaN (size (df), "single");
    v = m;
  else
    m = NaN (size (df));
    v = m;
  endif

  ## Compute mean and variance for valid parameter values.
  mk = df > 1;
  if (any (mk(:)))
    m(mk) = mu(mk) .* sqrt ((df(mk) / 2)) .* ...
            gamma ((df(mk) - 1) / 2) ./ gamma (df(mk) / 2);
  endif
  vk = df > 2;
  if (any (vk(:)))
    v(vk) = (df(vk) ./ (df(vk) - 2)) .* ...
            (1 + mu(vk) .^2) - 0.5 * (df(vk) .* mu(vk) .^ 2) .* ...
            exp (2 * (gammaln ((df(vk) - 1) / 2) - gammaln (df(vk) / 2)));
  endif

endfunction

## Input validation tests
%!error<nctstat: function called with too few input arguments.> nctstat ()
%!error<nctstat: function called with too few input arguments.> nctstat (1)
%!error<nctstat: DF and MU must be numeric.> nctstat ({}, 2)
%!error<nctstat: DF and MU must be numeric.> nctstat (1, "")
%!error<nctstat: DF and MU must not be complex.> nctstat (i, 2)
%!error<nctstat: DF and MU must not be complex.> nctstat (1, i)
%!error<nctstat: DF and MU must be of common size or scalars.> ...
%! nctstat (ones (3), ones (2))
%!error<nctstat: DF and MU must be of common size or scalars.> ...
%! nctstat (ones (2), ones (3))

## Output validation tests
%!shared df, mu
%! df = [2, 0, -1, 1, 4];
%! mu = [1, NaN, 3, -1, 2];
%!assert (nctstat (df, mu), [1.7725, NaN, NaN, NaN, 2.5066], 1e-4);
%!assert (nctstat ([df(1:2), df(4:5)], 1), [1.7725, NaN, NaN, 1.2533], 1e-4);
%!assert (nctstat ([df(1:2), df(4:5)], 3), [5.3174, NaN, NaN, 3.7599], 1e-4);
%!assert (nctstat ([df(1:2), df(4:5)], 2), [3.5449, NaN, NaN, 2.5066], 1e-4);
%!assert (nctstat (2, [mu(1), mu(3:5)]), [1.7725,5.3174,-1.7725,3.5449], 1e-4);
%!assert (nctstat (0, [mu(1), mu(3:5)]), [NaN, NaN, NaN, NaN]);
%!assert (nctstat (1, [mu(1), mu(3:5)]), [NaN, NaN, NaN, NaN]);
%!assert (nctstat (4, [mu(1), mu(3:5)]), [1.2533,3.7599,-1.2533,2.5066], 1e-4);
