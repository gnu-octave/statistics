## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} fstat (@var{df1}, @var{df2})
##
## Compute statistics of the @math{F}-distribution.
##
## @code{[@var{m}, @var{v}] = fstat (@var{df1}, @var{df2})} returns the mean and
## variance of the @math{F}-distribution with @var{df1} and @var{df2} degrees
## of freedom.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the @math{F}-distribution can be found at
## @url{https://en.wikipedia.org/wiki/F-distribution}
##
## @seealso{fcdf, finv, fpdf, frnd}
## @end deftypefn

function [m, v] = fstat (df1, df2)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("fstat: function called with too few input arguments.");
  endif

  ## Check for DF1 and DF2 being numeric
  if (! (isnumeric (df1) && isnumeric (df2)))
    error ("fstat: DF1 and DF2 must be numeric.");
  endif

  ## Check for DF1 and DF2 being real
  if (iscomplex (df1) || iscomplex (df2))
    error ("fstat: DF1 and DF2 must not be complex.");
  endif

  ## Check for common size of DF1 and DF2
  if (! isscalar (df1) || ! isscalar (df2))
    [retval, df1, df2] = common_size (df1, df2);
    if (retval > 0)
      error ("fstat: DF1 and DF2 must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = df2 ./ (df2 - 2);
  v = (2 .* (df2 .^ 2) .* (df1 + df2 - 2)) ./ ...
      (df1 .* ((df2 - 2) .^ 2) .* (df2 - 4));

  ## Continue argument check
  k = find (! (df1 > 0) | ! (df1 < Inf) | ! (df2 > 2) | ! (df2 < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

  k = find (! (df2 > 4));
  if (any (k))
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<fstat: function called with too few input arguments.> fstat ()
%!error<fstat: function called with too few input arguments.> fstat (1)
%!error<fstat: DF1 and DF2 must be numeric.> fstat ({}, 2)
%!error<fstat: DF1 and DF2 must be numeric.> fstat (1, "")
%!error<fstat: DF1 and DF2 must not be complex.> fstat (i, 2)
%!error<fstat: DF1 and DF2 must not be complex.> fstat (1, i)
%!error<fstat: DF1 and DF2 must be of common size or scalars.> ...
%! fstat (ones (3), ones (2))
%!error<fstat: DF1 and DF2 must be of common size or scalars.> ...
%! fstat (ones (2), ones (3))

## Output validation tests
%!test
%! df1 = 1:6;
%! df2 = 5:10;
%! [m, v] = fstat (df1, df2);
%! expected_mn = [1.6667, 1.5000, 1.4000, 1.3333, 1.2857, 1.2500];
%! expected_v = [22.2222, 6.7500, 3.4844, 2.2222, 1.5869, 1.2153];
%! assert (m, expected_mn, 0.001);
%! assert (v, expected_v, 0.001);
%!test
%! df1 = 1:6;
%! [m, v] = fstat (df1, 5);
%! expected_mn = [1.6667, 1.6667, 1.6667, 1.6667, 1.6667, 1.6667];
%! expected_v = [22.2222, 13.8889, 11.1111, 9.7222, 8.8889, 8.3333];
%! assert (m, expected_mn, 0.001);
%! assert (v, expected_v, 0.001);
