## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {statistics} {[@var{m}, @var{v}] =} stdrstat (@var{k}, @var{df})
##
## Compute statistics of the studentized range distribution.
##
## @code{[@var{m}, @var{v}] = stdrstat (@var{k}, @var{df})} returns the mean
## and variance of the studentized range distribution for @var{k} groups and
## @var{df} degrees of freedom.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the same
## size as the other input.  @var{k} must be an integer of at least 2 and
## @var{df} positive, @code{Inf} included; otherwise both are @code{NaN}.  The
## mean exists only for @var{df} greater than 1 and the variance only for
## @var{df} greater than 2; where either does not, it is @code{NaN}.  For
## @var{df} equal to @code{Inf} they are the mean and variance of the range of
## @var{k} standard normal variables.
##
## MATLAB has no public counterpart of this function.
##
## Further information about the studentized range distribution can be found
## at @url{https://en.wikipedia.org/wiki/Studentized_range_distribution}
##
## @seealso{stdrcdf, stdrinv, stdrpdf, stdrrnd}
## @end deftypefn

function [m, v] = stdrstat (k, df)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("stdrstat: function called with too few input arguments.");
  endif

  ## Check for K and DF being numeric
  if (! (isnumeric (k) && isnumeric (df)))
    error ("stdrstat: K and DF must be numeric.");
  endif

  ## Check for K and DF being real
  if (iscomplex (k) || iscomplex (df))
    error ("stdrstat: K and DF must not be complex.");
  endif

  ## Check for common size of K and DF
  if (! isscalar (k) || ! isscalar (df))
    [retval, k, df] = common_size (k, df);
    if (retval > 0)
      error ("stdrstat: K and DF must be of common size or scalars.");
    endif
  endif

  sz = size (k);
  k = double (k(:));
  df = double (df(:)) .* ones (numel (k), 1);
  m = NaN (numel (k), 1);
  v = NaN (numel (k), 1);
  valid = (k >= 2) & (k == fix (k)) & isfinite (k) & (df > 0);

  ## The moments of the range W of k standard normals, from its upper tail:
  ## E[W] = int U(w) dw and E[W^2] = 2 int w U(w) dw
  for kk = unique (k(valid & df > 1))'
    U = @(w) stdrcdf (w, kk, Inf, 'upper');
    EW = integral (U, 0, 40, 'RelTol', 1e-13, 'AbsTol', 0);
    EW2 = 2 * integral (@(w) w .* U (w), 0, 40, 'RelTol', 1e-13, 'AbsTol', 0);

    ## E[1/s] with s^2 a chi-square over df, from its asymptotic series once
    ## the direct gamma ratio would lose digits
    i = valid & (k == kk) & (df > 1);
    n = df(i);
    Es = ones (size (n));
    big = isfinite (n) & (n > 1e3);
    nb = n(big);
    Es(big) = 1 + 3 ./ (4 * nb) + 25 ./ (32 * nb .^ 2) + 105 ./ (128 * nb .^ 3);
    sml = (n <= 1e3);
    ns = n(sml);
    Es(sml) = sqrt (ns / 2) .* exp (gammaln ((ns - 1) / 2) - gammaln (ns / 2));
    m(i) = EW * Es;

    ## E[1/s^2] = df / (df - 2)
    j = i & (df > 2);
    Es2 = df(j) ./ (df(j) - 2);
    Es2(isinf (df(j))) = 1;
    v(j) = EW2 * Es2 - m(j) .^ 2;
  endfor

  m = reshape (m, sz);
  v = reshape (v, sz);

endfunction

## The range of two standard normals is sqrt (2) times a half-normal, and of
## three has mean 3 / sqrt (pi)
%!test
%! [m, v] = stdrstat (2, Inf);
%! assert_equal (m, 2 / sqrt (pi), -1e-13);
%! assert_equal (v, 2 - 4 / pi, -1e-12);
%!assert_equal (stdrstat (3, Inf), 3 / sqrt (pi), -1e-13)

## Two groups on DF degrees of freedom are sqrt (2) times the absolute value
## of a t variable
%!test
%! df = 5;
%! Et = 2 * sqrt (df) * exp (gammaln ((df + 1) / 2) - gammaln (df / 2)) ...
%!      / (sqrt (pi) * (df - 1));
%! [m, v] = stdrstat (2, df);
%! assert_equal (m, sqrt (2) * Et, -1e-12);
%! assert_equal (v, 2 * df / (df - 2) - 2 * Et ^ 2, -1e-11);

## The mean of the range of five standard normals, the control chart
## constant d2
%!assert_equal (stdrstat (5, Inf), 2.325929, 1e-6)

## Moments that do not exist, and invalid parameters
%!test
%! [m, v] = stdrstat (3, [1, 1.5, 2]);
%! assert_equal (isnan (m), [true, false, false]);
%! assert_equal (isnan (v), [true, true, true]);
%!test
%! [m, v] = stdrstat ([1, 2.5, NaN], 10);
%! assert_equal ([m, v], NaN (1, 6));

## Test input validation
%!error<stdrstat: function called with too few input arguments.> stdrstat ()
%!error<stdrstat: function called with too few input arguments.> stdrstat (3)
%!error<stdrstat: K and DF must be numeric.> stdrstat ({}, 10)
%!error<stdrstat: K and DF must be numeric.> stdrstat (3, '')
%!error<stdrstat: K and DF must not be complex.> stdrstat (i, 10)
%!error<stdrstat: K and DF must not be complex.> stdrstat (3, i)
%!error<stdrstat: K and DF must be of common size or scalars.> ...
%! stdrstat (ones (3), ones (2))
